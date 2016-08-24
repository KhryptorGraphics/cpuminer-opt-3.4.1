#include "miner.h"
#include "algo-gate-api.h"
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>

#include "sha3/sph_sha2.h"
#include "ripemd/sph_ripemd.h"

//#define DEBUG_ALGO

#define LBRY_NTIME_INDEX 25
#define LBRY_NBITS_INDEX 26
#define LBRY_NONCE_INDEX 27
#define LBRY_WORK_DATA_SIZE 192
#define LBRY_WORK_CMP_SIZE 76  // same as default

/* Move init out of loop, so init once externally, and then use one single memcpy with that bigger memory block */
typedef struct {
	sph_sha256_context	sha256;
	sph_sha512_context	sha512;
	sph_ripemd160_context	ripemd;
} lbryhash_context_holder;

/* no need to copy, because close reinit the context */
static  lbryhash_context_holder ctx;

void init_lbry_contexts(void *dummy)
{
	sph_sha256_init(&ctx.sha256);
	sph_sha512_init(&ctx.sha512);
	sph_ripemd160_init(&ctx.ripemd);
}

void lbry_hash(void* output, const void* input)
{
        sph_sha256_context      ctx_sha256;
        sph_sha512_context      ctx_sha512;
        sph_ripemd160_context   ctx_ripemd;
	uint32_t _ALIGN(64) hashA[16];
        uint32_t _ALIGN(64) hashB[16];
        uint32_t _ALIGN(64) hashC[16];

        sph_sha256_init(&ctx_sha256);
	sph_sha256 (&ctx_sha256, input, 112);
	sph_sha256_close(&ctx_sha256, hashA);

        sph_sha256_init(&ctx_sha256);
	sph_sha256 (&ctx_sha256, hashA, 32);
	sph_sha256_close(&ctx_sha256, hashA);

        sph_sha512_init(&ctx_sha512);
	sph_sha512 (&ctx_sha512, hashA, 32);
	sph_sha512_close(&ctx_sha512, hashA);

        sph_ripemd160_init(&ctx_ripemd);
	sph_ripemd160 (&ctx_ripemd, hashA, 32);
	sph_ripemd160_close(&ctx_ripemd, hashB);

        sph_ripemd160_init(&ctx_ripemd);
	sph_ripemd160 (&ctx_ripemd, hashA+8, 32);
	sph_ripemd160_close(&ctx_ripemd, hashC);

        sph_sha256_init(&ctx_sha256);
	sph_sha256 (&ctx_sha256, hashB, 20);
	sph_sha256 (&ctx_sha256, hashC, 20);
	sph_sha256_close(&ctx_sha256, hashA);

        sph_sha256_init(&ctx_sha256);
	sph_sha256 (&ctx_sha256, hashA, 32);
	sph_sha256_close(&ctx_sha256, hashA);

	memcpy(output, hashA, 32);
}

int scanhash_lbry( int thr_id, struct work *work, uint32_t max_nonce,
                   uint64_t *hashes_done)
{
  uint32_t *pdata = work->data;
  uint32_t *ptarget = work->target;
	uint32_t n = pdata[27] - 1;
	const uint32_t first_nonce = pdata[27];
	const uint32_t Htarg = ptarget[7];

	uint32_t hash64[8] __attribute__((aligned(32)));
	uint32_t endiandata[32];

	uint64_t htmax[] = {
		0,
		0xF,
		0xFF,
		0xFFF,
		0xFFFF,
		0x10000000
	};
	uint32_t masks[] = {
		0xFFFFFFFF,
		0xFFFFFFF0,
		0xFFFFFF00,
		0xFFFFF000,
		0xFFFF0000,
		0
	};

	// we need bigendian data...
        swab32_array( endiandata, pdata, 32 );

#ifdef DEBUG_ALGO
	printf("[%d] Htarg=%X\n", thr_id, Htarg);
#endif
	for (int m=0; m < sizeof(masks); m++) {
		if (Htarg <= htmax[m]) {
			uint32_t mask = masks[m];
			do {
				pdata[27] = ++n;
				be32enc(&endiandata[27], n);
				lbry_hash(hash64, &endiandata);
#ifndef DEBUG_ALGO
				if ((!(hash64[7] & mask)) && fulltest(hash64, ptarget)) {
					*hashes_done = n - first_nonce + 1;
					return true;
				}
#else
				if (!(n % 0x1000) && !thr_id) printf(".");
				if (!(hash64[7] & mask)) {
					printf("[%d]",thr_id);
					if (fulltest(hash64, ptarget)) {
						*hashes_done = n - first_nonce + 1;
						return true;
					}
				}
#endif
			} while (n < max_nonce && !work_restart[thr_id].restart);
			// see blake.c if else to understand the loop on htmax => mask
			break;
		}
	}

	*hashes_done = n - first_nonce + 1;
	pdata[27] = n;
	return 0;
}

void lbry_calc_network_diff(struct work *work)
{
        // sample for diff 43.281 : 1c05ea29
        // todo: endian reversed on longpoll could be zr5 specific...
//        uint32_t nbits = have_longpoll ? work->data[18] : swab32(work->data[18]);

   uint32_t nbits = swab32( work->data[ LBRY_NBITS_INDEX ] );
   uint32_t bits = (nbits & 0xffffff);
   int16_t shift = (swab32(nbits) & 0xff); // 0x1c = 28
   uint64_t diffone = 0x0000FFFF00000000ull;
   double d = (double)0x0000ffff / (double)bits;

   for (int m=shift; m < 29; m++) d *= 256.0;
   for (int m=29; m < shift; m++) d /= 256.0;
   if (opt_debug_diff)
      applog(LOG_DEBUG, "net diff: %f -> shift %u, bits %08x", d, shift, bits);

   net_diff = d;
}

// std_le should work but it doesn't
void lbry_le_build_stratum_request( char *req, struct work *work,
                                      struct stratum_ctx *sctx )
{
   unsigned char *xnonce2str;
   uint32_t ntime, nonce;
   char ntimestr[9], noncestr[9];

   le32enc( &ntime, work->data[ LBRY_NTIME_INDEX ] );
   le32enc( &nonce, work->data[ LBRY_NONCE_INDEX ] );
   bin2hex( ntimestr, (char*)(&ntime), sizeof(uint32_t) );
   bin2hex( noncestr, (char*)(&nonce), sizeof(uint32_t) );
   xnonce2str = abin2hex( work->xnonce2, work->xnonce2_len);
   snprintf( req, JSON_BUF_LEN,
        "{\"method\": \"mining.submit\", \"params\": [\"%s\", \"%s\", \"%s\", \"%s\", \"%s\"], \"id\":4}",
         rpc_user, work->job_id, xnonce2str, ntimestr, noncestr );
   free(xnonce2str);
}

void lbry_build_extraheader( struct work* work, struct stratum_ctx* sctx )
{
   for ( int i = 0; i < 8; i++ )
        work->data[17 + i] = ((uint32_t*)sctx->job.claim)[i];
   work->data[ LBRY_NTIME_INDEX ] = le32dec(sctx->job.ntime);
   work->data[ LBRY_NBITS_INDEX ] = le32dec(sctx->job.nbits);
   work->data[28] = 0x80000000;
}

void lbry_set_target( struct work* work, double job_diff )
{
 work_set_target( work, job_diff / (256.0 * opt_diff_factor) );
}

int64_t lbry_get_max64() { return 0x1ffffLL; }

bool register_lbry_algo( algo_gate_t* gate )
{
  gate->scanhash              = (void*)&scanhash_lbry;
  gate->hash                  = (void*)&lbry_hash;
  gate->hash_alt              = (void*)&lbry_hash;
  gate->calc_network_diff     = (void*)&lbry_calc_network_diff;
  gate->get_max64             = (void*)&lbry_get_max64;
  gate->build_stratum_request = (void*)&lbry_le_build_stratum_request;
  gate->build_extraheader     = (void*)&lbry_build_extraheader;
  gate->set_target            = (void*)&lbry_set_target;
  gate->ntime_index           = LBRY_NTIME_INDEX;
  gate->nbits_index           = LBRY_NBITS_INDEX;
  gate->nonce_index           = LBRY_NONCE_INDEX;
  gate->work_data_size        = LBRY_WORK_DATA_SIZE;
  return true;
}

