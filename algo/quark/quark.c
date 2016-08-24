#include "cpuminer-config.h"
#include "miner.h"
#include "algo-gate-api.h"

#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "algo/blake/sph_blake.h"
#include "algo/bmw/sph_bmw.h"
#include "algo/groestl/sph_groestl.h"
#include "algo/jh/sph_jh.h"
#include "algo/keccak/sph_keccak.h"
#include "algo/skein/sph_skein.h"

#include "algo/blake/sse2/blake.c"
#include "algo/bmw/sse2/bmw.c"
#include "algo/keccak/sse2/keccak.c"
#include "algo/skein/sse2/skein.c"
#include "algo/jh/sse2/jh_sse2_opt64.h"

#ifdef NO_AES_NI
  #include "algo/groestl/sse2/grso.h"
  #include "algo/groestl/sse2/grso-macro.c"
#else
 #include "algo/groestl/aes_ni/hash-groestl.h"
#endif

/*define data alignment for different C compilers*/
#if defined(__GNUC__)
      #define DATA_ALIGN16(x) x __attribute__ ((aligned(16)))
      #define DATA_ALIGNXY(x,y) x __attribute__ ((aligned(y)))

#else
      #define DATA_ALIGN16(x) __declspec(align(16)) x
      #define DATA_ALIGNXY(x,y) __declspec(align(y)) x
#endif

#ifndef NO_AES_NI
hashState_groestl quark_groestl_ctx;
#endif

void init_quark_ctx()
{
#ifndef NO_AES_NI
 init_groestl( &quark_groestl_ctx );
#endif
}

inline static void quarkhash(void *state, const void *input)
{
#ifdef NO_AES_NI
  grsoState sts_grs;
#else
  hashState_groestl ctx;
  memcpy(&ctx, &quark_groestl_ctx, sizeof(quark_groestl_ctx));
#endif

    /* shared  temp space */
    /* hash is really just 64bytes but it used to hold both hash and final round constants passed 64 */

    unsigned char hashbuf[128];
    size_t hashptr;
    sph_u64 hashctA;
    sph_u64 hashctB;

    int i;

    unsigned char hash[128];

    // Blake
    DECL_BLK;
    BLK_I;
    BLK_W;
    for(i=0; i<9; i++)
    {
    /* blake is split between 64byte hashes and the 80byte initial block */
    //DECL_BLK;
      switch (i+(16*((hash[0] & (uint32_t)(8)) == (uint32_t)(0))))
      {
        // Blake
        case 5 :
            BLK_I;
            BLK_U;
        case 0:
        case 16: 
            BLK_C;
            break;
        case 1:
        case 17:
        case 21:

            // BMW
            do
            { 
              DECL_BMW;
              BMW_I;
              BMW_U;
              /* bmw compress uses some defines */
              /* i havent gotten around to rewriting these */
              #define M(x)    sph_dec64le_aligned(data + 8 * (x))
              #define H(x)    (h[x])
              #define dH(x)   (dh[x])
              BMW_C;
              #undef M
              #undef H
              #undef dH
            } while(0); continue;;

        case 2:
            // dos this entry point represent a second groestl round?

        case 3:
        case 19:
          // Groestl 
          do
          {

#ifdef NO_AES_NI
           GRS_I;
           GRS_U;
           GRS_C;
#else
           reinit_groestl( &ctx );
           update_groestl(&ctx, (char*)hash,512);
           final_groestl(&ctx, (char*)hash);
#endif

          } while(0); continue;

        case 4:
        case 20:
        case 24:
            // JH
            do
            {
              DECL_JH;
              JH_H;
            } while(0); continue;

        case 6:
        case 22:
        case 8:
            // Keccak
            do
            {
              DECL_KEC;
              KEC_I;
              KEC_U;
              KEC_C;
            } while(0); continue;

        case 18:
        case 7:
        case 23:
            // Skein
            do
            {
              DECL_SKN;
              SKN_I;
              SKN_U;
              SKN_C; /* is a magintue faster than others, done */
            } while(0); continue;
 
       default:
            /* bad things happend, i counted to potato */
            abort();
    }
    /* only blake shouuld get here without continue */
    /* blake finishs from top split */
    //BLK_C;
 }
 

//    asm volatile ("emms");
  memcpy(state, hash, 32);
}

void quarkhash_alt(void *state, const void *input)
{
        sph_blake512_context    ctx_blake1,
                                ctx_blake2;
        sph_bmw512_context      ctx_bmw1,
                                ctx_bmw2;
        sph_groestl512_context  ctx_groestl1,
                                ctx_groestl2;
        sph_skein512_context    ctx_skein1,
                                ctx_skein2;
        sph_jh512_context       ctx_jh1,
                                ctx_jh2;
        sph_keccak512_context   ctx_keccak1,
                                ctx_keccak2;

        sph_blake512_init(&ctx_blake1);
        sph_bmw512_init(&ctx_bmw1);
        sph_groestl512_init(&ctx_groestl1);
        sph_skein512_init(&ctx_skein1);
        sph_groestl512_init(&ctx_groestl2);
        sph_jh512_init(&ctx_jh1);
        sph_blake512_init(&ctx_blake2);
        sph_bmw512_init(&ctx_bmw2);
        sph_keccak512_init(&ctx_keccak1);
        sph_skein512_init(&ctx_skein2);
        sph_keccak512_init(&ctx_keccak2);
        sph_jh512_init(&ctx_jh2);

        uint32_t _ALIGN(128) hash[16];
        uint32_t mask = 8;


        sph_blake512 (&ctx_blake1, input, 80);
        sph_blake512_close (&ctx_blake1, hash); //0

        sph_bmw512 (&ctx_bmw1, hash, 64);
        sph_bmw512_close(&ctx_bmw1, hash); //1

        if (hash[0] & mask) {
                sph_groestl512 (&ctx_groestl1, hash, 64);
                sph_groestl512_close(&ctx_groestl1, hash); //2
        } else {
                sph_skein512 (&ctx_skein1, hash, 64);
                sph_skein512_close(&ctx_skein1, hash); //2
        }

        sph_groestl512 (&ctx_groestl2, hash, 64);
        sph_groestl512_close(&ctx_groestl2, hash); //3

        sph_jh512 (&ctx_jh1, hash, 64);
        sph_jh512_close(&ctx_jh1, hash); //4

        if (hash[0] & mask) {
                sph_blake512 (&ctx_blake2, hash, 64);
                sph_blake512_close(&ctx_blake2, hash); //5
        } else {
                sph_bmw512 (&ctx_bmw2, hash, 64);
                sph_bmw512_close(&ctx_bmw2, hash); //5
        }

        sph_keccak512 (&ctx_keccak1, hash, 64);
        sph_keccak512_close(&ctx_keccak1, hash); //6

        sph_skein512 (&ctx_skein2, hash, 64);
        sph_skein512_close(&ctx_skein2, hash); //7

        if (hash[0] & mask) {
                sph_keccak512 (&ctx_keccak2, hash, 64);
                sph_keccak512_close(&ctx_keccak2, hash); //8
        } else {
                sph_jh512 (&ctx_jh2, hash, 64);
                sph_jh512_close(&ctx_jh2, hash); //8
        }

        memcpy(state, hash, 32);
}

int scanhash_quark( int thr_id, struct work *work, uint32_t max_nonce,
                    uint64_t *hashes_done)
{
        uint32_t endiandata[20] __attribute__((aligned(64)));
        uint32_t hash64[8] __attribute__((aligned(32)));
        uint32_t *pdata = work->data;
        uint32_t *ptarget = work->target;
	uint32_t n = pdata[19] - 1;
	const uint32_t first_nonce = pdata[19];

        swab32_array( endiandata, pdata, 20 );

	do {
		pdata[19] = ++n;
		be32enc(&endiandata[19], n); 
		quarkhash(hash64, &endiandata);
                if ((hash64[7]&0xFFFFFF00)==0)
                {
                  if (fulltest(hash64, ptarget)) 
                  {
                    *hashes_done = n - first_nonce + 1;
		    return true;
                  }
               }
	} while (n < max_nonce && !work_restart[thr_id].restart);
	
	*hashes_done = n - first_nonce + 1;
	pdata[19] = n;
	return 0;
}

bool register_quark_algo( algo_gate_t* gate )
{
  init_quark_ctx();
  gate->aes_ni_optimized = true;
  gate->optimizations = SSE2_OPT | AES_OPT | AVX_OPT;
  gate->scanhash         = (void*)&scanhash_quark;
  gate->hash             = (void*)&quarkhash;
  gate->hash_alt         = (void*)&quarkhash_alt;
  return true;
};

