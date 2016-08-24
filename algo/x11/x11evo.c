#include "cpuminer-config.h"
#include "miner.h"
#include "algo-gate-api.h"

#include <string.h>
#include <stdint.h>
#include <compat/portable_endian.h>

#include "algo/blake/sph_blake.h"
#include "algo/bmw/sph_bmw.h"
#include "algo/groestl/sph_groestl.h"
#include "algo/jh/sph_jh.h"
#include "algo/keccak/sph_keccak.h"
#include "algo/skein/sph_skein.h"
#include "algo/luffa/sph_luffa.h"
#include "algo/cubehash/sph_cubehash.h"
#include "algo/shavite/sph_shavite.h"
#include "algo/simd/sph_simd.h"
#include "algo/echo/sph_echo.h"

#ifdef NO_AES_NI
//  #include "algo/groestl/sse2/grso.h"
//  #include "algo/groestl/sse2/grso-macro.c"
#else
  #include "algo/groestl/aes_ni/hash-groestl.h"
  #include "algo/echo/aes_ni/hash_api.h"
#endif

#include "algo/luffa/sse2/luffa_for_sse2.h"
#include "algo/cubehash/sse2/cubehash_sse2.h"
#include "algo/simd/sse2/nist.h"

#define INITIAL_DATE 1462060800
#define HASH_FUNC_COUNT 11

typedef struct {
#ifdef NO_AES_NI
    sph_groestl512_context  groestl;
    sph_echo512_context     echo;
#else
    hashState_echo          echo;
    hashState_groestl       groestl;
#endif
    hashState_luffa         luffa;
    cubehashParam           cube;
    hashState_sd            simd;
    sph_blake512_context    blake;
    sph_bmw512_context      bmw;
    sph_skein512_context    skein;
    sph_jh512_context       jh;
    sph_keccak512_context   keccak;
    sph_shavite512_context  shavite;
} x11evo_ctx_holder;

static x11evo_ctx_holder x11evo_ctx;

void init_x11evo_ctx()
{
#ifdef NO_AES_NI
     sph_groestl512_init( &x11evo_ctx.groestl );
     sph_echo512_init( &x11evo_ctx.echo );
#else
     init_echo( &x11evo_ctx.echo, 512 );
     init_groestl( &x11evo_ctx.groestl );
#endif
     init_luffa( &x11evo_ctx.luffa, 512 );
     cubehashInit( &x11evo_ctx.cube, 512, 16, 32 );
     init_sd( &x11evo_ctx.simd, 512 );
     sph_blake512_init( &x11evo_ctx.blake );
     sph_bmw512_init( &x11evo_ctx.bmw );
     sph_skein512_init( &x11evo_ctx.skein );
     sph_jh512_init( &x11evo_ctx.jh );
     sph_keccak512_init( &x11evo_ctx.keccak );
     sph_shavite512_init( &x11evo_ctx.shavite );
}

uint32_t getCurrentAlgoSeq(uint32_t current_time, uint32_t base_time)
{
	return (current_time - base_time) / (60 * 60 * 24);
}

void swap( uint8_t *a, uint8_t *b )
{
	uint8_t __tmp = *a;
	*a = *b;
	*b = __tmp;
}

void initPerm( uint8_t n[], uint8_t count )
{
	int i;
	for ( i = 0; i<count; i++ )
		n[i] = i;
}

int nextPerm( uint8_t n[], uint32_t count )
{
	uint32_t tail, i, j;

	if (unlikely( count <= 1 ))
		return 0;

	for ( i = count - 1; i>0 && n[i - 1] >= n[i]; i-- );
           tail = i;

	if ( tail > 0 )
            for ( j = count - 1; j>tail && n[j] <= n[tail - 1]; j-- );
	         swap( &n[tail - 1], &n[j] );

	for ( i = tail, j = count - 1; i<j; i++, j-- )
		swap( &n[i], &n[j] );

	return ( tail != 0 );
}

void getAlgoString( char *str, uint32_t count )
{
	uint8_t algoList[HASH_FUNC_COUNT];
	char *sptr;
        int j;
        int k;
	initPerm( algoList, HASH_FUNC_COUNT );

	for ( k = 0; k < count; k++ )
		nextPerm( algoList, HASH_FUNC_COUNT );

	sptr = str;
	for ( j = 0; j < HASH_FUNC_COUNT; j++ )
        {
		if ( algoList[j] >= 10 )
			sprintf( sptr, "%c", 'A' + (algoList[j] - 10) );
		else
			sprintf( sptr, "%u", algoList[j] );
		sptr++;
	}
	*sptr = 0;

	//applog(LOG_DEBUG, "nextPerm %s", str);
}

void evocoin_twisted_code( char *result, char *code )
{
	uint32_t h32, *be32 = get_stratum_job_ntime();
	h32 = be32toh(*be32);
	
	uint32_t count = getCurrentAlgoSeq(h32, INITIAL_DATE);
	getAlgoString(code, count);
	sprintf(result, "_%d_%s_", count, code);
}

static inline void x11evo_hash( void *state, const void *input )
{
   uint32_t hash[16];
   char completeCode[64];
   char resultCode[HASH_FUNC_COUNT + 1];
   x11evo_ctx_holder ctx;
   memcpy( &ctx, &x11evo_ctx, sizeof(x11evo_ctx) );
   evocoin_twisted_code( completeCode, resultCode );

   int i;
   for ( i = 0; i < strlen(resultCode); i++ )
   {
	char elem = resultCode[i];
	uint8_t idx;
	if (elem >= 'A')
		idx = elem - 'A' + 10;
	else
		idx = elem - '0';

	int size = 64;

	switch (idx)
        {
           case 0:
	      sph_blake512( &ctx.blake, (char*)input, 80 );
	      sph_blake512_close( &ctx.blake, (char*)hash );
	      break;
	   case 1:
	      sph_bmw512( &ctx.bmw, (char*)hash, size );
	      sph_bmw512_close( &ctx.bmw, (char*)hash );
	      break;
	   case 2:
#ifdef NO_AES_NI
	      sph_groestl512( &ctx.groestl, (char*)hash, size );
	      sph_groestl512_close( &ctx.groestl, (char*)hash );
#else
              update_groestl( &ctx.groestl, (char*)hash, 512 );
              final_groestl( &ctx.groestl, (char*)hash );
#endif
	      break;
	    case 3:
	      sph_skein512( &ctx.skein, (char*)hash, size );
	      sph_skein512_close( &ctx.skein, (char*)hash );
	      break;
	    case 4:
	      sph_jh512( &ctx.jh, (char*)hash, size );
	      sph_jh512_close( &ctx.jh, (char*)hash );
	      break;
	    case 5:
	      sph_keccak512( &ctx.keccak, (char*)hash, size );
	      sph_keccak512_close( &ctx.keccak, (char*)hash );
	      break;
	    case 6:
              update_luffa( &ctx.luffa, (char*)hash, 512 );
              final_luffa( &ctx.luffa, (char*)hash );
	      break;
	    case 7:
              cubehashUpdate( &ctx.cube, (char*)hash, 64 );
              cubehashDigest( &ctx.cube, (char*)hash );
	      break;
	    case 8:
	      sph_shavite512( &ctx.shavite, (char*)hash, size );
	      sph_shavite512_close( &ctx.shavite, (char*)hash );
	      break;
	    case 9:
              update_sd( &ctx.simd, (char*)hash, 512 );
              final_sd( &ctx.simd, (char*)hash );
	      break;
	    case 10:
#ifdef NO_AES_NI
	      sph_echo512( &ctx.echo, (char*)hash, size );
	      sph_echo512_close( &ctx.echo, (char*)hash );
#else
              update_echo( &ctx.echo, (char*)hash, 512 );
              final_echo( &ctx.echo, (char*)hash );
#endif
	      break;
	}
    }
    memcpy( state, hash, 32 );
}

static const uint32_t diff1targ = 0x0000ffff;

int scanhash_x11evo( int thr_id, struct work* work, uint32_t max_nonce,
                     unsigned long *hashes_done )
{
        uint32_t endiandata[20] __attribute__((aligned(64)));
        uint32_t hash64[8] __attribute__((aligned(32)));
        uint32_t *pdata = work->data;
        uint32_t *ptarget = work->target;
	uint32_t n = pdata[19] - 1;
	const uint32_t first_nonce = pdata[19];
        const uint32_t Htarg = ptarget[7];

        swab32_array( endiandata, pdata, 20 );

        uint32_t hmask = 0xFFFFFFFF;
        if ( Htarg  > 0 )
         if ( Htarg <= 0xF )
            hmask = 0xFFFFFFF0;
         else if ( Htarg <= 0xFF )
            hmask = 0xFFFFFF00;
         else if ( Htarg <= 0xFFF )
            hmask = 0xFFFF000;
         else if ( Htarg <= 0xFFFF )
            hmask = 0xFFFF000;

        do
        {
          pdata[19] = ++n;
          be32enc( &endiandata[19], n );
          x11evo_hash( hash64, &endiandata );
          if ( ( hash64[7] & hmask ) == 0 )
          {
             if ( fulltest( hash64, ptarget ) )
             {
                 *hashes_done = n - first_nonce + 1;
                 return true;
             }
           }
        } while ( n < max_nonce && !work_restart[thr_id].restart );

	*hashes_done = n - first_nonce + 1;
	pdata[19] = n;
	return 0;
}

bool register_x11evo_algo( algo_gate_t* gate )
{
  gate->aes_ni_optimized = true;
  gate->optimizations = SSE2_OPT | AES_OPT | AVX_OPT | AVX2_OPT;
  gate->scanhash  = (void*)&scanhash_x11evo;
  gate->hash      = (void*)&x11evo_hash;
  gate->hash_alt  = (void*)&x11evo_hash;
  init_x11evo_ctx();
  return true;
};

