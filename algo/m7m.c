
#include "cpuminer-config.h"
#include "miner.h"
#include "algo-gate-api.h"

#include <gmp.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "algo/sha2/sph-sha2.h"
#include "algo/keccak/sph_keccak.h"
#include "algo/haval/sph-haval.h"
#include "algo/tiger/sph_tiger.h"
#include "algo/whirlpool/sph_whirlpool.h"
#include "algo/ripemd/sph_ripemd.h"

#define EPS1 DBL_EPSILON
#define EPS2 3.0e-11

inline double exp_n(double xt)
{
    if(xt < -700.0)
        return 0;
    else if(xt > 700.0)
        return 1e200;
    else if(xt > -0.8e-8 && xt < 0.8e-8)
        return (1.0 + xt);
    else
        return exp(xt);
}

// 1 / (1 + exp(x1-x2))
inline double exp_n2(double x1, double x2)
{
    double p1 = -700., p2 = -37., p3 = -0.8e-8, p4 = 0.8e-8, p5 = 37., p6 = 700.;
    double xt = x1 - x2;
    if (xt < p1+1.e-200)
        return 1.;
    else if (xt > p1 && xt < p2 + 1.e-200)
        return ( 1. - exp(xt) );
    else if (xt > p2 && xt < p3 + 1.e-200)
        return ( 1. / (1. + exp(xt)) );
    else if (xt > p3 && xt < p4)
        return ( 1. / (2. + xt) );
    else if (xt > p4 - 1.e-200 && xt < p5)
        return ( exp(-xt) / (1. + exp(-xt)) );
    else if (xt > p5 - 1.e-200 && xt < p6)
        return ( exp(-xt) );
    else if (xt > p6 - 1.e-200)
        return 0.;
    
    //return(1.0 / (1.0 + exp(x1 - x2)));
}

double swit_(double wvnmb)
{
    return pow( (5.55243*(exp_n(-0.3*wvnmb/15.762) - exp_n(-0.6*wvnmb/15.762)))*wvnmb, 0.5) 
	  / 1034.66 * pow(sin(wvnmb/65.), 2.);
}

double m7mGaussianQuad_N(const double x1, const double x2)
{
    double s=0.0;
    double x[6], w[6];
    //gauleg(a2, b2, x, w);
    
    int m,j;
    double z1, z, xm, xl, pp, p3, p2, p1;
    m=3;
    xm=0.5*(x2+x1);
    xl=0.5*(x2-x1);
    for(int i=1;i<=3;i++)
    {
		z = (i == 1) ? 0.909632 : -0.0;
		z = (i == 2) ? 0.540641 : z;
	    do
	    {
			p1 = z;
			p2 = 1;
			p3 = 0;
			
			p3=1;
			p2=z;
			p1=((3.0 * z * z) - 1) / 2;
			
			p3=p2;
			p2=p1;
			p1=((5.0 * z * p2) - (2.0 * z)) / 3;
			
			p3=p2;
			p2=p1;
			p1=((7.0 * z * p2) - (3.0 * p3)) / 4;
			
			p3=p2;
			p2=p1;
			p1=((9.0 * z * p2) - (4.0 * p3)) / 5;
		    
		    pp=5*(z*p1-p2)/(z*z-1.0);
		    z1=z;
		    z=z1-p1/pp;
	    } while (fabs(z-z1) > 3.0e-11);
	    
	    x[i]=xm-xl*z;
	    x[5+1-i]=xm+xl*z;
	    w[i]=2.0*xl/((1.0-z*z)*pp*pp);
	    w[5+1-i]=w[i];
    }
    
    for(int j=1; j<=5; j++) s += w[j]*swit_(x[j]);
    
    return s;
}

uint32_t m7m_sw_(int nnounce)
{
    double wmax = ((sqrt((double)(nnounce))*(1.+EPS1))/450+100);
    return ((uint32_t)(m7mGaussianQuad_N(0., wmax)*(1.+EPS1)*1.e6));
}

static void mpz_set_uint256(mpz_t r, uint8_t *u)
{
    mpz_import(r, 32 / sizeof(unsigned long), -1, sizeof(unsigned long), -1, 0, u);
}

static void mpz_get_uint256(mpz_t r, uint8_t *u)
{
    u=0;
    mpz_export(u, 0, -1, sizeof(unsigned long), -1, 0, r);
}

static void mpz_set_uint512(mpz_t r, uint8_t *u)
{
    mpz_import(r, 64 / sizeof(unsigned long), -1, sizeof(unsigned long), -1, 0, u);
}

static void set_one_if_zero(uint8_t *hash512) {
    for (int i = 0; i < 32; i++) {
        if (hash512[i] != 0) {
            return;
        }
    }
    hash512[0] = 1;
}

static bool fulltest_m7hash(const uint32_t *hash32, const uint32_t *target32)
{
    int i;
    bool rc = true;

    const unsigned char *hash = (const unsigned char *)hash32;
    const unsigned char *target = (const unsigned char *)target32;
    for (i = 31; i >= 0; i--) {
        if (hash[i] != target[i]) {
            rc = hash[i] < target[i];
            break;
        }
    }

    return rc;
}

typedef struct {
    sph_sha256_context       sha256;
    sph_sha512_context       sha512;
    sph_keccak512_context    keccak;
    sph_whirlpool_context    whirlpool;
    sph_haval256_5_context   haval;
    sph_tiger_context        tiger;
    sph_ripemd160_context    ripemd;
} m7m_ctx_holder;

sph_sha256_context m7m_ctx_final_sha256;
m7m_ctx_holder m7m_ctx;

void init_m7m_ctx()
{
    sph_sha256_init( &m7m_ctx_final_sha256 );
    sph_sha256_init( &m7m_ctx.sha256 );
    sph_sha512_init( &m7m_ctx.sha512 );
    sph_keccak512_init( &m7m_ctx.keccak );
    sph_whirlpool_init( &m7m_ctx.whirlpool );
    sph_haval256_5_init( &m7m_ctx.haval );
    sph_tiger_init( &m7m_ctx.tiger );
    sph_ripemd160_init( &m7m_ctx.ripemd );
}

#define BITS_PER_DIGIT 3.32192809488736234787
#define EPS (DBL_EPSILON)

#define NM7M 5
#define SW_DIVS 5
#define M7_MIDSTATE_LEN 76
int scanhash_m7m_hash( int thr_id, struct work* work,
                       uint64_t max_nonce, unsigned long *hashes_done )
{
    uint32_t *pdata = work->data;
    uint32_t *ptarget = work->target;
    uint32_t data[32] __attribute__((aligned(128)));
    uint32_t *data_p64 = data + (M7_MIDSTATE_LEN / sizeof(data[0]));
    uint32_t hash[8] __attribute__((aligned(32)));
    uint8_t bhash[7][64] __attribute__((aligned(32)));
    uint32_t n = pdata[19] - 1;
    const uint32_t first_nonce = pdata[19];
    char data_str[161], hash_str[65], target_str[65];
    //uint8_t *bdata = 0;
    uint8_t bdata[8192];
    int rc = 0;
    int bytes, nnNonce2;

    m7m_ctx_holder ctx1, ctx2;
    memcpy( &ctx1, &m7m_ctx, sizeof(m7m_ctx) );
    sph_sha256_context ctxf_sha256;
    memcpy( &ctxf_sha256, &m7m_ctx_final_sha256, sizeof(ctxf_sha256) );
    
    memcpy(data, pdata, 80);

    sph_sha256( &ctx1.sha256, data, M7_MIDSTATE_LEN );
    sph_sha512( &ctx1.sha512, data, M7_MIDSTATE_LEN );
    sph_keccak512( &ctx1.keccak, data, M7_MIDSTATE_LEN );
    sph_whirlpool( &ctx1.whirlpool, data, M7_MIDSTATE_LEN );
    sph_haval256_5( &ctx1.haval, data, M7_MIDSTATE_LEN );
    sph_tiger( &ctx1.tiger, data, M7_MIDSTATE_LEN );
    sph_ripemd160( &ctx1.ripemd, data, M7_MIDSTATE_LEN );

    mpz_t magipi, magisw, product, bns[8];
    mpf_t magifpi, mpt1, mpt2, mptmp, mpten;

    mpz_inits( magipi, magisw, bns[0], bns[1], bns[2], bns[3],
                               bns[4], bns[5], bns[6], bns[7], NULL );
    mpz_init2(product, 512);
	
    int digits=(int)((sqrt((double)(INT_MAX))*(1.+EPS))/9000+75);
    mpf_set_default_prec((long int)(digits*BITS_PER_DIGIT+16));
	
    mpf_init(magifpi);
    mpf_init(mpt1);
    mpf_init(mpt2);
    mpf_init(mptmp);
    mpf_init_set_ui(mpten, 10);
	
    do
    {
        data[19] = ++n;
        nnNonce2 = (int)(data[19]/2);
        memset(bhash, 0, 7 * 64);

        memcpy( &ctx2, &ctx1, sizeof(m7m_ctx) );

        sph_sha256( &ctx2.sha256, data_p64, 80 - M7_MIDSTATE_LEN );
        sph_sha256_close( &ctx2.sha256, (void*)(bhash[0]) );

        sph_sha512( &ctx2.sha512, data_p64, 80 - M7_MIDSTATE_LEN );
        sph_sha512_close( &ctx2.sha512, (void*)(bhash[1]) );
        
        sph_keccak512( &ctx2.keccak, data_p64, 80 - M7_MIDSTATE_LEN );
        sph_keccak512_close( &ctx2.keccak, (void*)(bhash[2]) );

        sph_whirlpool( &ctx2.whirlpool, data_p64, 80 - M7_MIDSTATE_LEN );
        sph_whirlpool_close( &ctx2.whirlpool, (void*)(bhash[3]) );
        
        sph_haval256_5( &ctx2.haval, data_p64, 80 - M7_MIDSTATE_LEN );
        sph_haval256_5_close( &ctx2.haval, (void*)(bhash[4])) ;

        sph_tiger( &ctx2.tiger, data_p64, 80 - M7_MIDSTATE_LEN );
        sph_tiger_close( &ctx2.tiger, (void*)(bhash[5]) );

        sph_ripemd160( &ctx2.ripemd, data_p64, 80 - M7_MIDSTATE_LEN );
        sph_ripemd160_close( &ctx2.ripemd, (void*)(bhash[6]) );

        for( int i=0; i < 7; i++ )
        {
            set_one_if_zero(bhash[i]);
            mpz_set_uint512(bns[i],bhash[i]);
        }

        //mpz_set_ui(bns[7],0);
        mpz_set(bns[7], bns[0]);

        for(int i=1; i < 7; i++)
	        mpz_add(bns[7], bns[7], bns[i]);

        //mpz_set_ui(product,1);
	mpz_set(product, bns[0]);
		
        for(int i=1; i < 8; i++)
            mpz_mul(product, product, bns[i]);
		
        //mpz_pow_ui(product, product, 2); 
	mpz_mul(product, product, product);
		
        bytes = mpz_sizeinbase(product, 256);
        //bdata = (uint8_t *)realloc(bdata, bytes);
        mpz_export((void *)bdata, NULL, -1, 1, 0, 0, product);

        sph_sha256( &ctxf_sha256, bdata, bytes );
        sph_sha256_close( &ctxf_sha256, (void*)(hash) );

        int digits=(int)((sqrt((double)(nnNonce2))*(1.+EPS))/9000+75);
        int iterations=20;
        //mpf_set_default_prec((long int)(digits*BITS_PER_DIGIT+16));
        mp_bitcnt_t prec = (long int)(digits*BITS_PER_DIGIT+16);
	mpf_set_prec_raw(magifpi, prec);
	mpf_set_prec_raw(mptmp, prec);
	mpf_set_prec_raw(mpt1, prec);
	mpf_set_prec_raw(mpt2, prec);

        uint32_t usw_ = m7m_sw_(nnNonce2);
        //if (usw_ < 1) usw_ = 1;
        usw_ = (usw_ == 0) ? 1 : usw_; 
        mpz_set_ui(magisw, usw_);
        //uint32_t mpzscale=mpz_size(magisw);
		
	uint32_t mpzscale = mpz_size(magisw);
	mpzscale = (mpzscale == 0) ? 1 : ((mpzscale > 1000) ? 1000 : mpzscale);
		
        for(int i = 0; i < 5; i++)
        {	
			//mpzscale = (mpzscale < 1) ? 1 : ((mpzscale > 1000) ? 1000 : mpzscale);
			
            mpf_set_d(mpt2, 0.25*mpzscale);
	    mpf_set_str(mpt1, "0.61f78a9abaa58b4698916152cf7eee1bbdf1f5b4ab3de24c3a3c159062718f716d12d59bba9c4881f@0", 16);
	    mpf_sub(mpt1, mpt2, mpt1);
			
	    mpf_set_str(mpt2, "0.28bb03352962950bc65974466fd2c021e28cc50addfa74e94a3ec525f1887d66104fe025ccb43562b@0", 16);
	    mpf_sub(mpt2, mpt1, mpt2);
	    mpf_swap(mpt1, mpt2);
			
	    mpf_set_str(mpt2, "0.3852f236e8bf02d5c6718fc775e190c575b20b91a45e1338761547573ad5e77e5f61fc760362c231@-1", 16);
	    mpf_sub(mpt2, mpt1, mpt2);
	    mpf_swap(mpt1, mpt2);
		
	    mpf_set_str(mpt2, "0.35da65929531617012a859be186b093a9f835b39b8f14f38ac8ce59915e85f4497b3fe99c7e44e44e@-3", 16);
	    mpf_sub(mpt2, mpt1, mpt2);
	    mpf_swap(mpt1, mpt2);
			
	    mpf_set_str(mpt2, "0.189daa22ce9d43d42f773cc46e4c97fc8f625df4c82a663970a52968d90df7aa3e05294f26eec633@-7", 16);
	    mpf_sub(mpt2, mpt1, mpt2);
	    mpf_swap(mpt1, mpt2);
			
	    mpf_set_str(mpt2, "0.87127377035941932784771446314900974221284468466736176188550809681111589251132779066239512337e-20", 10);
	    mpf_sub(mpt2, mpt1, mpt2);
	    mpf_swap(mpt1, mpt2);

            mpf_abs(mpt1, mpt1);
            mpf_set_str(magifpi, "0.717770011046129997821193223665779426657129889339984371989763663877269423125849866370161623127786e0", 10);
            mpf_div(magifpi, magifpi, mpt1);
			
            mpf_pow_ui(mptmp, mpten, digits >> 1);
            mpf_mul(magifpi, magifpi, mptmp);
			
			//mpf_mul_ui(magifpi, magifpi, pow(10.0, (double)(digits >> 1)));
			
            mpz_set_f(magipi, magifpi);
						
            mpz_add(product,product,magipi);
            mpz_add(product,product,magisw);
			
            mpz_set_uint256(bns[0], (void*)(hash));
            mpz_add(bns[7], bns[7], bns[0]);

            mpz_mul(product,product,bns[7]);
            mpz_cdiv_q (product, product, bns[0]);
            if (mpz_sgn(product) <= 0) mpz_set_ui(product,1);

            bytes = mpz_sizeinbase(product, 256);
            mpzscale=bytes;
            //bdata = (uint8_t *)realloc(bdata, bytes);
            mpz_export(bdata, NULL, -1, 1, 0, 0, product);
			
            sph_sha256( &ctxf_sha256, bdata, bytes );
            sph_sha256_close( &ctxf_sha256, (void*)(hash) );
	}

        rc = fulltest_m7hash(hash, ptarget);
        if (unlikely(rc)) {
            if (opt_debug) {
                bin2hex(hash_str, (unsigned char *)hash, 32);
                bin2hex(target_str, (unsigned char *)ptarget, 32);
                bin2hex(data_str, (unsigned char *)data, 80);
                applog(LOG_DEBUG, "DEBUG: [%d thread] Found share!\ndata   %s\nhash   %s\ntarget %s", thr_id, 
                    data_str,
                    hash_str,
                    target_str);
            }

            pdata[19] = data[19];

            goto out;
        }
    } while (n < max_nonce && !work_restart[thr_id].restart);

    pdata[19] = n;
	
out:
	digits = (int)((sqrt((double)(INT_MAX))*(1.+EPS))/9000+75);
	mp_bitcnt_t prec = (long int)(digits*BITS_PER_DIGIT+16);
	mpf_set_prec_raw(magifpi, prec);
	mpf_set_prec_raw(mptmp, prec);
	mpf_set_prec_raw(mpt1, prec);
	mpf_set_prec_raw(mpt2, prec);
	mpf_clear(magifpi);
	mpf_clear(mpten);
	mpf_clear(mptmp);
	mpf_clear(mpt1);
	mpf_clear(mpt2);
    
	mpz_clears(magipi, magisw, product, bns[0], bns[1], bns[2], bns[3], bns[4], bns[5], bns[6], bns[7], NULL);

    //free(bdata);

    *hashes_done = n - first_nonce + 1;
    return rc;
}

void m7m_reverse_endian( struct work *work )
{
   swab32_array( work->data, work->data, 20 );
}

bool register_m7m_algo( algo_gate_t *gate )
{
  gate->aes_ni_optimized      = true;
  gate->optimizations = SSE2_OPT | AES_OPT | AVX_OPT;
  init_m7m_ctx();
  gate->scanhash              = (void*)scanhash_m7m_hash;
  gate->build_stratum_request = (void*)&std_be_build_stratum_request;
  gate->set_target            = (void*)&scrypt_set_target;
  gate->get_max64             = (void*)&get_max64_0x1ffff;
  gate->set_work_data_endian  = (void*)&m7m_reverse_endian;
  gate->work_data_size        = 80;
  return true;
}

