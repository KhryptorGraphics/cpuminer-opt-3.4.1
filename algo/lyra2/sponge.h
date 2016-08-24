/**
 * Header file for Blake2b's internal permutation in the form of a sponge.
 * This code is based on the original Blake2b's implementation provided by
 * Samuel Neves (https://blake2.net/)
 *
 * Author: The Lyra PHC team (http://www.lyra-kdf.net/) -- 2014.
 *
 * This software is hereby placed in the public domain.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef SPONGE_H_
#define SPONGE_H_

#include <stdint.h>

/* Blake2b IV Array */
static const uint64_t blake2b_IV[8] =
{
	0x6a09e667f3bcc908ULL, 0xbb67ae8584caa73bULL,
	0x3c6ef372fe94f82bULL, 0xa54ff53a5f1d36f1ULL,
	0x510e527fade682d1ULL, 0x9b05688c2b3e6c1fULL,
	0x1f83d9abfb41bd6bULL, 0x5be0cd19137e2179ULL
};

/* Blake2b's rotation */
static __inline uint64_t rotr64(const uint64_t w, const unsigned c) {
#ifdef _MSC_VER
	return _rotr64(w, c);
#else
	return ( w >> c ) | ( w << ( 64 - c ) );
#endif
}

#if defined __AVX2__

// _m256i
#define  mm256_rotr_64(w,c) _mm256_or_si256(_mm256_srli_epi64(w, c), \
                                            _mm256_slli_epi64(w, 64 - c))

// Rotate uint64 by one uint64
// __m256i
#define mm256_rotl256_1x64(s) _mm256_permute4x64_epi64( s, 0x39 )
#define mm256_rotr256_1x64(s) _mm256_permute4x64_epi64( s, 0x93 )

// swap hi and lo 128 bits in 256 bit vector
// _m256i
#define mm256_swap128(s) _mm256_permute2f128_si256( s, s, 1 )

// void
#define G_4X64(a,b,c,d) \
   a = _mm256_add_epi64( a, b ); \
   d = mm256_rotr_64( _mm256_xor_si256( d, a), 32 ); \
   c = _mm256_add_epi64( c, d ); \
   b = mm256_rotr_64( _mm256_xor_si256( b, c ), 24 ); \
   a = _mm256_add_epi64( a, b ); \
   d = mm256_rotr_64( _mm256_xor_si256( d, a ), 16 ); \
   c = _mm256_add_epi64( c, d ); \
   b = mm256_rotr_64( _mm256_xor_si256( b, c ), 63 );

#elif defined __AVX__

// _m128i
#define  mm_rotr_64(w,c) _mm_or_si128(_mm_srli_epi64(w, c), \
                                      _mm_slli_epi64(w, 64 - c))

// swap 128 bit source vectors
// void
#define mm128_swap128(s0, s1) s0 = _mm_xor_si128(s0, s1); \
                              s1 = _mm_xor_si128(s0, s1); \
                              s0 = _mm_xor_si128(s0, s1);
 
// swap uint64 in source vector
// __m128i
#define mm128_swap64(s) _mm_or_si128( _mm_slli_si128( s, 8 ), \
                                      _mm_srli_si128( s, 8 ) )

// rotate 2 128 bit vectors as one 256 vector by 1 uint64
//void
#define mm128_rotl256_1x64(s0, s1) do { \
   __m128i t; \
   s0 = mm128_swap64( s0); \
   s1 = mm128_swap64( s1); \
   t = _mm_or_si128( _mm_and_si128( s0, _mm_set_epi64x(0ull,0xffffffffffffffffull) ), \
                     _mm_and_si128( s1, _mm_set_epi64x(0xffffffffffffffffull,0ull) ) ); \
   s1 = _mm_or_si128( _mm_and_si128( s0, _mm_set_epi64x(0xffffffffffffffffull,0ull) ), \
                      _mm_and_si128( s1, _mm_set_epi64x(0ull,0xffffffffffffffffull) ) ); \
   s0 = t; \
} while(0)

#define mm128_rotr256_1x64(s0, s1) do { \
   __m128i t; \
   s0 = mm128_swap64( s0); \
   s1 = mm128_swap64( s1); \
   t = _mm_or_si128( _mm_and_si128( s0, _mm_set_epi64x(0xffffffffffffffffull,0ull) ), \
                        _mm_and_si128( s1, _mm_set_epi64x(0ull,0xffffffffffffffffull) ) ); \
   s1 = _mm_or_si128( _mm_and_si128( s0, _mm_set_epi64x(0ull,0xffffffffffffffffull) ), \
                      _mm_and_si128( s1, _mm_set_epi64x(0xffffffffffffffffull,0ull) ) ); \
   s0 = t; \
} while(0)


#define G_2X64(a,b,c,d) \
   a = _mm_add_epi64( a, b ); \
   d = mm_rotr_64( _mm_xor_si128( d, a), 32 ); \
   c = _mm_add_epi64( c, d ); \
   b = mm_rotr_64( _mm_xor_si128( b, c ), 24 ); \
   a = _mm_add_epi64( a, b ); \
   d = mm_rotr_64( _mm_xor_si128( d, a ), 16 ); \
   c = _mm_add_epi64( c, d ); \
   b = mm_rotr_64( _mm_xor_si128( b, c ), 63 );

#endif   // AVX2

/* Blake2b's G function */
#define G(r,i,a,b,c,d) do { \
	a = a + b; \
	d = rotr64(d ^ a, 32); \
	c = c + d; \
	b = rotr64(b ^ c, 24); \
	a = a + b; \
	d = rotr64(d ^ a, 16); \
	c = c + d; \
	b = rotr64(b ^ c, 63); \
  } while(0)


/*One Round of the Blake2b's compression function*/
#define ROUND_LYRA(r) \
	G(r,0,v[ 0],v[ 4],v[ 8],v[12]); \
	G(r,1,v[ 1],v[ 5],v[ 9],v[13]); \
	G(r,2,v[ 2],v[ 6],v[10],v[14]); \
	G(r,3,v[ 3],v[ 7],v[11],v[15]); \
	G(r,4,v[ 0],v[ 5],v[10],v[15]); \
	G(r,5,v[ 1],v[ 6],v[11],v[12]); \
	G(r,6,v[ 2],v[ 7],v[ 8],v[13]); \
	G(r,7,v[ 3],v[ 4],v[ 9],v[14]);


//---- Housekeeping
void initState(uint64_t state[/*16*/]);

//---- Squeezes
void squeeze(uint64_t *state, unsigned char *out, unsigned int len);
void reducedSqueezeRow0(uint64_t* state, uint64_t* row, const uint32_t nCols);

//---- Absorbs
void absorbBlock(uint64_t *state, const uint64_t *in);
void absorbBlockBlake2Safe(uint64_t *state, const uint64_t *in);

//---- Duplexes
void reducedDuplexRow1(uint64_t *state, uint64_t *rowIn, uint64_t *rowOut, const uint32_t nCols);
void reducedDuplexRowSetup(uint64_t *state, uint64_t *rowIn, uint64_t *rowInOut, uint64_t *rowOut, const uint32_t nCols);
void reducedDuplexRow(uint64_t *state, uint64_t *rowIn, uint64_t *rowInOut, uint64_t *rowOut, const uint32_t nCols);

//---- Misc
void printArray(unsigned char *array, unsigned int size, char *name);

#endif /* SPONGE_H_ */
