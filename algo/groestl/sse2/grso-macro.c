/* hash.c     January 2011
 *
 * Groestl-512 implementation with inline assembly containing mmx and
 * sse instructions. Optimized for Opteron.
 * Authors: Krystian Matusiewicz and Soeren S. Thomsen
 *
 * This code is placed in the public domain
 */

//#include "grso.h"
//#include "grso-asm.h"
// #include "grsotab.h"

#define DECL_GRS

/* load initial constants */
#define GRS_I \
do { \
  int i; \
  /* set initial value */ \
  for (i = 0; i < grsoCOLS-1; i++) sts_grs.grsstate[i] = 0; \
  sts_grs.grsstate[grsoCOLS-1] = grsoU64BIG((u64)(8*grsoDIGESTSIZE)); \
 \
  /* set other variables */ \
  sts_grs.grsbuf_ptr = 0; \
  sts_grs.grsblock_counter = 0; \
} while (0); \

/* load hash */
#define GRS_U \
do { \
    unsigned char* in = hash; \
  unsigned long long index = 0; \
 \
  /* if the buffer contains data that has not yet been digested, first \
     add data to buffer until full */ \
  if (sts_grs.grsbuf_ptr) { \
    while (sts_grs.grsbuf_ptr < grsoSIZE && index < 64) { \
      hashbuf[(int)sts_grs.grsbuf_ptr++] = in[index++]; \
    } \
    if (sts_grs.grsbuf_ptr < grsoSIZE) continue; \
 \
    /* digest buffer */ \
    sts_grs.grsbuf_ptr = 0; \
    grsoTransform(&sts_grs, hashbuf, grsoSIZE); \
  } \
 \
  /* digest bulk of message */ \
  grsoTransform(&sts_grs, in+index, 64-index); \
  index += ((64-index)/grsoSIZE)*grsoSIZE; \
 \
  /* store remaining data in buffer */ \
  while (index < 64) { \
    hashbuf[(int)sts_grs.grsbuf_ptr++] = in[index++]; \
  } \
 \
} while (0);

/* groestl512 hash loaded */
/* hash = groestl512(loaded) */
#define GRS_C \
do { \
    char *out = hash; \
  int i, j = 0; \
  unsigned char *s = (unsigned char*)sts_grs.grsstate; \
 \
  hashbuf[sts_grs.grsbuf_ptr++] = 0x80; \
 \
  /* pad with '0'-bits */ \
  if (sts_grs.grsbuf_ptr > grsoSIZE-grsoLENGTHFIELDLEN) { \
    /* padding requires two blocks */ \
    while (sts_grs.grsbuf_ptr < grsoSIZE) { \
      hashbuf[sts_grs.grsbuf_ptr++] = 0; \
    } \
    /* digest first padding block */ \
    grsoTransform(&sts_grs, hashbuf, grsoSIZE); \
    sts_grs.grsbuf_ptr = 0; \
  } \
  while (sts_grs.grsbuf_ptr < grsoSIZE-grsoLENGTHFIELDLEN) { \
    hashbuf[sts_grs.grsbuf_ptr++] = 0; \
  } \
 \
  /* length padding */ \
  sts_grs.grsblock_counter++; \
  sts_grs.grsbuf_ptr = grsoSIZE; \
  while (sts_grs.grsbuf_ptr > grsoSIZE-grsoLENGTHFIELDLEN) { \
    hashbuf[--sts_grs.grsbuf_ptr] = (unsigned char)sts_grs.grsblock_counter; \
    sts_grs.grsblock_counter >>= 8; \
  } \
 \
  /* digest final padding block */ \
  grsoTransform(&sts_grs, hashbuf, grsoSIZE); \
  /* perform output transformation */ \
  grsoOutputTransformation(&sts_grs); \
 \
  /* store hash result in output */ \
  for (i = grsoSIZE-grsoDIGESTSIZE; i < grsoSIZE; i++,j++) { \
    out[j] = s[i]; \
  } \
 \
  /* zeroise relevant variables and deallocate memory */ \
  for (i = 0; i < grsoCOLS; i++) { \
    sts_grs.grsstate[i] = 0; \
  } \
  for (i = 0; i < grsoSIZE; i++) { \
    hashbuf[i] = 0; \
  } \
} while (0); 
 

