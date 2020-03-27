/*
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 Copyright (C) 2017 -- Martin Schlather

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/


//  #define LANDSLIDE 1
//#include "chb.h"
// #include <inttypes.h> // uintptr_t
#include "chaploblocker.h"
#include "zzz_CHaploblocker.h"
#include <R_ext/Rdynload.h>
#include "error.h"
#include <zzz_RandomFieldsUtils.h>

CALL4(void, orderingInt, int*, data, int, len, int, dim, int *, pos)


#define LOG std::log
#define STRNCOPY std::strncpy // strncpy(dest, src, n
#define CEIL(X) std::ceil((double) X) // keine Klammern um X!
//#define STRNCMP std::strncmp

#define BitsPerUnit 32
#define Int int  // uint32_T
// #define Int uint64_t
#define BitsPerBlock 16
#define MaxUserValues 256
#define MaxChar 8

#define START 0
#define END 1
#define CODENROW 2
#define STARTREMAINDER 3
#define ENDREMAINDER 4
#define CODENCOL 5
#define SNPxIND 6
#define POSLAST SNPxIND


#define Uint unsigned int
Uint bitsPerCode, codesPerUnit, leadingPattern, leadingShift,
  MaxUserCodesPerBlock;
typedef enum {numeric, character, local} coding_type_def;
coding_type_def coding_type;
bool BigEndian, // true iff 0x8000 -> [0] is != 0
  allowexceptions = true;

Uint
  ncodes = 0,
  ONES = 0xFFFFFFFF;
int *codes = NULL;
typedef char ZK[MaxChar +1];
static ZK *strcodes = NULL,
  ANY[2] = {"@", "ANY"}, //"ANY",
  DUMMY = "!;";

SEXP Position = R_NilValue,
  Where = R_NilValue,
  Counts = R_NilValue,
  Codings = R_NilValue;

#define BitsPerBlockValues (1 << BitsPerBlock)
int unequal[BitsPerBlockValues];

const char namecodes[] = "strings", namestrcodes[] = "numbers";
#define ENSURE(CODING) if (CODING == NULL) ERR1("Coding expected that contains %s only.", name##CODING)


void fixcodingSub(int nval) {
  assert(BitsPerBlockValues >= MaxUserValues);
  if (sizeof(int) != 4 ||
      sizeof(short int) != 2 ||
      BitsPerUnit != 32 ||
      BitsPerBlock != 16
      ) BUG;
  if (nval > MaxUserValues)
    ERR2("maximum number of 'values' is %d. Got %d.", MaxUserValues, nval);

 // bitsPerCode sollen Bitlange von Int exakt teilen
  // potentielle Verlust sollte irrelvant sein, da er nur fuer
  // Zahl der values groesser 4 auftreten kann
  if (nval < 2) {ERR("not enough values");}
  bitsPerCode =
    1L << ((Uint) CEIL( LOG(LOG((double) nval) / LOG2) / LOG2));
  codesPerUnit = (Uint) (BitsPerUnit / bitsPerCode);
  assert(codesPerUnit * bitsPerCode == BitsPerUnit);
  leadingShift = bitsPerCode * (codesPerUnit - 1);
  leadingPattern = (ONES >> leadingShift) << leadingShift;

  MaxUserCodesPerBlock = (int) (BitsPerBlock / bitsPerCode);
  if (MaxUserCodesPerBlock * bitsPerCode != BitsPerBlock) BUG;

  // determine how many unequals are in each code
  for (Uint i = 0; i<BitsPerBlockValues; i++) {
    Uint code = i << (BitsPerUnit - BitsPerBlock);
    unequal[i] = 0;
    for (Uint k=0; k<MaxUserCodesPerBlock; k++) {
      unequal[i] += (code & leadingPattern) != 0;
      code <<= bitsPerCode;
    }
    assert(unequal[i] <= (int) MaxUserCodesPerBlock);
  }

  if (Position == R_NilValue) Position = install("position");
  if (Codings == R_NilValue) Codings = install("codings");
  if (Counts == R_NilValue) Counts = install("counts");
  if (Where == R_NilValue) Where = install("where.to.find");
}

void fixcodingIntern(int *values, int nval) {
  if (codes != NULL) FREE(codes);
  if (nval == 1) {
    // strings !!!!
    coding_type = local;
    ncodes = std::abs(*values);
    allowexceptions = *values > 0;
  } else {
    coding_type = numeric;
    ncodes = nval;
    codes = (int*) MALLOC(sizeof(int) * ncodes);
    for (int i=0; i<(int) ncodes; i++) codes[i] = values[i];
  }
  fixcodingSub(ncodes);
}


void fixcodingIntern(SEXP values) {
  coding_type = character;
  int nval = length(values);
  fixcodingSub(nval);
  ncodes = nval;
  if (strcodes != NULL) FREE(strcodes);
  strcodes = (ZK*) MALLOC(sizeof(ZK) * nval);
  for (int i=0; i<nval; i++) {
    if (STRNCMP(ANY[1], CHAR(STRING_ELT(values, i)), MaxChar) == 0)
      STRNCOPY(strcodes[i], ANY[0], MaxChar);
    else STRNCOPY(strcodes[i], CHAR(STRING_ELT(values, i)), MaxChar);
    strcodes[i][MaxChar] = '\n';
  }
}

SEXP fixcoding(SEXP values) {
  if (TYPEOF(values) == STRSXP) fixcodingIntern(values);
  else fixcodingIntern(INTEGER(values), length(values));
  return R_NilValue;
}

void initHaploBlocker() {
  union {
    unsigned short a;
    unsigned char b[2];
  } ab;
  ab.a = 0xFF00;
  BigEndian = ab.b[0] != 0;
  //  assert(!BigEndian);
  int value = -2;
  fixcodingIntern(&value, 1);
  assert(codes != NULL || coding_type == local);
}

void redoCoding(SEXP M, int nrow, int ncol) {
  //  printf("ncodes = %d\n", ncodes);
  if (strcodes != NULL) FREE(strcodes);
  strcodes = (ZK *) MALLOC(sizeof(ZK) * nrow * ncodes);
  ZK *pcodes = strcodes;
  int count[MaxUserValues];
  for (int nr=0; nr<nrow; nr++, pcodes += ncodes) {
    ZK temp[MaxUserValues];
    int idx = nr,
      order[MaxUserValues],
      n=0;
    for (int nc=0; nc<ncol; nc++, idx+=nrow) {
      const char *m = CHAR(STRING_ELT(M, idx));
      int i;
      for (i=0; i<n; i++) {
	if (STRNCMP(temp[i], m, MaxChar) == 0) {
	  count[i]--;
	  break;
	}
      }
      if (i == n) { // not found
	if (n < MaxUserValues) {
	  STRNCOPY(temp[n], m, MaxChar);
	  temp[n][MaxChar] = '\n';
	  count[n] = -1;
	} else {
	  if (!allowexceptions)
	    ERR1("maximal number of different entries per SNP reached. Setting 'fixcoding(%d)' is a way out.", ncodes);
	  // else ignore what is coming up
	}
	n++;
     }
    }
    RU_orderingInt(count, n, 1, order);
    int i,
      endfor = n >= (int) ncodes ? (int) ncodes : n;
    for (i=0; i<endfor; i++) {
      STRNCOPY(pcodes[i], temp[order[i]], MaxChar);
      pcodes[i][MaxChar] = '\n';
    }
    int k;
    if (allowexceptions) {
      k = (n >= (int) ncodes ? (int) ncodes - 1 : n);
      STRNCOPY(pcodes[k], ANY[0], MaxChar);
      pcodes[k][MaxChar] = '\n';
      k++;
    } else k = i;
    for (i=k; i<(int) ncodes; i++) STRCPY(pcodes[i], DUMMY);
    //if (nr < 4) {for (i=0; i<k; i++) PRINTF("%s", pcodes[i]); PRINTF("\n");}
  }
}



void codeInner(SEXP M, Uint nrow, Uint ncol, Uint start, Uint *code, bool sXi
	       , int total
	       ) {

  Uint
    startremainder = start % codesPerUnit,
    codenrow = 1 + (nrow + startremainder - 1) / codesPerUnit,//sxp x indiv
     codencol = 1 + (ncol + startremainder - 1) / codesPerUnit,//not sxp x indiv
    endk0 = codesPerUnit - startremainder;
  ZK *pcodes = strcodes;
  if (coding_type == local) pcodes += start * ncodes;

#define LOOP(EQUAL, ENDOUTER, ENDINNER, CODEINCR, IINCR, IDXFACTOR, LASTI, NUM)\
    for (Uint nc=0; nc<ENDOUTER; nc++, code+=CODEINCR) {		\
      int idx = nc * IDXFACTOR;						\
      for (Uint i=0, n=0, endk = endk0; i<LASTI; i+=IINCR, endk=codesPerUnit){ \
	assert(i < (Uint) total);					\
	code[i] = 0;							\
	for (Uint k=0; k<endk && n < ENDINNER; k++, n++, idx+=IINCR) {	\
	  Uint l;							\
	  for (l=0; l<ncodes; l++) {					\
	    assert(idx < length(M));					\
	    bool equal = EQUAL;						\
	    if (equal) {						\
	      code[i] <<= bitsPerCode;					\
	      code[i] |= l;						\
	      break;							\
	    }								\
	  }								\
	  if (l >= ncodes) {						\
	    if (NUM) { ERR2("value '%d' is outside range of {0,..., %d}%", \
			    (int) m[idx], ncodes - 1);			\
	    } else {							\
	      Uint n_ncodes = n * ncodes;				\
	      if (coding_type == local &&				\
		!STRNCMP(pcodes[(ncodes - 1) + n_ncodes], DUMMY, MaxChar)) { \
		assert(!allowexceptions);				\
		l=0;							\
		while(STRNCMP(pcodes[l + n_ncodes], DUMMY, MaxChar)) l++; \
		STRNCOPY(pcodes[l+n_ncodes], CHAR(STRING_ELT(M,idx)),MaxChar); \
		pcodes[l + n_ncodes][MaxChar] = '\n';			\
		code[i] <<= bitsPerCode;				\
		code[i] |= l;						\
	      } else {							\
		ERR3("Not all SNP values (%dx%d) are recognized in %s coding.",\
		     n + 1, nc + 1, coding_type == local ? "local"	\
		     : coding_type == numeric ? "numeric" : "character" ); \
	      }								\
	    }								\
	  }								\
	}								\
      }									\
       Uint nn = LASTI * codesPerUnit - ENDINNER - startremainder;	\
      code[(LASTI - 1) * IINCR] <<=  (nn * bitsPerCode);		\
    }

#define EQU !STRNCMP(CHAR(STRING_ELT(M, idx)), pcodes[l], MaxChar) || !STRNCMP(ANY[0], pcodes[l], MaxChar);

#define EQULOCAL !STRNCMP(CHAR(STRING_ELT(M, idx)), pcodes[l + n * ncodes], MaxChar) || !STRNCMP(ANY[0], pcodes[l + n * ncodes], MaxChar); //if (ncol==1 || true) { printf("nr=%d l=%d %s %s idx=%d\n", n, l, CHAR(STRING_ELT(M, idx)), pcodes[l + n * ncodes], l + n * ncodes); }

  if (TYPEOF(M) == REALSXP) {
    ENSURE(codes);
    double *m = REAL(M);
      if (sXi) LOOP(m[idx] == codes[l], ncol, nrow, codenrow, 1, nrow, codenrow,
                    true)
      else LOOP(m[idx] == codes[l], nrow, ncol, 1, nrow, 1, codencol, true)
  } else if (TYPEOF(M) == INTSXP) {
    ENSURE(codes);
    int *m = INTEGER(M);
    if (sXi) LOOP(m[idx] == codes[l], ncol, nrow, codenrow, 1, nrow, codenrow,
                 true)
      else LOOP(m[idx] == codes[l], nrow, ncol, 1, nrow, 1, codencol, true)
  } else if (TYPEOF(M) == LGLSXP) {
    ENSURE(codes);
    int *m = LOGICAL(M);
    if (sXi) LOOP(m[idx] == codes[l], ncol, nrow, codenrow, 1, nrow, codenrow,
             true)
    else LOOP(m[idx] == codes[l], nrow, ncol, 1, nrow, 1, codencol, true)
  } else if (TYPEOF(M) == STRSXP) {
    int *m = NULL;
    ENSURE(strcodes);
    if (coding_type == character) {
      if (sXi) LOOP( EQU , ncol, nrow, codenrow, 1, nrow, codenrow, false)
      else LOOP( EQU , nrow, ncol, 1, nrow, 1, codencol, false)
    } else if (coding_type == local) {
      if (sXi) LOOP( EQULOCAL , ncol, nrow, codenrow, 1, nrow, codenrow, false)
      else LOOP( EQULOCAL, nrow, ncol, 1, nrow, 1, codencol, false)
    } else ERR("coding type mismatch.");
  } else ERR("incompatible type of 'M'");
}

/*
+ }
.ncodes = 7 start=0
ncodes = 7 start=5
 [1] 1364538113  337912914 1158873376  622998036 1141122100  859841602
 [7] 1376780595  621941522 1079136802 1428443714
[1]        532 1141122100  859841602 1376780595  620756992
attr(,"position")
[1]  5 33  5  5  1
E
*/

SEXP codeSNPs(SEXP M, SEXP Start, SEXP RedoCoding, SEXP SNPxINDIV) {
  // M always starts from the very beginning
  // start gives the shift for the coded M
   if (length(M) == 0) return R_NilValue;
  Uint ncol, nrow, codenrow, codencol,
    start = INTEGER(Start)[0],
    startremainder = start % codesPerUnit;
  int size;
  SEXP Code;
  bool ismatrix = isMatrix(M),
    snpxind = LOGICAL(SNPxINDIV)[0];

  if (ismatrix) {
    nrow = nrows(M);
    ncol = ncols(M);
  } else {
    snpxind = true;
    nrow = length(M);
    ncol = 1;
  }

  if (snpxind) {
    size = nrow;
    codenrow = 1 + (nrow + startremainder - 1) / codesPerUnit;
    codencol = ncol;
    if (coding_type == local && LOGICAL(RedoCoding)[0] &&
	(strcodes == NULL || ncol > 1)) {
      if (TYPEOF(M) != STRSXP) ERR("'redo = TRUE' (currently) only allowd for string valued matrices. If the option 'redo = TRUE' is needed for integer values matrices please contact the author.");
      //if (ncol == 1) ERR("'redo = TRUE', but a single SNP sequence is given");
      if (start != 0) ERR("recoding must start from the beginning (start=0).");
      redoCoding(M, nrow, ncol);
    }
  } else {
    ERR("'idx x snp' currently disabled. Please contact author if needed.");
    size = ncol;
    codencol = 1 + (ncol + startremainder - 1) / codesPerUnit;
    codenrow = nrow;
    if (coding_type == local && LOGICAL(RedoCoding)[0] && nrow > 1)
      ERR("local recoding is not defined for individual x SNP matrices");
  }

  if (ismatrix) PROTECT(Code = allocMatrix(INTSXP, codenrow, codencol));
  else PROTECT(Code = allocVector(INTSXP, codenrow));
  for (Uint i=0; i<codenrow * codencol; INTEGER(Code)[i++] = 0);

  codeInner(M, nrow, ncol, start, (Uint*) INTEGER(Code), snpxind, codenrow * codencol);

  SEXP Posvec;
  PROTECT(Posvec = allocVector(INTSXP, POSLAST + 1));
  int *posvec = INTEGER(Posvec);
  posvec[START] = start; // first SNP in SNP sequence
  //  posvec[END] = start + size - 1; // last SNP in SNP sequence
  posvec[END] = start + size - 1; // last SNP in SNP sequence
  posvec[CODENROW] = codenrow;// number of filled codes
  posvec[CODENCOL] = codencol;// number of filled codes
  posvec[SNPxIND] = snpxind;// number of filled codes
  posvec[STARTREMAINDER] = startremainder;// first used position in code
  posvec[ENDREMAINDER] = posvec[END] % codesPerUnit;
  setAttrib(Code, Position, Posvec);
  UNPROTECT(2);
  return Code;
}


SEXP decodeSNPs(SEXP CM) {
  SEXP Ans;
  SEXP Pos = getAttrib(CM, Position);
  int incr, cv_factor, k_factor,
    k=0,
    *pos = INTEGER(Pos),
    type = coding_type == numeric ? INTSXP : STRSXP;
  Uint  endnc,
    size = pos[END] - pos[START] + 1,
    codenrow = pos[CODENROW], // =nrows(CM) bzw length(CM)
    codencol = pos[CODENCOL], // =nrows(CM) bzw length(CM)
    startremainder = pos[STARTREMAINDER],
    *cm  = (Uint*) INTEGER(CM);
  bool snpxind = pos[SNPxIND];

  if (snpxind) {
    incr = 1;
    endnc = codencol;
    cv_factor = codenrow;
    k_factor = size;
  } else {
    incr = codenrow;
    endnc = codenrow;
    cv_factor = 1;
    k_factor = 1;
  }

  if (isMatrix(CM)) PROTECT(Ans = allocMatrix(type, snpxind ? size : codenrow,
					      snpxind ? codencol : size));
  else PROTECT(Ans = allocVector(type, size));

  if (coding_type == numeric) { ENSURE(codes); } else {  ENSURE(strcodes); }

  for (Uint nc = 0; nc < endnc; nc++) {
    Uint *cv = cm + nc * cv_factor,
      pattern = leadingPattern >> (startremainder * bitsPerCode),
      shift = leadingShift - startremainder * bitsPerCode;
    ZK *pcodes = strcodes + pos[START] * ncodes;
    k = nc * k_factor;
    for (Uint i=0; i<size; i++, k+=incr) {
      Uint idx = (*cv & pattern) >> shift;
      if (idx >= ncodes) BUG;
      if (coding_type == numeric) INTEGER(Ans)[k] = codes[idx];
      else {
	SET_STRING_ELT(Ans, k, mkChar(pcodes[idx]));
	if (coding_type == local) pcodes += ncodes;
      }
      if (shift > 0) {
	pattern >>= bitsPerCode;
	shift -= bitsPerCode;
      } else {
	pattern = leadingPattern;
	shift = leadingShift;
	cv+=incr;
      }
    }
  }

  UNPROTECT(1);
  return Ans;
}


bool notequal(Uint *a, Uint *b, Uint c_start, Uint c_end, Uint startremainder,
	      Uint endremainder) {
  Uint i = c_start;
  if (c_start == c_end) {
    //   printf(":");
    Uint shift = startremainder * bitsPerCode;
    return (((a[i] ^ b[i]) << shift)
	    >> (shift + (codesPerUnit - 1 - endremainder) * bitsPerCode)) != 0;
  }
  if ((a[i] ^ b[i]) << (startremainder * bitsPerCode)) return true;
  for (i++; i<c_end; i++) if (a[i] != b[i]) return true;
  if ((a[i] ^ b[i]) >> ((codesPerUnit - 1 - endremainder) * bitsPerCode))
    return true;
  return false;
}

typedef struct pattern {
  int where, ith, count;
  pattern *next;
} pattern;




SEXP factorSNPs(SEXP CM, SEXP Start, SEXP End) {
  SEXP Pos = getAttrib(CM, Position);
  int *pos = INTEGER(Pos);
  if (!pos[SNPxIND]) ERR("'factorSNPs' only defined for SNPxINDIVID matrices");
  if (INTEGER(Start)[0] < 0) ERR("value of 'start' must be positive.");
  Uint addtofirst,
    start = (Uint) INTEGER(Start)[0],
    end =  (Uint) INTEGER(End)[0],
    c_start = start / codesPerUnit,
    c_end = end / codesPerUnit,
    startremainder = start % codesPerUnit,
    endremainder = end % codesPerUnit,
    *cm = ((Uint*) INTEGER(CM)),
    *cur_cm = cm,
    nrow = nrows(CM),
    ncol = ncols(CM)
    ;
  short Uint *first,
    SHORT_ONES = 0xFFFF,
    firstpattern = SHORT_ONES;

  SEXP Factor;
  PROTECT(Factor = allocVector(INTSXP, ncol));

  int
     *factor = INTEGER(Factor),
    ith = 0;

#define tablelength (1 << (sizeof(*first) * 8))
  pattern *patternTable[tablelength];
  for (Uint i=0; i < tablelength; i++) patternTable[i] = NULL;
  if (start > end) ERR("value of 'start' smaller than that of 'end'.");
  if ((int) end > pos[END]) ERR("value of 'end' outside the matrix 'CM'.");

  // where to base in hashing on
  // c_end > c_start + 1 : easy, at least one filled word (word = short Uint)
  // c_end == c_start : only the choice of two (lower or upper word)
  // else : 4 choices (lower and upper word of c_start and c_end
  if (c_end > c_start + 1) {
    addtofirst = sizeof(Uint) / sizeof(short Uint);
    assert(addtofirst == 2);
  } else {
    Uint shift,
      in_start = codesPerUnit - startremainder,
      in_end = endremainder + 1,
      half = codesPerUnit / 2;
    if (c_start == c_end) {
      firstpattern = SHORT_ONES;
      if (in_start < in_end) {
	addtofirst = BigEndian;
	shift = (half - in_start) * bitsPerCode;
	if (shift > 0)
	  firstpattern = ((firstpattern << shift) & SHORT_ONES) >> shift;
	assert(in_end >= half);
	shift = (codesPerUnit - in_end) * bitsPerCode;
	firstpattern = (firstpattern >> shift) << shift;
      } else {
	//
	assert(in_start >= half);
	addtofirst = 1 - BigEndian;
	shift = (codesPerUnit - in_start) * bitsPerCode;
	firstpattern = ((firstpattern << shift) & SHORT_ONES) >> shift;
	shift = (half - in_end) * bitsPerCode;
	if (shift > 0) firstpattern = (firstpattern >> shift) << shift;
      }
    } else {
      if (in_start >= half) addtofirst = BigEndian;
      else if (in_end >= half) addtofirst = 2 + 1 - BigEndian;
      else if (in_start >= in_end) {
	addtofirst = BigEndian;
	shift = (half - in_start) * bitsPerCode;
	firstpattern = ((firstpattern << shift)  & SHORT_ONES) >> shift;
      } else {
	addtofirst = 2 + 1 - BigEndian;
	shift = (half - in_end) * bitsPerCode;
	firstpattern = (firstpattern >> shift) << shift;
      }
    }
  }
  for (Uint nc=0; nc<ncol; nc++, cur_cm += nrow) {
    first  = ((short Uint*) (cur_cm + c_start) ) + addtofirst;
    //    assert(*first < tablelength);
    pattern **p = patternTable + (*first & firstpattern);
    while (*p != NULL && notequal(cur_cm, cm + (*p)->where * nrow,
				  c_start, c_end,startremainder,endremainder)) {
      p = &((*p)->next);
    }
    if (*p == NULL) {
      *p = (pattern *) MALLOC(sizeof(pattern));
      (*p)->where = nc;
      (*p)->ith = ith++;
      (*p)->count = 1;
      (*p)->next = NULL;
    } else ((*p)->count)++;
    factor[nc] = (*p)->ith + 1;
  }

  SEXP WhereVec, CountsVec;
  PROTECT(WhereVec = allocVector(INTSXP, ith));
  PROTECT(CountsVec = allocVector(INTSXP, ith));
  int *where = INTEGER(WhereVec),
    *count = INTEGER(CountsVec);
  for (Uint i=0; i < tablelength; i++) {
    pattern *p = patternTable[i];
    while (p != NULL) {
      where[p->ith] = p->where + 1;
      count[p->ith] = p->count;
      pattern *q = p->next;
      FREE(p);
      p = q;
    }
  }
  setAttrib(Factor, Where, WhereVec);
  setAttrib(Factor, Counts, CountsVec);
  UNPROTECT(3);
  return Factor;
}




SEXP colSumsEqualSNPs(SEXP  CM, SEXP Start, SEXP CV, SEXP Select) {
  SEXP Pos = getAttrib(CV, Position);
  int *pos = INTEGER(Pos);
  if (!pos[SNPxIND])
    ERR("'colSumsEqualSNPs' only defined for SNPxINDIVID matrices");
  Uint *cm_orig = (Uint*)  INTEGER(CM),
    *cv = (Uint*)  INTEGER(CV),
    ncol_cv = isMatrix(CV) ? ncols(CV) : 1,
    nrow_cv = isMatrix(CV) ? nrows(CV) : length(CV),
    len = length(Select),
    *select = len == 0 ? NULL : (Uint*) INTEGER(Select),
    nrow = nrows(CM),
    ncol = len == 0 ? ncols(CM) : len,
    size =  (Uint) (pos[END] - pos[START] + 1),
    start = INTEGER(Start)[0],
    end =  start + size - 1,
    c_start = start / codesPerUnit,
    c_end = end / codesPerUnit,
    startremainder = start % codesPerUnit,
    endremainder = end % codesPerUnit,
    startremainderBpC = startremainder * bitsPerCode,
    endshift = ((codesPerUnit - 1 - endremainder) * bitsPerCode);
  if (c_end > nrow)
    ERR3("'V' too long (debugging info: length=%d start=%d nrow=%d)",
	 endremainder, start, c_end);

  if (false) {
    for (Uint nr=0; nr<nrow; nr++) {
      for (Uint nc=0; nc<ncol; nc++) PRINTF("%8x ", cm_orig[nr + nc * nrow]);
      PRINTF(" %8x\n", cv[nr]);//
    }
    PRINTF("%d %d %d %d end=%d %d %d %d %d -- %d \n", startremainder, //
	   (Uint) pos[STARTREMAINDER],  endremainder, (Uint) pos[ENDREMAINDER],
	   end, codesPerUnit, start, size, pos[END], pos[START]);
  }

  if (startremainder != (Uint) pos[STARTREMAINDER] ||
      endremainder != (Uint) pos[ENDREMAINDER]) BUG;

  assert(CM != CV || ncol == ncol_cv);

  union {
    Uint d;
    unsigned short int x[2];
  } difference;
  // do not delete next comment: started work 20 March 2018 with Torsten, see
  // email to him
  // bool ANYvalue = GLOBAL.blocker.ANY_diff_value != 0.0 ||
  //    GLOBAL.blocker.ANY_allequal_value != 0.0;
  //  int SXP = ANYvalue ? REALSXP : INTSXP;
  int SXP = INTSXP;
  SEXP Ans;
  if (isMatrix(CV)) PROTECT(Ans = allocMatrix(SXP, ncol, ncol_cv));
  else PROTECT(Ans = allocVector(SXP, ncol));

  //do not delete next comment: started work 20 March 2018 with Torsten, see
  // email to him
  /*
  if (ANYvalue) {
    double *ans = REAL(Ans),
      sizeD = (double) size;
    for (Uint nc_cv = 0; nc_cv<ncol_cv; nc_cv++, cv += nrow_cv, ans += ncol_cv){
      Uint *cm = cm_orig;
      for (Uint nc = 0; nc < ncol; nc++, cm += nrow) {
	if (CM == CV && nc_cv > nc) { // quadratic & symmetric matrix !
	  ans[nc] = REAL(Ans)[nc * ncol + nc_cv];
	  continue;
	}
	Uint sum_unequal,
	  i = c_start,
	  j = 0;
	if (select != NULL) // selection given
	  cm = cm_orig + (nrow * (select[nc] - 1)); // ueberschreibt cm+=nrow
	difference.d = (cv[j] ^ cm[i]) << startremainderBpC;
	if (i == c_end) { // might happen for very small sizes
	  difference.d >>= startremainderBpC + endshift;
	  sum_unequal = unequal[difference.x[0]] + unequal[difference.x[1]];
	  ans[nc] = sizeD - sum_unequal;
	  continue;
	}

	sum_unequal = unequal[difference.x[0]] + unequal[difference.x[1]];
	for (i++, j++ ; i<c_end; i++, j++) {
	  difference.d = cv[j] ^ cm[i];
	  sum_unequal += unequal[difference.x[0]] + unequal[difference.x[1]];
	}
	difference.d = cv[j] ^ cm[i];
	difference.d >>= endshift;
	sum_unequal += unequal[difference.x[0]] + unequal[difference.x[1]];
	ans[nc] = sizeD - sum_unequal;
      }
    }
 } else
    */
  {
    int *ans = INTEGER(Ans);
    for (Uint nc_cv = 0; nc_cv<ncol_cv; nc_cv++, cv += nrow_cv, ans += ncol_cv){
      Uint *cm = cm_orig;
      for (Uint nc = 0; nc < ncol; nc++, cm += nrow) {
	if (CM == CV && nc_cv > nc) { // quadratic & symmetric matrix !
	  ans[nc] = INTEGER(Ans)[nc * ncol + nc_cv];
	  continue;
	}
	Uint sum_unequal,
	  i = c_start,
	  j = 0;
	if (select != NULL) // selection given
	  cm = cm_orig + (nrow * (select[nc] - 1)); // ueberschreibt cm+=nrow
	difference.d = (cv[j] ^ cm[i]) << startremainderBpC;
	if (i == c_end) { // might happen for very small sizes
	  difference.d >>= startremainderBpC + endshift;
	  sum_unequal = unequal[difference.x[0]] + unequal[difference.x[1]];
	  ans[nc] = size - sum_unequal;
	  continue;
	}

	sum_unequal = unequal[difference.x[0]] + unequal[difference.x[1]];
	for (i++, j++ ; i<c_end; i++, j++) {
	  difference.d = cv[j] ^ cm[i];
	  sum_unequal += unequal[difference.x[0]] + unequal[difference.x[1]];
	}
	difference.d = cv[j] ^ cm[i];
	difference.d >>= endshift;
	sum_unequal += unequal[difference.x[0]] + unequal[difference.x[1]];
	ans[nc] = size - sum_unequal;
      }
    }
  }

  UNPROTECT(1);
  return Ans;
}


#define INTERSECT(INT1, INT2, INTEGER1, INTEGER2)			\
  INT1 *a = INTEGER1(A);						\
  INT2 *b = INTEGER2(B);						\
  for (int i=0; i<lenA; i++) {						\
    INT1 aa = a[i];							\
    if (i > 0 && aa <= a[i-1]) { PRINTF("%d %d\n", (int) a[i-1], (int) aa); BUG;} \
    while (iB < lenB && b[iB] < aa) iB++;				\
    if (iB == lenB) break;						\
    if (aa == b[iB]) {							\
      intersct[n++] = aa;						\
      iB++;								\
    }									\
  }


SEXP intersect(SEXP A, SEXP B){
  //printf("%d %d %d %d\n", TYPEOF(A), TYPEOF(B), length(A), length(B));
  int n = 0,
    lenA = length(A),
    lenB = length(B),
    maxlen = lenA <= lenB ? lenA : lenB,
    iB = 0;
  if (lenA == 0 || lenB == 0) return allocVector(INTSXP, 0);
  int
    *intersct = (int*) MALLOC(sizeof(int) * maxlen);
  SEXP Ans;
  if (TYPEOF(A) == INTSXP) {
    if (TYPEOF(B) == INTSXP) { INTERSECT(int, int, INTEGER, INTEGER); }
    else { INTERSECT(int, double, INTEGER, REAL); }
  } else {
    if (TYPEOF(B) == INTSXP) { INTERSECT(double, int, REAL, INTEGER); }
    else { INTERSECT(double, double, REAL, REAL); }
  }
  PROTECT(Ans = allocVector(INTSXP, n));
  int *ans = INTEGER(Ans);
  for (int i=0; i<n; i++) ans[i] = intersct[i];
  FREE(intersct);
  UNPROTECT(1);
  return Ans;
}
