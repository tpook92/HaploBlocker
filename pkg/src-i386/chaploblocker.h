
#ifndef chaploblocker_H
#define chaploblocker_H 1

//
#include <R.h>
#include <Rinternals.h>

//extern "C" {
#ifdef __cplusplus
extern "C" {
#endif
  SEXP fixcoding(SEXP values);
  SEXP codeSNPs(SEXP M, SEXP Start, SEXP RedoCoding, SEXP SNPxIND);
  SEXP decodeSNPs(SEXP CM);
  SEXP factorSNPs(SEXP M, SEXP Start, SEXP End);
  SEXP colSumsEqualSNPs(SEXP  CM, SEXP start, SEXP CV, SEXP Select);
  SEXP intersect(SEXP A, SEXP B);

  void loadoptions();
  SEXP attachoptions();
  void detachoptions();
  //  SEXP copyoptions();

  
#ifdef __cplusplus
}
#endif

#endif /* _chaploblocker_H */

