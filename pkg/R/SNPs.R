
## M =  SNPx x animals
## CM is coded matrix
## V is a vector

START_C <- 0

fixcoding <- function(values)
  .Call(C_fixcoding, if (is.numeric(values)) as.integer(values) else values)

codeSNPs <- function(M, redo=is.character(M), SNPxINDIVID=TRUE)
  .Call(C_codeSNPs, M, 0L, as.logical(redo), SNPxINDIVID)

decodeSNPs <- function(CM) .Call(C_decodeSNPs, CM)

factorSNPs <- function(CM, start=1, end=attr(CM, "position")[2] + 1) {
  CM <- as.matrix(CM)
  if (is.null(attr(CM, "position")))
    CM <- .Call(C_codeSNPs, CM, 0L, FALSE, TRUE)
  .Call(C_factorSNPs, CM, as.integer(start - 1), as.integer(end - 1))
}

colSumsEqualSNPs <- function(CM, start = 1, V, select = NULL) {
  start <- as.integer(start - 1)
  if (is.null(attr(CM, "position")))
    CM <- .Call(C_codeSNPs, CM, start, FALSE,TRUE)
  pos <- attr(CM, "position")
  stopifnot(!is.null(pos))
  if (missing(V)) V <- CM
  else if (is.null(attr(V, "position")))
    V <- .Call(C_codeSNPs, V, start, FALSE, TRUE)
 .Call(C_colSumsEqualSNPs, CM, start, V, select)
}

intersect <- function(a, b) .Call(C_intersect, a, b)


if (FALSE) {

  Equal <- function(A, B) {
    attr(A, "position") <- attr(A, "position")[1:5]
#   Print(A, B)
    stopifnot(all.equal(A, B))
  }
  
  fixcoding <- function(values) {
    ans <- .Call(C_fixcoding, if (is.numeric(values)) as.integer(values) else values)
    Equal(ans, CHaploBlockerAlt::fixcoding(values))
    return(ans)
  }
  
  codeSNPs <- function(M, redo=is.character(M), SNPxINDIVID=TRUE) {
##    Print(list(C_codeSNPs, M, 0L, as.logical(redo), SNPxINDIVID))
    ans <- .Call(C_codeSNPs, M, 0L, as.logical(redo), SNPxINDIVID)
    Equal(ans, CHaploBlockerAlt::codeSNPs(M, redo))
    return(ans)
  }
  
  decodeSNPs <- function(CM) {
   ans <- .Call(C_decodeSNPs, CM)
     Equal(ans, CHaploBlockerAlt::decodeSNPs(CM))
   return(ans)
  }
  
  factorSNPs <- function(CM, start=1, end=attr(CM, "position")[2] + 1) {
    CM <- as.matrix(CM)
    if (is.null(attr(CM, "position")))
      CM <- .Call(C_codeSNPs, CM, 0L, FALSE, TRUE)
    ans <-  .Call(C_factorSNPs, CM, as.integer(start - 1), as.integer(end - 1))
     Equal(ans, CHaploBlockerAlt::factorSNPs(CM, start, end))
   return(ans)
  }

  colSumsEqualSNPs <- function(CM, start = 1, V, select = NULL) {
    xx <- CHaploBlockerAlt::colSumsEqualSNPs(CM, start, V, select)
    
##     Print(start)
    start <- as.integer(start - 1)
    if (is.null(attr(CM, "position")))
      CM <- .Call(C_codeSNPs, CM, start, FALSE,TRUE)
  
     if (missing(V)) V <- CM
     else if (is.null(attr(V, "position")))
      V <- .Call(C_codeSNPs, V, start, FALSE, TRUE)
      
    ##    Print("neu", CM, start, V, select);
  ##  Print(C_colSumsEqualSNPs, CM, start, V, select)
    ans <- .Call(C_colSumsEqualSNPs, CM, start, V, select)
    Equal(ans, xx)
    return(ans)
 }

  intersect <- function(a, b) {
   ans <-  .Call(C_intersect, a, b)
    Equal(ans, CHaploBlockerAlt::intersect(a,b))
   return(ans)
 }

}
