/*
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 library for simulation of random fields

 Copyright (C) 2001 -- 2016 Martin Schlather,

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.
RO
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

//#include <stdio.h>
//#include <stdlib.h>
//#include <unistd.h>
//#include <string.h>


// ACHTUNG: Reihenfolge nicht aendern!
#include "chaploblocker.h"
#include "options.h"
//#include "def.h"
#include "xport_import.h"
#include <Basic_utils.h>
#include "error.h"
#include <zzz_RandomFieldsUtils.h>



CALL1(void, getErrorString, errorstring_type, errorstring)
CALL1(void, setErrorLoc, errorloc_type, errorloc)



const char * prefixlist[prefixN] =
  {"blocker"};


// IMPORTANT: all names of general must be at least 3 letters long !!!
const char *blocker[blockerN] =
  {"unused!", "unused!"};



int PL=C_PRINTLEVEL,
  CORES=INITCORES;
globalparam GLOBAL = {
blocker_START
};
utilsparam *GLOBAL_UTILS;


const char **all[prefixN] = {blocker};
int allN[prefixN] = {blockerN};


void setparameter(int i, int j, SEXP el, char name[200],
		  bool VARIABLE_IS_NOT_USED isList, int local) {
#ifdef DO_PARALLEL
  if (local != isGLOBAL) ERR1("Options specific to haploblocker, here '%s', can be set only via 'RFoptions' outside any parallel code.", name);
#endif  
  globalparam *options = &GLOBAL;
  switch(i) {
  case 0: {// blocker
    blocker_param *gp;
    gp = &(options->blocker);
    switch(j) {
    case 0: gp->ANY_diff_value = POS0NUM; break;
    case 1: gp->ANY_allequal_value = POS0NUM; break;
     default: BUG;
    }}
    break;
  default: BUG;
  }
}


void finalparameter(int VARIABLE_IS_NOT_USED local) {
  PL = GLOBAL_UTILS->basic.Cprintlevel;
  CORES = GLOBAL_UTILS->basic.cores;
}

void getparameter(SEXP sublist, int i, int VARIABLE_IS_NOT_USED local) {
  int  k = 0;
  globalparam *options = &GLOBAL;
  switch(i) {
  case 0 : {
    blocker_param *p = &(options->blocker);
    ADD(ScalarReal(p->ANY_diff_value));
    ADD(ScalarReal(p->ANY_allequal_value));
  }
    break;
  default : BUG;
  }
    assert (i == prefixN - 1);
}


void attachRFoptionsHaploBlocker() {
  includeXport();
  Ext_getUtilsParam(&GLOBAL_UTILS);
  //  GLOBAL_UTILS->solve.max_chol = 8192;
  //  GLOBAL_UTILS->solve.max_svd = 6555;
/*
  spam.min.n = as.integer(400) # 400
  spam.tol = 1e-8
  spam.min.p = 0.8 # 0.8
  spam.n = as.integer(1:500)
  spam.factor = as.integer(4294967)# prime, so that n * factor is still integer
  silent = TRUE
*/
  finalparameter(isGLOBAL);
  Ext_attachRFoptions(prefixlist, prefixN, all, allN,
		      setparameter, finalparameter, getparameter, NULL,
		      -10, false);
  finalparameter(isGLOBAL);
  
}

void detachRFoptionsHaploBlocker() {
  Ext_detachRFoptions(prefixlist, prefixN);
}

