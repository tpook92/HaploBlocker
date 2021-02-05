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
#include "intrinsics.h"
#include <Basic_utils.h>
#include <zzz_RandomFieldsUtils.h>
#include "chaploblocker.h"
#include "options.h"
//#include "def.h"
#include "xport_import.h"
#include "error.h"
#include "kleinkram.h"




const char * prefixlist[prefixN] = {"blocker", "blocker_messagesx"};


// IMPORTANT: all names of general must be at least 3 letters long !!!
const char *blocker[blockerN] = {"(unused)" };
const char *messages[messagesN] = 
  {"(none)"}; 
  

globalparam GLOBAL = {
blocker_START,
  messages_START
};


const char **all[prefixN] = {blocker};
int allN[prefixN] = {blockerN};


void setoptions(int i, int j, SEXP el, char name[200],
		  bool VARIABLE_IS_NOT_USED isList, bool local) {

 if (!local && parallel())
    ERR("'RFoptions' may not be set from a parallel process.");
      
  globalparam *options = WhichOptionList(local);
 
  switch(i) {
  case 0: {// blocker
    blocker_param *gp;
    gp = &(options->blocker);
    switch(j) {
    case 0: gp->ANY_diff_value = POS0NUM; break;
     default: BUG;
    }}
    break;
  case 1: {
    messages_param *gp = &(options->messages);
    switch(j) {
    case 0: gp->warn_dummy = LOGI; break;
    default: BUG;
    }}
    break;   
  default: BUG;
  }
}



void getoptions(SEXP sublist, int i, bool local) {
  int  k = 0;
  globalparam *options = WhichOptionList(local);
  switch(i) {
  case 0 : {
    blocker_param *p = &(options->blocker);
    ADD(ScalarReal(p->ANY_diff_value));
  }
    break;

    case 1 : {
   messages_param *p = &(options->messages);
   ADD(ScalarLogical(p->warn_dummy));
  }
    break;   
 
  default : BUG;
  }
    assert (i == prefixN - 1);
}


