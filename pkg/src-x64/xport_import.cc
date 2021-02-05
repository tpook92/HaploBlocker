/*
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 Copyright (C) 2015 -- 2019 Martin Schlather, 

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

#include <R_ext/Rdynload.h>
#include <General_utils.h>
#include "intrinsics.h"
#include "options.h"
#include "xport_import.h"
#include "kleinkram.h"
#include "chaploblocker.h"



#define importfrom "RandomFieldsUtils"

#ifdef CALL
#undef CALL
#endif
#define CALL(what) what##_type Ext_##what = NULL
UTILSCALLS;

#undef CALL
#define CALL(what) Ext_##what = (what##_type) R_GetCCallable(importfrom, #what)
void includeXport() {
  UTILSCALLS;
} // export C

bool ToFalse[1] = { false };

#define PLoffset -10
int PL = C_PRINTLEVEL,
  CORES = INITCORES; //  return;  TO DO: replace by KEYT->global_utils

SEXP Information = NULL,
  Coding = NULL;
bool HAS_CUDA = false;


utilsparam *GLOBAL_UTILS;
KEY_type *PIDKEY[PIDMODULUS];
int parentpid=0;
bool parallel() {
  int mypid;
  //  printf("pid\n");
  Ext_pid(&mypid);
  //  printf("pid = %d %d\n", mypid, parentpid);
  return mypid != parentpid;
}



void globalparam_NULL(KEY_type *KT, bool copy_messages) {
  //  printf("%d %d\n", blockerN, messagesN);
  assert(blockerN==6 && messagesN == 1);
  messages_param m;
  if (!copy_messages)
    MEMCOPY(&m, &(KT->global.messages), sizeof(messages_param));

  MEMCOPY(&(KT->global), &GLOBAL, sizeof(globalparam));
  // pointer auf NULL setzten
  if (!copy_messages)
    MEMCOPY(&(KT->global.messages), &m, sizeof(messages_param));

 
  Ext_utilsparam_NULL(&(KT->global_utils));
}

void globalparam_NULL(KEY_type *KT) {
  globalparam_NULL(KT, true);
}

void globalparam_DELETE(KEY_type *KT) {
   // pointer loeschen
  blocker_param *gp = &(KT->global.blocker);
  Ext_utilsparam_DELETE(&(KT->global_utils));
}


void KEY_type_NULL(KEY_type *KT) {
  // ACHTUNG!! setzt nur die uninteressanten zurueck. Hier also gar ncihts.
  KT->next = NULL; // braucht es eigentlich nicht
  globalparam_NULL(KT);
}

void KEY_type_DELETE(KEY_type **S) {
  KEY_type *KT = *S;
  globalparam_DELETE(KT);
   
  UNCONDFREE(*S);
}

KEY_type *KEYT() {
  int mypid;
  Ext_pid(&mypid);
  //  printf("entering KEYT %d %d \n", mypid, parentpid);
  // for (int i=0; i<PIDMODULUS; i++) printf("%ld ", PIDKEY[i]);
  KEY_type *p = PIDKEY[mypid % PIDMODULUS];
  //  printf("%d %d %ld\n", mypid,  PIDMODULUS, p);
  if (p == NULL) {
    KEY_type *neu = (KEY_type *) XCALLOC(1, sizeof(KEY_type));
    assert(neu != NULL);
    PIDKEY[mypid % PIDMODULUS] = neu;
    neu->visitingpid = mypid;    
    if (PIDKEY[mypid % PIDMODULUS] != neu) { // another process had the
      //                                        same idea
      FREE(neu);
      return KEYT(); // ... and try again
    }
    neu->pid = mypid;
    //    printf("neu %d %d\n", mypid);
    neu->visitingpid = 0;
    neu->ok = true;
    if (PIDKEY[mypid % PIDMODULUS] != neu) BUG;
    KEY_type_NULL(neu);    
    if (GLOBAL_UTILS->basic.warn_parallel && mypid == parentpid) {
      PRINTF("Do not forget to run 'RFoptions(storing=FALSE)' after each call of a parallel command (e.g. from packages 'parallel') that calls a function in 'RandomFields'. (OMP within RandomFields is not affected.) This message can be suppressed by 'RFoptions(warn_parallel=FALSE)'.\n"); // ok
    }
   return neu;
  }
  while (p->pid != mypid && p->next != NULL) {
    //    printf("pp = %d\n", p->pid);
    p = p->next;
  }
  //  printf("pp m = %d %d\n", p->pid, mypid);
  if (p->pid != mypid) {
    if (!p->ok || p->visitingpid != 0) {
      if (PL >= PL_ERRORS) {
	PRINTF("pid collision %d %d\n",  p->ok, p->visitingpid);
      }
      //
      BUG;
      return KEYT();
    }
    p->visitingpid = mypid;
    p->ok = false;
    if (p->visitingpid != mypid || p->ok) {
      return KEYT();
    }
    KEY_type *neu = (KEY_type *) XCALLOC(1, sizeof(KEY_type));
    neu->pid = mypid;
    if (!p->ok && p->visitingpid == mypid) {
      p->next = neu;
      p->visitingpid = 0;
      p->ok = true;      
      return neu;
    }
    FREE(neu);
    p->visitingpid = 0;
    p->ok = true;
    KEY_type_NULL(neu); 
   return KEYT();
  }
  return p;
}



SEXP copyoptions() {
  KEY_type *KT = KEYT();
  globalparam_NULL(KT, false);
  return R_NilValue;
}

SEXP setlocalRFutils(SEXP seed, SEXP printlevel) {
  KEY_type *KT = KEYT();
  assert(KT != NULL);
  utilsparam *global_utils = &(KT->global_utils);
  assert(global_utils != NULL);
  if (length(seed) > 0)
    global_utils->basic.seed = Integer(seed, (char *) "seed", 0);
  if (length(printlevel) > 0) {
    PL = global_utils->basic.Rprintlevel =
      Integer(printlevel, (char *) "printlevel", 0);
    global_utils->basic.Cprintlevel = global_utils->basic.Rprintlevel +PLoffset;
  }
  return R_NilValue;
}

void finalizeoptions() {
  utilsparam *global_utils = GLOBAL_UTILS;
  PL = global_utils->basic.Cprintlevel - PLoffset;
  CORES = global_utils->basic.cores;
}

void setoptions(int i, int j, SEXP el, char name[LEN_OPTIONNAME], bool isList, bool local);
void getoptions(SEXP sublist, int i, bool local);


void initHaploBlocker();
void loadoptions() {  
  //assert(Haplo == 30);
  for (int i=0; i<PIDMODULUS; i++) PIDKEY[i] = NULL;
  includeXport();
  Ext_pid(&parentpid);
  Ext_getUtilsParam(&GLOBAL_UTILS);
  // utilsparam *global_utils = GLOBAL_UTILS; // OK
  //  global_utils->solve.max_chol = 8192;
  //  global_utils->solve.max_svd = 6555; 
/*
  spam.min.n = as.integer(400) # 400
  spam.tol = 1e-8
  spam.min.p = 0.8 # 0.8
  spam.n = as.integer(1:500)
  spam.factor = as.integer(4294967)# prime, so that n * factor is still integer
  silent = T R U E
*/
  Ext_attachRFoptions(prefixlist, prefixN, all, allN,
		      setoptions, finalizeoptions, getoptions,
		      NULL, -10, false);
  finalizeoptions();  
  initHaploBlocker();
}




#define NEED_AVX2 true
#define NEED_AVX false
#define NEED_SSSE3  true
#define NEED_SSE2 true
#define NEED_SSE false
SEXP attachoptions() { // no print commands!!!
  #ifdef SCHLATHERS_MACHINE
  PRINTF("floating point double precision: %s\n",
#if defined DO_FLOAT
	 "no"
#else
	 "yes"
#endif
	 );
#endif
  ReturnAttachMessage(HaploBlocker, true);
}


globalparam *WhichOptionList(bool local) {  
  if (local) {
    KEY_type *KT = KEYT();
    if (KT == NULL) BUG;
    //    printf("wol = %ld\n", &(KT->global));
    return &(KT->global);
  }
  return &GLOBAL;
}

void PIDKEY_DELETE() {
  for (int kn=0; kn<PIDMODULUS; kn++) {
    KEY_type *KT = PIDKEY[kn];
    while (KT != NULL) {
      KEY_type *q = KT;
      KT = KT->next;
      KEY_type_DELETE(&q);
    }
    PIDKEY[kn] = NULL;
  }
}


void detachoptions() {
  PIDKEY_DELETE();
  Ext_detachRFoptions(prefixlist, prefixN);
  //  freeGlobals();
}


