


/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2016 Martin Schlather

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  
*/


#ifndef CHaploBlocker_options_H
#define CHaploBlocker_options_H 1 

#define blockerN 2
typedef struct haploblocker_param {
  double ANY_diff_value, ANY_allequal_value;
} blocker_param;

#define blocker_START { 0.0, 0.0 }


typedef struct globalparam{
  blocker_param blocker;
} globalparam;
extern globalparam GLOBAL;

#define prefixN 1
extern const char * prefixlist[prefixN], **all[prefixN];
extern int allN[prefixN];
void setparameter(int i, int j, SEXP el, char name[200], bool isList,
		  int local);
void getRFoptions(SEXP *sublist);
void finalparameter(int local);


#endif
