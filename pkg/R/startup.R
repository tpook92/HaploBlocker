'#
  Authors
Torsten Pook, torsten.pook@uni-goettingen.de

Copyright (C) 2017 -- 2020  Torsten Pook

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
'#

.onAttach <- function(libname, pkgname) {

  #version_t <- as.numeric(as.Date("2020-05-12"))
  #today_t <- as.numeric(Sys.Date())

    packageStartupMessage(
"#############################################################
###################### HaploBlocker #########################
#############################################################
################ Version: 1.6.04 (22-03-2021) ###############
######## To update to the most recent stable version: #######
## devtools::install_github('tpook92/HaploBlocker', subdir='pkg') ##
#############################################################")


}
