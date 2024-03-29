#***************************************************************************#
# cuDFNsys - simulating flow and transport in 3D fracture networks          #
# Copyright (C) 2022, Tingchang YIN, Sergio GALINDO-TORRES                  #
#                                                                           #
# This program is free software: you can redistribute it and/or modify      #
# it under the terms of the GNU Affero General Public License as            #
# published by the Free Software Foundation, either version 3 of the        #
# License, or (at your option) any later version.                           #
#                                                                           #
# This program is distributed in the hope that it will be useful,           #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU Affero General Public License for more details.                       #
#                                                                           #
# You should have received a copy of the GNU Affero General Public License  #
# along with this program.  If not, see <https://www.gnu.org/licenses/>.    #
#***************************************************************************#


# SET(UMFPACK_INCLUDE_SEARCH_PATH
#   $ENV{HOME}/pkg/SuiteSparse-master/include
# )
# 
# SET(UMFPACK_LIBRARY_SEARCH_PATH
#   $ENV{HOME}/pkg/SuiteSparse-master/lib
#   $ENV{HOME}/pkg/SuiteSparse-master/lib
# )

FIND_PATH(UMFPACK_AMD_H      amd.h      ${UMFPACK_INCLUDE_SEARCH_PATH})
FIND_PATH(UMFPACK_UMFPACK_H  umfpack.h  ${UMFPACK_INCLUDE_SEARCH_PATH})

FIND_LIBRARY(UMFPACK_UMFPACK NAMES umfpack PATHS ${UMFPACK_LIBRARY_SEARCH_PATH})
FIND_LIBRARY(UMFPACK_AMD     NAMES amd     PATHS ${UMFPACK_LIBRARY_SEARCH_PATH})

SET(UMFPACK_FOUND 1)
FOREACH(var UMFPACK_UMFPACK_H UMFPACK_AMD_H UMFPACK_UMFPACK UMFPACK_AMD)
  IF(NOT ${var})
	SET(UMFPACK_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(UMFPACK_FOUND)
  SET(UMFPACK_INCLUDE_DIRS ${UMFPACK_AMD_H} ${UMFPACK_UMFPACK_H})
  SET(UMFPACK_LIBRARIES    ${UMFPACK_UMFPACK} ${UMFPACK_AMD})
ENDIF(UMFPACK_FOUND)
