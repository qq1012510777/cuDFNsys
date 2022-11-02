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

# do not use it
SET(MATLAB_INCLUDE_SEARCH_PATH
    /usr/local/MATLAB/R2020b/extern/include
)

SET(MATLAB_LIBRARY_SEARCH_PATH
    /usr/local/MATLAB/R2020b/bin/glnxa64
)

FIND_PATH(MATLAB_H    mat.h    ${MATLAB_INCLUDE_SEARCH_PATH})
FIND_LIBRARY(MATLAB_MAT_L    NAMES mat    PATHS ${MATLAB_LIBRARY_SEARCH_PATH})
FIND_LIBRARY(MATLAB_MX_L     NAMES mx     PATHS ${MATLAB_LIBRARY_SEARCH_PATH})
FIND_LIBRARY(MATLAB_MEX_L    NAMES mex    PATHS ${MATLAB_LIBRARY_SEARCH_PATH})
FIND_LIBRARY(MATLAB_ENG_L    NAMES eng    PATHS ${MATLAB_LIBRARY_SEARCH_PATH})

SET(MATLAB_FOUND 1)

FOREACH(var MATLAB_H MATLAB_MAT_L MATLAB_MX_L MATLAB_MEX_L MATLAB_ENG_L)
    IF(NOT ${var})
    SET(MATLAB_FOUND 0)
    ENDIF(NOT ${var})
ENDFOREACH(var)

IF(MATLAB_FOUND)
  SET(MATLAB_INCLUDE_DIR  ${MATLAB_H})
  # MESSAGE(${MATLAB_INCLUDE_DIR}"+++++++++++")
  SET(MATLAB_LIBRARIES    ${MATLAB_MAT_L} ${MATLAB_MX_L} ${MATLAB_MEX_L} ${MATLAB_ENG_L})
ENDIF(MATLAB_FOUND)