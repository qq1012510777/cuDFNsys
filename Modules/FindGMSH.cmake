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

SET(GMSH_INCLUDE_SEARCH_PATH
    $ENV{HOME}/pkg/gmsh-4.8.4-source/MY_GMSH/include
)

SET(GMSH_LIBRARY_SEARCH_PATH
    $ENV{HOME}/pkg/gmsh-4.8.4-source/MY_GMSH/lib
)

FIND_PATH(GMSH_H    gmsh.h    ${GMSH_INCLUDE_SEARCH_PATH})
FIND_LIBRARY(GMSH_L    NAMES gmsh    PATHS ${GMSH_LIBRARY_SEARCH_PATH})

SET(GMSH_FOUND 1)

FOREACH(var GMSH_H GMSH_L)
    IF(NOT ${var})
    SET(GMSH_FOUND 0)
    ENDIF(NOT ${var})
ENDFOREACH(var)

IF(GMSH_FOUND)
  SET(GMSH_INCLUDE_DIR  ${GMSH_H})
  SET(GMSH_LIBRARIES    ${GMSH_L})
ENDIF(GMSH_FOUND)