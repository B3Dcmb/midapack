/**
@file midapack.h, version 1.2b, November 2012
@brief Header file with main definitions and declarations for the Midapack library
@author Pierre Cargemel, Frederic Dauvergne, Maude Le Jeune, Antoine Rogier, Radek Stompor
**
** Project:  Midapack library, ANR MIDAS'09
** Purpose:  Provide algebra tools suitable for Cosmic Microwave Background (CMB)
**           data analysis.
**
***************************************************************************
@note Copyright (c) 2010-2012 APC CNRS Université Paris Diderot
@note
@note This program is free software; you can redistribute it and/or modify it under the terms
@note of the GNU Lesser General Public License as published by the Free Software Foundation;
@note either version 3 of the License, or (at your option) any later version. This program is
@note distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
@note the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
@note Lesser General Public License for more details.
@note
@note You should have received a copy of the GNU Lesser General Public License along with this
@note program; if not, see http://www.gnu.org/licenses/lgpl.html
@note For more information about ANR MIDAS'09 project see :
@note http://www.apc.univ-paris7.fr/APC_CS/Recherche/Adamis/MIDAS09/index.html
@note ACKNOWLEDGMENT: This work has been supported in part by the French National Research
@note Agency (ANR) through COSINUS program (project MIDAS no. ANR-09-COSI-009).
***************************************************************************
**
**
*/


#include <fftw3.h>
#include <mpi.h>

//=========================================================================
//=======================  Toeplitz algebra module ========================
//=========================================================================
#include "../../src/toeplitz/toeplitz.h"



//=========================================================================
//============================  Mapmat module =============================
//=========================================================================

//============================ mapmat =============================
#include "../../src/mapmat/mapmat.h"
//============================ mapmat coarse ======================
#include "../../src/mapmat/mapmatc.h"
