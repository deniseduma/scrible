/***************************************************************************
 *   Copyright (C) 2011 by Marco Beccuti   *
 *   beccuti@di.unito.it   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

//!Disable (0) the reverse in Hash Table
#define REVERSE 0
//!if 1 the input file extension is .fasta otherwise .txt
#define FASTAEXFILE 0
//!Enable detailed output 
#define OUTPUT 0
//!Enable Pools output  per scripts
#define OUTPUTPOOLS 0
//!Enable Pools output  per deconvolution
#define OUTPUTPOOLSDEC 1
//! If 1  we enable counting output
#define COUNTING 0
//!if 1 we use in the hash table integer instead of pointer and  a char vector to encodes the pools bitvector
#define HASHINT  1
//! if 1 enable the hash analisys (YACC&LEX)
#define ANALYSIS 0
//! if 1 enable reduce the HASH DIMENSION (the frequency vector is reduced of dimension). E.g. 1 for RICE
#define REDUCEDHASH 1
//! if 1 enable the output for Denise deconvolution tool 
#define OUTPUTXDENISE 0
//! it filters k-mers that do not appear in more than F pools out of the hashtable 
#define FILTERKMERS 1
