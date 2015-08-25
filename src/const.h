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

#ifndef __CNF_H__
	#define __CNF_H__
	#include "conf.h"
#endif

//!k-mer size.
#define DIM 31 
//!number of pools.
#define POOLS 91
//!number of char needed to encoded all pools.
#define POOLS2C 12
//!max BAC 
#define LAYERS 7
#define BACS 2197
#define BACS_IN_POOL 169
#define POOLS_PER_LAYER 13

//!max buffer's size / read size.
#define MAXSIZE 250
#define SIZEMAX 1000000

//!max HASH BAC
#define HASHBACS 34999001

//!Hash function alpha coefficient.
#define A 48271
//!Hash function prime number
#define L 2147483647

//!Max read size
#define READSIZE 200
//!Read threshold
#define TREAD 91
//!Window threshold
#define TWINDOW 28
//!Window thresholds on pools
#define TW1 6
#define TW2 7
#define TW3 12
#define TW4 13
#define TW5 14
#define TW6 19
#define TW7 20
#define TW8 21

//!Pool frequency threshold
#define TCOR  0.05
//!Default value used for  null pointer
#define DEFAULTP  4294967295U

//!For statistic, frequencies cut.
#define CUT 100

//!Size of counting vector (used only for counting windows' frequencies in the pools).
#define NUMC 91 

//!Size of each cell in the counting vector
//#if REDUCEDHASH
//	#define MAXFREQ 255U //per rice
//#else
	#define MAXFREQ 65534
//#endif

//! Maximum number of overlappings
#define MAXOVERLAP 3

//!Size of each cell in BAC vector
#define MAXSHORT  65535U
//!Remove from hashtable all k-mers that do not appear in more than F pools
#define F 2
//!Remove from hash table all kmers that appear in more than MAXREP pools 
#define MAXREP 91
//!Remove from hash table all kmers that appear in less than MINREP pools 
#define MINREP 3
//!The min ration between hash table size and hash table capacity
#define MIN 3
//!The max ration between hash table size and hash table capacity
#define MAX 7
//!The number of prime hash sizes we consider when resizing the hash table up and down 
#define SIZES 34

//max BACs returned by a search in the trie
#define MAX_BACS 2197
//!Min number of pools to attempt search in the trie
#define LOW 4
#define LOW1 6
#define LOW2 6
//!Max number of pools to attempt search in the trie
#define HIGH 30
#define HIGH1 30
#define HIGH2 91
