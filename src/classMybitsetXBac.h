
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


/*#ifndef __BIT_H__
	#define __BIT_H__
	#include <bitset>

#endif
#ifndef __MTH_H__
	#define __MTH_H__
	#include <math.h>
#endif

#ifndef __SET_H__
	#define __SET_H__
	#include <list>
#endif

#ifndef __CON_H__
	#define __CON_H__
	#include "const.h"
#endif*/

#ifndef __CMB_H__
	#define __CMB_H__
	#include "classMybitsetX.h"
#endif

#include <assert.h>

extern unsigned char c[8];
extern unsigned char r[8];
 
namespace MYBIT{

using namespace std;

class mybitsetxbac: public mybitsetx {
	
protected: 
//unsigned short numBacs;
//!Corresponding Bacs
bitset<MAX_BACS> bacs;

public:

//!Empty Constructor.
mybitsetxbac(): mybitsetx()/*, numBacs(0)*/ {};

//!Constructor
mybitsetxbac(const bitset<2*DIM>& b): mybitsetx(b)/*, numBacs(0)*/ {}

//!Constructor
mybitsetxbac(const mybitset& b): mybitsetx(b)/*, numBacs(0)*/ {}

//!Constructor
mybitsetxbac(const mybitsetx& b): mybitsetx(b) {
	//cc=b.cc;
	
	//for (unsigned int i=0; i<POOLS; i++) 
	//	cw[i]=b.cw[i];

}

//!Copy constructor
mybitsetxbac(const mybitsetxbac& b) {
	cc=b.cc;
	
	for (unsigned int i=0; i<POOLS; i++) 
		cw[i]=b.cw[i];
	
	//numBacs=b.numBacs;
	bacs=b.bacs;
}

//!Assignmenet operator
mybitsetxbac& operator=(const mybitsetxbac& b) {
	if (this!=&b) {
		cc=b.cc;
		
		for (int i=0; i<POOLS; i++) 
			cw[i]=b.cw[i];
		
		//numBacs=b.numBacs;
		bacs=b.bacs;
	}
	return *this;
}

~mybitsetxbac() {}

//inline unsigned short getBac(unsigned short i) {
//	return bacs[i];
//}

inline unsigned short getNumBacs() const {
	unsigned short numBacs=0;
	for (unsigned short i=0;i<BACS;i++) {
		if (bacs[i]==1)
			numBacs++;
	}
	return numBacs;
}


inline unsigned short getBacs(unsigned short res[]) const {
	unsigned short numBacs=0;
	for (unsigned short i=0;i<BACS;i++) {
		if (bacs[i]==1)
			res[numBacs++]=(i+1);
	}
	return numBacs;
}

inline void setBacs(const unsigned short input[], unsigned short numBacs) {
	//this->numBacs=numBacs;
	for (unsigned short i=0;i<numBacs;i++) {
		assert(input[i]>0&&input[i]<=BACS);
		bacs.set(input[i] - 1);
	}	
}

//!Copy operator for pools' bitvector.
inline void copy(mybitsetxbac& p){
	for (unsigned short jj=0; jj<POOLS; jj++)
		p.cw[jj]=cw[jj];
		
		//p.numBacs=numBacs;
		p.bacs=bacs;
};


void print(ostream& out) const {
	
	mybitset::print(out); 
	out << "%" << endl;
	for (unsigned short i=0;i<POOLS;i++) {
		if (getPool(i)==1)
			out<<i<<"("<<getPoolFreq(i)<<"), ";
	};
	out << endl << "%%" << endl;
	for (unsigned short i=0;i<BACS;i++) {
		if (bacs[i]) 
			out<<(i + 1)<<", ";
	};
	out << endl;
};

};
};
