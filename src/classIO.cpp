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

#ifndef _CIO_H_
	#define _CIO_H_
	#include "classIO.h"
#endif
#include<unistd.h>

namespace IOFn {

template <class A_Type>
void IOF<A_Type>::compareTables() {
	using namespace general;
	//cmpTables(shOld, sh);
};

template <class A_Type>
void IOF<A_Type>::Write(){
ofstream pFile;

time_t time_1,time_4;


time(&time_1);


pFile.open(OFname.c_str(),ofstream::out);
if (!pFile)
	{
	cerr << "\n*****Error opening ouput file *****" << endl;
	exit(EXIT_FAILURE);
	}

//set<class mybitset>::iterator iter=sh.begin();
pFile<<"W:"<<kmerSize<<" |W|:"<<sh.size()<<" R:"<<numReads<<endl;
pFile.close();
time(&time_4);
cout<<"Time to format input REFERENCE: "<<(time_4-time_1)<<endl;
};


template <class A_Type>
inline void IOF<A_Type>::convertN(char tmp[], const int& num)
{
char inv[MAXSIZE];
string kk;
for (int i=0;i<num;i++)
		{
		switch (tmp[i])
			{
			case 'A':
				inv[num-i-2]='T';
				kk+='T';
			break;
			case 'C':
				inv[num-i-2]='G';
				kk+='G';
			break;
			case 'G':
				inv[num-i-2]='C';
				kk+='C';
			break;
			case 'T':  
				inv[num-i-2]='A';
				kk+='A';
			break;
			case '\0':
 				inv[num-1]='\0';
			break;
			default:
				inv[num-i-2]='N';	
			}
		}
strcpy(tmp,inv);
};


/*void IOF::Average(vector<class OPERATOR>::iterator it,int i){

ostringstream of;
of<<OFname<<i;
ofstream out(of.str().c_str(),ofstream::out);
if (!out)
	{
	cerr << "\n*****Error opening ouput file *****" << endl;
	exit(EXIT_FAILURE);
	}

sh.Average(out,it,count);

// #if REPET
// //sh.CountR(out,count);
// sh.Corretion(out,count);
// #endif
cout<<"\t Output  operation "<<i<<" save in "<<of.str().c_str()<<endl;
out.close();
}*/

/*void IOF::Count(int i){
ostringstream of;
of<<OFname<<i;
ofstream out(of.str().c_str(),ofstream::out);
if (!out)
	{
	cerr << "\n*****Error opening ouput file *****" << endl;
	exit(EXIT_FAILURE);
	}

sh.Count(out);
cout<<"\t Output  operation "<<i<<" save in "<<of.str().c_str()<<endl;
out.close();
}*/


/*void IOF::Search(int l,int u){

using namespace general;

HashTable<mybitset> rd;
rd.updateHash(READSIZE-kmerSize);

char buffer[MAXSIZE],name[MAXSIZE],name2[MAXSIZE];

ifstream in;
ofstream out;

int num=0;
clock_t startGlobal,endGlobal;
mybitsetBAC p;

char delimC[] = ".";
Parser parser;
parser.update(delimC,ReadFname);
string fileExt="";
if (parser.size()>1)
	{
	fileExt=parser.get(1);
	}
for (int i=l;i<=u;i++)
{
startGlobal=clock();
ostringstream of,of1;

#if FASTAEXFILE
of<<parser.get(0)<<i<<"."<<fileExt;
#else
of<<parser.get(0)<<i<<"."<<fileExt;
#endif
//#pragma omp parallel
//{
//int i = omp_get_thread_num();
cout << "[Search] " << fileExt << " in file " << of.str() << endl;
//}

in.open(of.str().c_str(),ofstream::in);
if(!in) 
	{
	cerr << "\n*****Error opening input file *****" << endl;
	exit(EXIT_FAILURE);
	}
of1<<OFname<<i;
out.open(of1.str().c_str(),ofstream::out);
if (!out)
	{
	cerr << "\n*****Error opening ouput file *****" << endl;
	exit(EXIT_FAILURE);
	}

int dec=0;
int pdec=0;
while (!in.eof())
	{
	buffer[0]='\0';
	in.getline(buffer,MAXSIZE);
	num=in.gcount();
	if (buffer[num-1]!='\0')
		{
		buffer[num]='\0';
		num++;
		}
	if ((num>=kmerSize+1)&&((buffer[0]=='A')||(buffer[0]=='G')||(buffer[0]=='C')||(buffer[0]=='T')||(buffer[0]=='N'))) {
			p=search(rd,buffer,num,out);
	#if OUTPUTPOOLSDEC
			out<<name<<"\n";
			out<<buffer<<"\n";
	#else
		#if  OUTPUTPOOLS
			out<<name<<":$";
		#endif
	#endif
	} else
		if ((num>0)&&((buffer[0]=='@')||(buffer[0]=='>'))) {
			strcpy(name,buffer);
			p.clear();
		}	
	}

in.close();
out.close();
endGlobal=clock();
}//end for
};*/


/*mybitsetBAC  IOF::search(class HashTable<mybitset>& rd,const char buffer[],const int& num,ofstream& out )
{
mybitsetBAC p;
class mybitset b,d;
for(int k=0; k<=num-(kmerSize)-1; k++)
	{
	int z=0;
	int j=DIM*2-1;
	bool findN=false;
	for (int c = 0; c<kmerSize; c++) 
		{
              	switch (buffer[c+k])
			{
			case 'A':
				b.set(z,0);
				z++;
				b.set(z,0);
				z++;
				d.set(j,1);
				j--;
				d.set(j,1);
				j--;
			break;
			case 'C':
				b.set(z,0);
				z++;
				b.set(z,1);
				z++;
				d.set(j,0);
				j--;
				d.set(j,1);
				j--;
			break;
			case 'G':
				b.set(z,1);
				z++;
				b.set(z,0);
				z++;
				d.set(j,1);
				j--;
				d.set(j,0);
				j--;
			break;
			case 'T':
				b.set(z,1);
				z++;
				b.set(z,1);
				z++;
				d.set(j,0);
				j--;
				d.set(j,0);
				j--;
			break;
			default:
				findN=true;
				c=kmerSize;
			}
		}
	if (!findN)
		{
		z=0;	
		if (d<b)
			{
			rd.insert(d);
			}
		else
			{
			rd.insert(b);
			}
		}
	}
	if (rd.size()!=0)
 		p=rd.search(sh,shBAC);
#if OUTPUT
	out<<p<<endl;
#endif
return p;
}*/

};
