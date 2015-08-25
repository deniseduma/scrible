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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef _CIO_H_
	#define _CIO_H_
	#include "classIO.h"
#endif

#ifndef THR_H_
	#define THR_H_
	#include <pthread.h>
#endif

#ifndef __OPER_H__
	#define __OPER_H__
	#include "operation.h"
#endif


#ifndef __VCT_H__
	#define __VCT_H__
	#include <vector>
#endif

#if ((COUNTW || REPET)&&(ANALYSIS))
	#ifndef __PRS_H__
		#define __PRS_H__
		#include "parserYacc.h"
	#endif
#endif

#include <sys/time.h>
#include <sys/resource.h>

#ifndef __CON_H__
	#define __CON_H__
	#include "const.h"
#endif

#include "Timer.h"
#include <unistd.h>
//#include <thread>
//#include<atomic>
//#include<future>

using namespace IOFn;
using namespace Cl_HASH;

IOF<mybitsetx>* p_ref; 

struct thread_data* td_array;

unsigned char c[8]={ 0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80};
unsigned char r[8]={ 0xFE, 0xFD, 0xFB, 0xF7, 0xEF, 0xDF, 0xBF, 0x7F};

unsigned int hashSizes[SIZES]={5000549, 10000643, 20000003, 30000001, 40000003, 50000017, 60000011,70000027, 80000023, 90000049, 100000567, 150000001, 200000477, 250000013, 300000641, 350000041, 400000651, 450000007, 500000461, 550000001, 600000607, 650000011, 700000747, 750000007, 800000533, 850000013, 900000637, 950000017, 1000000531, 1250000027, 1500000713, 1750000027, 2000000579, 2147483647}; 
//unsigned int hashSizes[SIZES]={30000001, 40000003, 50000017, 60000011,70000027, 80000023, 90000049, 100000567, 150000001, 200000477, 250000013, 300000641, 350000041, 400000651, 450000007, 500000461, 550000001, 600000607, 650000011, 700000747, 750000007, 800000533, 850000013, 900000637, 950000017, 1000000531, 1250000027, 1500000713, 1750000027, 2000000579, 2147483647}; 

#if ANALYSIS
using namespace OPER;
extern int initParser(string filenamem,vector <class OPERATOR>& );
#endif


/*void executeFile(int i, int l, int u, IOF* p_ref) {
	p_ref->OverlapToKmersLight(i, l, u);

}*/

int main(int argc, char *argv[])
{
clock_t startGlobal,endGlobal;
startGlobal=clock();

cout<<"\n\n =========================================================\n";
cout<<"|	      	        HASHFILTER        	          |\n";
cout<<"|	      	   FILTER OVERLAPPINGS     	          |\n";
cout<<" =========================================================\n";
cout<<"\n If you find any bug, send an email to beccuti@di.unito.it\n";

if (argc!=7)
	{
	std::cerr<<"\n\nUSE: hashfilter <path/ReadFile> <path/MapPoolsToBacsFile> <path/OutputFile> <k-mer_size> <input_file_number> <thread>\n\n";
		exit(EXIT_FAILURE);
	}

unsigned int qgram=atoi(argv[4]);

if (qgram != DIM)
	{
	cout<<"Error: the size of constant DIM (const.h) is different than the kmer size"<<endl;
	}
unsigned int cores=atoi(argv[6]); 
unsigned int files=atoi(argv[5]); 

cout<<"\n\n___________________________________________________________\n\n\t\t\tInput parameters \n___________________________________________________________\n\n";
cout<<"\tkmer size: "<<qgram<<"\n";
cout<<"\tNumber of input files: "<<files;
if (FASTAEXFILE)
cout<<" with file extension \".fasta\"\n";
else
cout<<" with file extension \".txt\"\n";
cout<<"\tNumber of threads for mapping: "<<cores<<"\n";
cout<<"___________________________________________________________\n";

cout<<"\n\nSTART EXECUTION...\n"<<endl;

//string ReadFname=argv[1];
//string IOread=argv[3];
//string PoolsToBac=argv[2];
int who = RUSAGE_SELF;  
struct rusage usage;  
/*int ret;  
ret=getrusage(who,&usage); 
cout<<"\n\n\tTotal memory used: "<<usage.ru_maxrss<<"KB"<<endl;*/

//load BAC sigs into trie
Trie* trie=new Trie();
unsigned short* bacPools=(unsigned short*)malloc(sizeof(unsigned short)*BACS*LAYERS);
readBacMappings(bacPools);
for (unsigned short i=0;i<BACS;i++)
	trie->insert(&bacPools[i*LAYERS], 7, i+1);

//class IOF ref(ReadFname,PoolsToBac,IOread,qgram,files);
IOF<mybitsetx> ref(argv[1], argv[2], argv[3], qgram, files);

startGlobal=clock();

printf("size of size_t %lu\n", sizeof(std::size_t));
printf("size of short %lu\n", sizeof(short));
printf("size of unsigned short %lu\n", sizeof(unsigned short));
printf("size of unsigned char %lu\n", sizeof(unsigned char));
printf("size of char %lu\n", sizeof(char));

bool storedhash=ref.ReadB();
if (!storedhash)
	{
	Timer* pTimer=new Timer("Building k-mer hashtable");
	//for rice
	ref.ReadLight();
	ref.AllocHash();
	ref.ReadIntermediate();
	ref.Read();
	//for barley
	//ref.ReadBarley();
	//ref.RemoveRepKmers();
	ref.ReadBac();
	//ref.RetrieveBacs(trie);
	delete pTimer;
	}

getrusage(who,&usage); 
cout<<"\n\n\tTotal memory used: "<<usage.ru_maxrss/((double)1024*1024)<<"GB"<<endl;
cout<<"\n\nEND EXECUTION\n"<<endl;

endGlobal=clock();
cout<<"\n=========================== TIME ===========================\n\n";
cout<<"\tTime to create HASH TABLE: "<<((double)(endGlobal-startGlobal))/CLOCKS_PER_SEC<<"s."<<endl;
cout<<"\n============================================================\n\n";

if (!storedhash)
	{
	cout<<"\n---------------------- SAVING HASH TABLE ---------------------\n\n";
	ref.WriteB();
	cout<<"\n------------------------------------------------------------\n\n";
	cout << "\nDone saving. \n";
	}

//DEBUG
//while (1);

time_t time1,time2;
time1=time(NULL);
startGlobal=clock();

p_ref=&ref;

td_array = (thread_data*) malloc (cores*sizeof(thread_data));
unsigned int last=0,dim=files/cores;
if (cores>files) {
	cores=files;
	cout<<"Reseting the cores used to "<<cores<<endl;
}

omp_set_num_threads(cores);

for (unsigned int i=0;i<cores;i++) {
	td_array[i].l=last;
	td_array[i].u=last+dim-1;
	if ((td_array[i].u>files-1)||(i==cores-1))
		td_array[i].u=files-1;
	else
		last=last+dim;
	cout<<"\t Thread "<<i<<" works on the pools from "<<td_array[i].l<<" to "<<td_array[i].u<<endl;
}
cout<<endl;

//int i, id;
//#pragma omp parallel private(id, i) shared(td_array) 
//{
//id = omp_get_thread_num();
//p_ref->RemoveLowComplexity(id, td_array[id].l, td_array[id].u);
//p_ref->OverlapToKmers(id, td_array[id].l, td_array[id].u);
//p_ref->OverlapToKmersLight(id, td_array[id].l,td_array[id].u);
//}

/*******************C++11 threading********************************/

/*thread* threadarray[cores];
for (unsigned int i=0;i<cores;i++) {
	//async(launch::async,executeFile,i,td_array[i].l, td_array[i].u, p_ref);
	//threadarray[i] = new thread([&]() {
	//	p_ref->OverlapToKmersLight(i, td_array[i].l, td_array[i].u);
	//});

	threadarray[i] = new thread(executeFile, i, td_array[i].l, td_array[i].u, p_ref);
}

for (unsigned int i=0;i<cores;i++) 
	threadarray[i]->join();
	
cout << "hardware reported concurrency support = " << thread::hardware_concurrency() << endl;
*/

//remove all kmers from hash table which appear in less than MINREP pools or more than MAXREP pools 
//p_ref->RemoveRepKmers();

//check correctness of the hash table
//p_ref->CheckCorrectBacs();
p_ref->CheckCorrect();
//p_ref->print();

free(td_array);

time2=time(NULL);
endGlobal=clock();

cout<<"\n=========================== TIME ===========================\n\n";
cout<<"\tReal Time to split reads into kmers: "<<time2-time1<<"s."<<endl;
cout<<"\tClock Time to split reads into kmers: "<<((double)(endGlobal-startGlobal))/CLOCKS_PER_SEC<<"s."<<endl;
cout<<"\n============================================================\n\n";

return 0;

}


