#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include<unistd.h>
#include <bitset>
#include <inttypes.h>
#include <cstring>
#include <assert.h>


#define SIZE1 1000000000
#define SIZE2 2000000000

using namespace std;

void allocBigMem(int size) {
	uint64_t* bigArray;
	cout << "[allocBigMem]Before allocating array" << endl;
	//sleep(30000);
	bigArray = (uint64_t*)malloc(size*sizeof(uint64_t));
	memset(bigArray, 0, size*sizeof(uint64_t));
	cout << "[allocBigMem]After allocating array" << endl;
	sleep(30);
	//free(bigArray);
}

int main() {
	cout << "" << endl;
	allocBigMem(SIZE1);
	cout << "After returning from 1st function call" << endl;
	sleep(30);
	cout << "Entering the function again" << endl;
	allocBigMem(SIZE2);
	cout << "After returning from 2nd function call" << endl;
	sleep(30);
	
	cout <<"Done." << endl;
	
	return 0;
}
