#include "trie.h"
#include <iostream>
#include "constEigen.h"

using namespace std;

int main() 
{	
	unsigned short* bacPools=(unsigned short*)malloc(sizeof(unsigned short)*BACS*LAYERS);
	readBacMappings(bacPools);
	Trie trie;
	//insert all BACs in the trie
	for (unsigned short i=0;i<BACS;i++)
		trie.insert(&bacPools[i*LAYERS], 7, i+1);
	
	//test the trie by doing some searches
	int i, j, len=9, numRes=0;
	unsigned short res[7];
	/*when the sig lenght is 10-11, 14 or 18-25 don't consider matches of length 6 because it can create false positives!!!*/
	unsigned short sig[] = {10,15,16,34,42,62,71,81,90};
	trie.search(sig, len, res, numRes);
	cout << "Input sig ";
	for (i=0;i<len;i++) 
		cout << sig[i] << " ";
	cout << endl;
	cout << "Bacs found " << endl;
	for (i=0;i<numRes;i++) {
		cout << "\t" << res[i] << ": ";
		for (j=0;j<LAYERS;j++)
			cout << bacPools[(res[i]-1)*LAYERS+j] << " ";
		cout << endl;
	}
	delete bacPools;
	return 0;
}
