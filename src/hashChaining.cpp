#ifndef HASH_H_
	#define HASH_H_
	#include "hashChaining.h"
#endif

namespace Cl_HASH
{

//!It inserts all the elements from the first hash table into the second hash table. 
/*template <class A_Type>
void HashTable<A_Type>::insertAll(const HashTable<mybitset>& sh) {
	for (unsigned int i=0;i<sh.capacity();i++)
	{
		HashNode<mybitset>* tmpHashNode = sh.HashVec[i];
		while(tmpHashNode != NULL) {
		 	//mybitsetx mb(tmpHashNode->b.getBitset());
		 	A_Type mb(tmpHashNode->b);
			insert(mb);
			tmpHashNode = tmpHashNode->next;
    		};
	};
	assert(size()==sh.size());
};*/

/*template <class A_Type>
void HashTable<A_Type>::cmpTables(const HashTable<mybitset>& sh) {
	for  (unsigned int i=0;i<sh.capacity();i++) 
	{
		HashNode<mybitset>* tmpHashNode = sh.HashVec[i];
		while (tmpHashNode != NULL) {
			A_Type mb(tmpHashNode->b);
			assert(search(mb)==false);
			tmpHashNode = tmpHashNode->next;
		}
	}
};*/

}
