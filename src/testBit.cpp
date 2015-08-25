#include <stdio.h>
#include <stdlib.h>
#include <bitset>
//#include <unistd.h>
//#include <sys/types.h>
//#include <sys/stat.h>
//#include <limits.h>
#include <inttypes.h>

#define CODE_A 0x00
#define CODE_C 0x01
#define CODE_G 0x02
#define CODE_T 0x03

using namespace std;

int main() {

int  bits=2;
uint64_t gram =0;           /* Succinct representation of a q-gram */
//uint64_t gram_reverse;   /* Succinct representation of the reverse of a q-gram */
uint64_t mask =((uint64_t)1 << (2*26))-1;

printf("sizeof(uint64_t) is %lu\n", sizeof(uint64_t));
printf("sizeof(unsigned long long int) is %lu\n", sizeof(unsigned long long));
bitset<2*1> b;
printf("sizeof(bitset) is %lu\n", sizeof(b));
gram = gram << bits | CODE_A;
gram = gram << bits | CODE_C;
gram = gram << bits | CODE_G;
gram = gram << bits | CODE_T;
printf("gram before applying the mask is %lu \n", gram);
printf("gram before applying the mask is %lx \n", gram);
gram = gram & mask;
printf("1 is %lu \n", (uint64_t)1);
printf("1 is %lx \n", (uint64_t)1);
printf("mask is %lu \n", mask);
printf("mask is %lx \n", mask);
printf("gram after applying the mask is %lu \n", gram);
printf("gram after applying the mask is %lx \n", gram);

return 0;
}
