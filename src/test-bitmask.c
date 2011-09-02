#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <bitmask.h>


int main(int argc, char* argv[]) {
#define LEN 127
	char * mask = bitmask_alloc(LEN);
	bitmask_clear_all(mask);
	int i, j;
	//printf("%p\n", mask);
	for(i = LEN - 1; i >= 0; i--) {
		if(bitmask_test_and_set(mask, i)) {
			abort();
		}
		for(j = 0; j < LEN; j++) {
			printf("test bit %d against %d\n", i, j);
			if(i!=j && bitmask_test(mask, j)) {
				abort();
			}
			if(i==j && !bitmask_test(mask, j)) {
				abort();
			}
		}
		printf("test bit %d\n", i);
		if(!bitmask_test_and_clear(mask, i)) {
			abort();
		}
	}

}
