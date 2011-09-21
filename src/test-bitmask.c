#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#define BITMASK_SHUFFLE
#include <bitmask.h>


int * buffer;
bitmask_t * mask;
int main(int argc, char* argv[]) {
#define LEN (65543)
#define NT (1024)
	buffer = calloc(sizeof(int), LEN);
	mask = bitmask_alloc(LEN);
	bitmask_clear_all(mask);
	int i;
	int js[NT];
	#pragma omp parallel for private(i)
	for(i = 0; i < NT; i++) {
		int j;
		for(j = 0; j < LEN; j++) {
			js[i] = j;
			while(bitmask_test_and_set(mask, j)) {
				continue;
			}
			buffer[j]++;
			buffer[j]++;
			buffer[j]--;
			bitmask_clear(mask, j);
		}
	}
	for(i = 0; i < LEN; i++) {
		if(buffer[i] != NT) {
			printf("failed on %d\n", i);
			abort();
		}
	}
	intptr_t tests[] = {
28250585,
28250789,
	};
	for(i = 0; i < sizeof(tests) / sizeof(intptr_t); i++) {
		intptr_t idx = bitmask_shuffle(tests[i]);
		printf("%lX %lX %lX ", tests[i], idx, bitmask_shuffle(idx));
		printf("ele= %ld", idx >> 5);
		printf("offset = %d\n", idx & 31);
	}
}
