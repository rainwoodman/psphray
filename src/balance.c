#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <printf.h>
#include <mpi.h>
#include <array.h>
#include <messages.h>
#include "mortonkey.h"
#include "psystem.h"

#define N 32
#define BITS 8
#define THRESH 1

static void pos_hash_key(const void *p, void * k) {
	morton_key_t key = morton_key(p);
	memcpy(k, &key, sizeof(morton_key_t));
}

typedef struct {
	uint64_t key;
	int order;
	short int child_length;
	short int incomplete; /* if a node is incomplete, we can't return the ray tracing result as one segment, because some of the particles in this square are hosted on other processors. */
	union {
		intptr_t child[8];
		struct {
			intptr_t first_par;
			intptr_t npar;
		};
	};
	intptr_t parent;
} FNode;

typedef struct {
	FNode * pool;
	size_t size;
	size_t used;
} Framework;


Framework fwork;
MPI_Datatype MPI_PARTICLE;

size_t nearest_smaller_value(void* array, intptr_t output[], size_t nmemb, size_t bsize, int (*cmp)(const void* p1, const void * p2), int dir);

int main(int argc, char* argv[]) {
	PRINTF_BINARY('b');
	MPI_Init(&argc, &argv);
	morton_key_init(BITS);
	MPI_Type_contiguous(sizeof(Particle), MPI_BYTE, &MPI_PARTICLE);
	MPI_Type_commit(&MPI_PARTICLE);
	mockdata();

	int rank, size;
	intptr_t i;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	printf("rank = %d, psys.used/psys.size= %ld/%ld\n", rank, psys.used, psys.size);
	char * fname = NULL;
	asprintf(&fname, "balance-orig-%03d", rank);
	FILE * fp = fopen(fname, "w");
	for(i = 0; i < psys.used; i++) {
		fprintf(fp, "%08lo %g %g %g\n", morton_key(PSYS(pos,i)),
			PSYS(pos,i)[0],
			PSYS(pos,i)[1],
			PSYS(pos,i)[2]);
	}
	fclose(fp);

	balance();
	build_framework();


	asprintf(&fname, "balance-out-%03d", rank);
	fp = fopen(fname, "w");
	for(i = 0; i < psys.used; i++) {
		fprintf(fp, "%08lo %g %g %g\n", morton_key(PSYS(pos, i)),
			PSYS(pos, i)[0],
			PSYS(pos, i)[1],
			PSYS(pos, i)[2]);
	}
	fclose(fp);
	
	printf("rank = %d, fwork.used/fwork.size = %ld/%ld\n", rank, fwork.used, fwork.size);
	asprintf(&fname, "balance-fwork-%03d", rank);
	fp = fopen(fname, "w");
	for(i = 0; i < fwork.used; i++) {
		int j;
		fprintf(fp, "%04ld %08lo %02d %04ld %1d ", i, fwork.pool[i].key, fwork.pool[i].order, fwork.pool[i].parent, fwork.pool[i].child_length);
		if(fwork.pool[i].child_length ==0) {
			for(j = 0; j < fwork.pool[i].npar; j++) {
				fprintf(fp, "%08lo ", morton_key(psys.pool[j + fwork.pool[i].first_par].pos));
			}
		} else {
			for(j = 0; j < fwork.pool[i].child_length && j < 8; j++) {
				fprintf(fp, "%04ld ", fwork.pool[i].child[j]);
			}
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	MPI_Finalize();
	return 0;
}


int mockdata() {
	intptr_t i;
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	srandom(rank * 100033);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	psys.size = N *1.2;
	psys.used = N * (0.5 * random() / RAND_MAX + 0.3);
	psys.pool = malloc(sizeof(Particle) * psys.size);
	for(i = 0; i < N; i++) {
		int d;
		for(d = 0; d < 3; d++) {
			PSYS(pos, i)[d] = fmin((float) random() / RAND_MAX, 1.0);
		}
	}
	intptr_t *arg = malloc(sizeof(intptr_t) * psys.used);
	argsort(psys.pool, psys.used, sizeof(Particle), arg, morton_key_t_compare, particle_hash_key, sizeof(morton_key_t));
	argpermute(psys.pool, NULL, psys.used, sizeof(Particle), arg);
	free(arg);
	return 1;
}

int MPI_Gatherv2(void *sendbuf, int sendcount, MPI_Datatype type,
            void **recvbuf, int * recvcount, int root, MPI_Comm comm) {
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	int recvcounts[size];
	int displs[size];

	MPI_Gather(&sendcount, 1, MPI_INT, recvcounts, 1, MPI_INT, root, comm);

	displs[0] = 0;
	int i;
	for(i = 1; i < size; i++) {
		displs[i] = displs[i - 1] + recvcounts[i - 1];
	}

	if(rank == root) {
		int bsize = 0;
		MPI_Type_size(type, &bsize);
		*recvcount = displs[size - 1] + recvcounts[size - 1];
		*recvbuf = malloc(bsize * *recvcount);
	} else {
		*recvbuf = NULL;
		*recvcount = 0;
	}

	int rt = MPI_Gatherv(sendbuf, sendcount, type, *recvbuf, recvcounts, displs, type, root, comm);
	return rt;
}

int balance() {
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	intptr_t i;
	intptr_t j;

	size_t gnpar = 0;
	MPI_Allreduce(&psys.used, &gnpar, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	
/* step 1: find the percentiles from iteration */
	morton_key_t pctl_min[size];
	morton_key_t pctl_max[size];
	morton_key_t lmax = 0;
	morton_key_t lmin = 0;
	maxmin(psys.pool, psys.used, sizeof(Particle), &lmax, &lmin, morton_key_t_compare, particle_hash_key, sizeof(morton_key_t));
	morton_key_t gmax = 0;
	morton_key_t gmin = 0;
	MPI_Allreduce(&lmax, &gmax, 1, MPI_TYPE_MORTON_KEY, MPI_OP_MORTON_KEY_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(&lmin, &gmin, 1, MPI_TYPE_MORTON_KEY, MPI_OP_MORTON_KEY_MIN, MPI_COMM_WORLD);

	for(j = 0; j < size; j++) {
		pctl_max[j] = gmax;
		pctl_min[j] = gmin;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	intptr_t lcount[size];
	morton_key_t pctl_mid[size];

	while(1) {
		for(j = 0; j < size; j++) {
			pctl_mid[j] = (pctl_max[j] + pctl_min[j]) / 2;
		}
		for(j = 0; j < size; j++) {
			if(j > 0 && pctl_mid[j] == pctl_mid[j -1]) {
				/* omit some countings on the first few rounds because all pctls start the same */
				lcount[j] = lcount[j - 1];
				continue;
			}
			lcount[j] = searchsorted(pctl_mid + j, psys.pool, psys.used, sizeof(Particle), SEARCH_DIR_LEFT, morton_key_t_compare, particle_hash_key, sizeof(morton_key_t));
		}
		
		intptr_t gcount[size];
		MPI_Allreduce(lcount, gcount, size, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
		int bad = size;
		/* each chunk has about gnpar / size pars, and pctl_mid[j] is the end of the chunk */
		for(j = 0; j < size; j++) {
			if(gcount[j] > (j + 1) * gnpar / size) {
				if(pctl_max[j] == pctl_mid[j]) {
					bad--;
				}
				pctl_max[j] = pctl_mid[j];
			} else
			if(gcount[j] < (j + 1) * gnpar / size) {
				if(pctl_min[j] == pctl_mid[j]) {
					bad--;
				}
				pctl_min[j] = pctl_mid[j];
			} else
				bad --;
		}
		if(bad == 0) break;
	}

	for(i = 0; i < size; i++) {
		MESSAGE("rank = %d, pctl_mid[%d] = %08lo lcount[%d] = %ld", rank, i, pctl_mid[i], i, lcount[i]);
	}
/* Step 2: calculate the displacement and counts for sending and receiving */
	int senddispls[size];
	int sendcounts[size];

	senddispls[0] = 0;
	for (j = 1; j < size; j++) {
		senddispls[j] = searchsorted(pctl_mid + j - 1, psys.pool, psys.used, sizeof(Particle), SEARCH_DIR_LEFT, morton_key_t_compare, particle_hash_key, sizeof(morton_key_t));
		sendcounts[j - 1] = senddispls[j] - senddispls[j - 1];
	}
	sendcounts[j - 1] = psys.used - senddispls[j - 1];

/* Step 3: exchange particles */
	int recvdispls[size];
	int recvcounts[size];

	MPI_Alltoall(sendcounts, 1, MPI_INT, recvcounts, 1, MPI_INT, MPI_COMM_WORLD);
	recvdispls[0] = 0;
	for(j = 1; j < size; j++) {
		recvdispls[j] = recvdispls[j - 1] + recvcounts[j - 1];
	}
	int sum_recvcount = recvdispls[size - 1] + recvcounts[size - 1];

	Particle * pool2 = malloc(sizeof(Particle) * sum_recvcount * 1.2);
	
	MPI_Alltoallv(psys.pool, sendcounts, senddispls, MPI_PARTICLE,
			pool2, recvcounts, recvdispls, MPI_PARTICLE, MPI_COMM_WORLD);

	free(psys.pool);
	psys.pool = pool2;
	psys.size = sum_recvcount * 1.2;
	psys.used = sum_recvcount;

	intptr_t *arg = malloc(sizeof(intptr_t) * psys.used);
	argsort(psys.pool, psys.used, sizeof(Particle), arg, morton_key_t_compare, particle_hash_key, sizeof(morton_key_t));
	argpermute(psys.pool, NULL, psys.used, sizeof(Particle), arg);

	free(arg);

	size_t total_exchange = 0;

	for(j = 0; j < size; j++) {
		if(j != rank) total_exchange += recvcounts[j];
	}
	size_t all_total_exchange = 0;
	MPI_Allreduce(&total_exchange, &all_total_exchange, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

	size_t gnpar_new;
	MPI_Allreduce(&psys.used, &gnpar_new, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	if(rank == 0) {
		printf("total exchange = %ld, all_psys_used = %ld(new) %ld(old)\n", all_total_exchange, gnpar_new, gnpar);
	}

	return 0;
}

int in_square(uint64_t sqkey, int order, uint64_t k2) {
	return 0 == ((sqkey ^ k2) >> (order * 3));
}
int disjoint_squares(uint64_t k1, int o1, uint64_t k2, int o2) {
	int o = o1;
	if(o2 > o1) o = o2;
	return 0 != ((k1 ^ k2) >> (o * 3));
}
uint64_t derived_square(uint64_t k1, uint64_t k2, int * order) {
	/* returns the label of the square and the size in log2 in *size. or squaresize = 1.0 / (1<<size)*/
	/* when *order is nonzero, it is the size of the derived square labeld by k1 */
	/* the returning *order is always >= input *order, and the lower bits in return value are all cleared  */
	uint64_t diff = k1 ^ k2;
	uint64_t diff_save = diff;
	int r = 0;
	while(diff != 0) {
		diff >>= 3;
		r++;
	}
	if( r < *order) r = *order;
	else *order = r;
	uint64_t mask = 0;
	while(r > 0) {
		mask = (mask << 3) + 7;
		r--;
	}
	mask = ~mask;
	if((k1 & mask) != (k2 & mask)) abort();
	return k1 & mask;
}

intptr_t create_child(intptr_t first_par, intptr_t parent) {
	/* creates a child of parent from first_par, returns the new child */
	fwork.pool[fwork.used].first_par = first_par;
	fwork.pool[fwork.used].npar = 1;
	fwork.pool[fwork.used].parent = parent;
	fwork.pool[fwork.used].incomplete = 0;
	fwork.pool[fwork.used].order = fwork.pool[parent].order - 1;
	/* the lower bits of a sqkey is cleared off but I don't think it is necessary */
	fwork.pool[fwork.used].key = (morton_key(psys.pool[first_par].pos) >> (fwork.pool[fwork.used].order * 3)) << (fwork.pool[fwork.used].order * 3);

	fwork.pool[parent].child[fwork.pool[parent].child_length] = fwork.used;
	fwork.pool[parent].child_length ++;
	if(fwork.pool[parent].child_length == 9) abort();
	intptr_t rt = fwork.used;
	fwork.used++;
	if(fwork.used == fwork.size) {
		fwork.size *= 1.2;
		fwork.pool = realloc(fwork.pool, sizeof(FNode) * fwork.size);
	}
	return rt;
}

void split_child(intptr_t child) {
	intptr_t i = fwork.pool[child].first_par;
	size_t last = fwork.pool[child].npar + i;
	intptr_t j = create_child(i, child);
	/* i is already in the child */
	i++;
	while(i < last) {
		while(in_square(fwork.pool[j].key, fwork.pool[j].order, morton_key(PSYS(pos, i)))) {
			fwork.pool[j].npar ++;
			i++;
		}
		if(i < last) {
			j = create_child(i, child);
			i++;
		}
	}
}

int build_framework() {
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	intptr_t i;
	intptr_t j;
	fwork.size = 1024;
	fwork.pool = malloc(sizeof(FNode) * fwork.size);
	fwork.pool[0].key = 0;
	fwork.pool[0].first_par = 0;
	fwork.pool[0].npar = 0;
	fwork.pool[0].order = BITS;
	fwork.pool[0].incomplete = 0;
	fwork.pool[0].parent = -1;
	fwork.pool[0].child_length = 0;
	j = 0;
	fwork.used = 1;
	for(i = 0; i < psys.used; i++) {
		while(!in_square(fwork.pool[j].key, fwork.pool[j].order, morton_key(PSYS(pos, i)))) {
			j = fwork.pool[j].parent;
			if(j == -1) abort();
		}
		/* because we are on a morton key ordered list, no need to deccent into children */
		if(fwork.pool[j].child_length > 0 || (fwork.pool[j].npar >= THRESH && fwork.pool[j].order > 0)) {
			if(fwork.pool[j].child_length == 0) {
				/* splitting a leaf */
				i = fwork.pool[j].first_par;
				/* otherwise appending to a new child */
			}
			j = create_child(i, j);
		} else {
			/* put the particle into the leaf. */
			fwork.pool[j].npar ++;
		}
	}

	/* Step 2: refine the head and tail nodes to remove leaf overlappings */
    /* 
     * 1 Compare the last node of task i and the first node of task i+1. 
     * 2 If they are disjoint, good to go.
     * 3 Otherwise: If task are of different orders,
     *     refine the higher order node, goto 1.
     * 4 Now we have two identical nodes on both processors and both are nonempty. 
     *   Move particles from i + 1 to i, remove the node on rank i + 1 from the tree.
     *    (node still in memory but there is no way to access it from the tree.
     *     also remove all particles from the node. in case it gets accessed the particles
     *     won't be accessed twice. particles also still reside in the memory because we do 
     *     not want to update all indices  )
     *  */
    /* how to refine */
	int ggood = 0;
	while(ggood < size * 2) {
		intptr_t first_leaf;
		intptr_t last_leaf;
		first_leaf = 0;
		while(fwork.pool[first_leaf].child_length > 0) {
			first_leaf = fwork.pool[first_leaf].child[0];
		}
		last_leaf = 0;
		while(fwork.pool[last_leaf].child_length > 0) {
			last_leaf = fwork.pool[last_leaf].child[fwork.pool[last_leaf].child_length-1];
		}
		printf("%d %04ld %04ld\n", rank, first_leaf, last_leaf);
		uint64_t next_leaf_key;
		uint64_t prev_leaf_key;
		int next_leaf_order;
		int prev_leaf_order;
		MPI_Sendrecv(&fwork.pool[first_leaf].key, 1, MPI_UNSIGNED_LONG, (rank - 1 + size) % size, 0,
			&next_leaf_key, 1, MPI_UNSIGNED_LONG, (rank + 1) % size, 0, 
			MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Sendrecv(&fwork.pool[first_leaf].order, 1, MPI_INT, (rank - 1 + size) % size, 1,
			&next_leaf_order, 1, MPI_INT, (rank + 1) % size, 1, 
			MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Sendrecv(&fwork.pool[last_leaf].key, 1, MPI_UNSIGNED_LONG, (rank + 1 + size) % size, 2,
			&prev_leaf_key, 1, MPI_UNSIGNED_LONG, (rank - 1) % size, 2, 
			MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Sendrecv(&fwork.pool[last_leaf].order, 1, MPI_INT, (rank + 1 + size) % size, 3,
			&prev_leaf_order, 1, MPI_INT, (rank - 1) % size, 3, 
			MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		int next_good = 0;
		int prev_good = 0;
		if(disjoint_squares(fwork.pool[last_leaf].key, fwork.pool[last_leaf].order,
						next_leaf_key, next_leaf_order)) {
			next_good = 1;
		} else {
			if(next_leaf_order <= fwork.pool[last_leaf].order) {
				split_child(last_leaf);
			} /* the other case is handled at the first, prev of next process */
		}
		if(disjoint_squares(fwork.pool[first_leaf].key, fwork.pool[first_leaf].order,
						prev_leaf_key, prev_leaf_order)
		) {
			prev_good = 1;
		} else {
			if(!prev_good && prev_leaf_order <= fwork.pool[first_leaf].order) {
				split_child(first_leaf);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		int lgood = next_good + prev_good;
		MPI_Allreduce(&lgood, &ggood, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	}
	/* after this point some of the nodes will have order of 0. Be aware what is the size of
     * a order zero node? */
}

#if 0
		MPI_Request request = {0};
		size_t nrecv = 0;
		if(next_leaf_order == fwork.pool[last_leaf].order) {
			MPI_Irecv(&nrecv, 1, MPI_UNSIGNED_LONG, (rank + 1) % size, 0,
				MPI_COMM_WORLD, &request);
		}
		if(prev_leaf_order == fwork.pool[first_leaf].order) {
			MPI_Send(&fwork.pool[first_leaf].npar, 1, MPI_UNSIGNED_LONG, (rank - 1 + size) % size, 0,
				MPI_COMM_WORLD);
		}
		if(next_leaf_order == fwork.pool[last_leaf].order) {
			MPI_Wait(&request, MPI_STATUS_IGNORE);
			MPI_Irecv(&psys.pool[psys.used], nrecv, MPI_PARTICLE, (rank + 1) % size, 1, 
				MPI_COMM_WORLD, &request);
		}
		if(prev_leaf_order == fwork.pool[first_leaf].order) {
			MPI_Send(&psys.pool[fwork.pool[first_leaf].first_par],
				fwork.pool[first_leaf].npar, MPI_PARTICLE, (rank - 1 + size) % size, 1,
				MPI_COMM_WORLD);
			}
		if(next_leaf_order == fwork.pool[last_leaf].order) {
			MPI_Wait(&request, MPI_STATUS_IGNORE);
		}
		if(next_leaf_order == fwork.pool[last_leaf].order) {
			psys.used += nrecv;
		}
		if(prev_leaf_order == fwork.pool[first_leaf].order) {
			intptr_t j = first_leaf;
			while(fwork.pool[j].child_length == 0) {
				j = fwork.pool[j].parent;
				/* always the first child can be empty */
				for(i = 1; i < fwork.pool[j].child_length; i++) {
					fwork.pool[j].child[i - 1] = fwork.pool[j].child[i];
				}
				fwork.pool[j].child_length --;
			}
		}
#endif
size_t nearest_smaller_value(void* array, intptr_t output[], size_t nmemb, size_t bsize, int (*cmp)(const void* p1, const void * p2), int dir) {
/*
S = new empty stack data structure
for x in the input sequence:
    while S is nonempty and the top element of S is greater than or equal to x:
        pop S
    if S is empty:
        x has no preceding smaller value
    else:
        the nearest smaller value to x is the top element of S
    push x onto S

*/

	intptr_t i, j;
	intptr_t top = -1;
	size_t stacksize = nmemb / 1024;
	if(stacksize < 8) stacksize = 8;
	intptr_t * stack = malloc(sizeof(intptr_t) * stacksize);
	if(dir == 0) abort();
	for(j = 0; j < nmemb; j++) {
		if(dir > 0) i = nmemb - j - 1;
		if(dir < 0) i = j;
		while(top >= 0 && cmp(array + stack[top] * bsize, array + i * bsize) >=0) {
			top--;
		}
		if(top < 0) {
			output[i] = -1;
		} else {
			output[i] = stack[top];
		}
		stack[++top] = i;
		if(top == stacksize - 1) {
			stacksize *= 2;
			if(stacksize > nmemb) stacksize = nmemb;
			stack = realloc(stack, stacksize * sizeof(intptr_t));
		}
	}
	free(stack);
	return stacksize;
}
