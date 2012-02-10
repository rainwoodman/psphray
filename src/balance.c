#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <printf.h>
#include <mpi.h>
#define N 32
#define BITS 12
#define THRESH 4
typedef struct {
	float pos[3];
	uint64_t key;
} Particle;

typedef struct {
	Particle * pool;
	size_t size;
	size_t used;
} PSystem;

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

const uint64_t morton_key(const float pos[3]) {
	int mask;
	uint64_t morton = 0;
	int ix,iy,iz;
	const int bits = BITS;
	ix = pos[0] * (1 << bits);
	iy = pos[1] * (1 << bits);
	iz = pos[2] * (1 << bits);
	for(mask = 1 << (bits - 1); mask > 0; mask >>= 1)
	{
		morton <<= 3;
		morton += ((iz & mask) ? 4 : 0) + ((iy & mask) ? 2 : 0) + ((ix & mask) ? 1 : 0);
	}

	return morton;
}


PSystem psys;
Framework fwork;
MPI_Datatype MPI_PARTICLE;

int printf_binary(FILE * stream, const struct printf_info * info, const void * const *args);
int printf_binary_arginfo (const struct printf_info *info, size_t n, int *argtypes);
size_t nearest_smaller_value(void* array, intptr_t output[], size_t nmemb, size_t bsize, int (*cmp)(const void* p1, const void * p2), int dir);

int main(int argc, char* argv[]) {

	register_printf_function('b', printf_binary, printf_binary_arginfo);
	MPI_Init(&argc, &argv);
	MPI_Type_contiguous(sizeof(Particle), MPI_BYTE, &MPI_PARTICLE);
	MPI_Type_commit(&MPI_PARTICLE);
	mockdata();
	balance();
	build_framework();

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	printf("rank = %d, psys.used/psys.size= %ld/%ld\n", rank, psys.used, psys.size);
	char * fname = NULL;
	asprintf(&fname, "balance-out-%03d", rank);
	FILE * fp = fopen(fname, "w");
	intptr_t i;
	for(i = 0; i < psys.used; i++) {
		fprintf(fp, "%08lo %g %g %g\n", psys.pool[i].key,
			psys.pool[i].pos[0],
			psys.pool[i].pos[1],
			psys.pool[i].pos[2]);
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
				fprintf(fp, "%08lo ", psys.pool[j + fwork.pool[i].first_par].key);
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


static int uint64_compare(const void * p1, const void * p2) {
	return (*(uint64_t*)p1 > *(uint64_t*)p2) - (*(uint64_t*)p1 < *(uint64_t*)p2);
}
static int key_compare(const void* p1, const void* p2) {
	return uint64_compare(&((Particle*)p1)->key, &((Particle*)p2)->key);
}
int mockdata() {
	psys.size = N *1.2;
	psys.used = N * (0.5 * random() / RAND_MAX + 0.3);
	psys.pool = malloc(sizeof(Particle) * psys.size);
	intptr_t i;
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	for(i = 0; i < N; i++) {
		int d;
		for(d = 0; d < 3; d++) {
			psys.pool[i].pos[d] = fmin((float) random() / RAND_MAX, 1.0);
		}
		psys.pool[i].key = morton_key(psys.pool[i].pos);
	}
	qsort(psys.pool, psys.used, sizeof(Particle), key_compare);
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

#define DIR_LEFT -1 
#define DIR_RIGHT 1
intptr_t searchsorted(void * data, void * array, size_t nmemb, size_t bsize, int (*cmp)(const void * p1, const void * p2), int flag) {
	/**/
	intptr_t mid = nmemb / 2;
	int c = 0;
	c = cmp(data, ((char*)array)) ;
	if(c < 0 || (flag == DIR_LEFT && c == 0)) return 0;
	c = cmp(data, ((char*)array + (nmemb - 1) * bsize)) ;
	if(c > 0 || (flag == DIR_RIGHT && c == 0)) return nmemb;
	c = cmp(data, ((char*)array + mid * bsize));
	if(c > 0 || ((flag == DIR_RIGHT) && c == 0))
		return searchsorted(data, (char*) array + mid * bsize, nmemb - mid, bsize, cmp, flag) + mid;
	if(c < 0 || (flag == DIR_LEFT && c == 0))
		return searchsorted(data, (char*) array, mid, bsize, cmp, flag);
	abort();
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
	uint64_t pctl_min[size];
	uint64_t pctl_max[size];
	intptr_t lmax = psys.pool[0].key;
	intptr_t lmin = psys.pool[0].key;
	for(i = 0; i < psys.used; i++) {
		if(psys.pool[i].key > lmax) lmax = psys.pool[i].key;
		if(psys.pool[i].key < lmin) lmin = psys.pool[i].key;
	}
	intptr_t gmax = 0;
	intptr_t gmin = 0;
	MPI_Allreduce(&lmax, &gmax, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(&lmin, &gmin, 1, MPI_UNSIGNED_LONG, MPI_MIN, MPI_COMM_WORLD);

	for(j = 0; j < size; j++) {
		pctl_max[j] = gmax;
		pctl_min[j] = gmin;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	intptr_t lcount[size];
	uint64_t pctl_mid[size];

	while(1) {
		for(j = 0; j < size; j++) {
			pctl_mid[j] = (pctl_max[j] + pctl_min[j]) / 2;
		}
		for(j = 0; j < size; j++) {
			Particle p;
			p.key = pctl_mid[j];
			if(j > 0) {
				/* omit some countings on the first few rounds because all pctls start the same */
				if(pctl_mid[j] != pctl_mid[j -1]) {
					lcount[j] = searchsorted(&p, psys.pool, psys.used, sizeof(Particle), key_compare, DIR_LEFT);
				} else {
					lcount[j] = lcount[j - 1];
				}
			} else {
				lcount[j] = searchsorted(&p, psys.pool, psys.used, sizeof(Particle), key_compare, DIR_LEFT);
			}
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

/* Step 2: calculate the displacement and counts for sending and receiving */
	int senddispls[size];
	int sendcounts[size];

	senddispls[0] = 0;
	for (j = 1; j < size; j++) {
		Particle p;
		p.key = pctl_mid[j - 1];
		senddispls[j] = searchsorted(&p, psys.pool, psys.used, sizeof(Particle), key_compare, DIR_LEFT);
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

	qsort(psys.pool, psys.used, sizeof(Particle), key_compare);
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
	fwork.pool[fwork.used].key = (psys.pool[first_par].key >> (fwork.pool[fwork.used].order * 3)) << (fwork.pool[fwork.used].order * 3);

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
		while(in_square(fwork.pool[j].key, fwork.pool[j].order, psys.pool[i].key)) {
			fwork.pool[j].npar ++;
			i++;
		}
		j = create_child(i, child);
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
		while(!in_square(fwork.pool[j].key, fwork.pool[j].order, psys.pool[i].key)) {
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

	/* Step 2: fillin the incomplete flags */
	/* The incomplete nodes are only in the end or at the front of the leaf nodes */
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
	while(ggood < size) {
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
		}
		if(disjoint_squares(fwork.pool[first_leaf].key, fwork.pool[first_leaf].order,
						prev_leaf_key, prev_leaf_order)
		) {
			prev_good = 1;
		}
		if(prev_leaf_order < fwork.pool[first_leaf].order) {
			/* ::=> last_leaf_order < next_leaf_order on previous rank */
			split_child(first_leaf);
		}
		if(next_leaf_order < fwork.pool[last_leaf].order) {
			/* ::=> last_leaf_order < next_leaf_order on previous rank */
			split_child(last_leaf);
		}
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
		MPI_Barrier(MPI_COMM_WORLD);
	}
}

int printf_binary(FILE * stream, const struct printf_info * info, const void * const *args) {
	const int bits = info->width;
	int buflen = bits + 1;
	if(info->prec != -1) {
		buflen == bits / info->prec + 1;
	}
	char buf[buflen];
	int p = 0;
	unsigned int mask = (1 << (bits - 1));
	intptr_t l;
	if(info->is_long) {
		l = *(long * )(args[0]);
	} else if(info->is_short) {
		l = *(short *)(args[0]);
	} else if(info->is_char) {
		l = *(char *)(args[0]);
	} else {
		l = *(int *)(args[0]);
	}
	int j = 0;
	while(p < bits) {
		buf[j] = (l & mask)?'1':'0';
		mask >>= 1;
		p++;
		j++;
		if((info->prec != -1) && ((bits - p) % info->prec == 0)) {
			buf[j] = ' ';
			j++;
		}
	}
	buf[j] = 0;
	fputs(buf, stream);
	return bits;
}
int printf_binary_arginfo (const struct printf_info *info, size_t n, int *argtypes) {
	if (n > 0)
		argtypes[0] = PA_POINTER;
	return 1;
}

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
