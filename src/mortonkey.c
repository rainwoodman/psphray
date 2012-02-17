#include <stdio.h>

#include <mpi.h>
#include "mortonkey.h"

MPI_Datatype MPI_TYPE_MORTON_KEY;
MPI_Op MPI_OP_MORTON_KEY_MIN;
MPI_Op MPI_OP_MORTON_KEY_MAX;

int morton_key_bits = 0;

static void _morton_key_mpi_min(void * invec, void * inoutvec, int* len, MPI_Datatype *type) {
	int i;
	for(i = 0; i < *len; i++) {
		if(morton_key_t_compare(((morton_key_t*)invec) + i, ((morton_key_t*)inoutvec) + i) < 0) {
			((morton_key_t*)inoutvec)[i] = ((morton_key_t*)invec)[i];
		}
	}
}
static void _morton_key_mpi_max(void * invec, void * inoutvec, int* len, MPI_Datatype *type) {
	int i;
	for(i = 0; i < *len; i++) {
		if(morton_key_t_compare(((morton_key_t*)invec) + i, ((morton_key_t*)inoutvec) + i) > 0) {
			((morton_key_t*)inoutvec)[i] = ((morton_key_t*)invec)[i];
		}
	}
}
void morton_key_init(int bits) {
	morton_key_bits = bits;
#ifdef MPI_INCLUDED
	MPI_Type_contiguous(sizeof(morton_key_t), MPI_BYTE, &MPI_TYPE_MORTON_KEY);
	MPI_Type_commit(&MPI_TYPE_MORTON_KEY);
	MPI_Op_create(_morton_key_mpi_min, 1, &MPI_OP_MORTON_KEY_MIN);
	MPI_Op_create(_morton_key_mpi_max, 1, &MPI_OP_MORTON_KEY_MAX);
#endif
}

