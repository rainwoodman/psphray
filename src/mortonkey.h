#include <stdint.h>

typedef uint64_t morton_key_t;

#ifdef MPI_INCLUDED
extern MPI_Datatype MPI_TYPE_MORTON_KEY;
extern MPI_Op MPI_OP_MORTON_KEY_MIN;
extern MPI_Op MPI_OP_MORTON_KEY_MAX;
#endif

extern int morton_key_bits;
static int morton_key_t_compare(const void * p1, const void * p2) {
	return (*(morton_key_t*)p1 > *(morton_key_t*)p2) - (*(morton_key_t*)p1 < *(morton_key_t*)p2);
}

static const morton_key_t morton_key(const float pos[3]) {
	int mask;
	morton_key_t morton = 0;
	int ix,iy,iz;
	ix = pos[0] * (1 << morton_key_bits);
	iy = pos[1] * (1 << morton_key_bits);
	iz = pos[2] * (1 << morton_key_bits);
	for(mask = 1 << (morton_key_bits - 1); mask > 0; mask >>= 1)
	{
		morton <<= 3;
		morton += ((iz & mask) ? 4 : 0) + ((iy & mask) ? 2 : 0) + ((ix & mask) ? 1 : 0);
	}

	return morton;
}

