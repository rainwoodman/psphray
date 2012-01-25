#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#define max(a, b) ((a) > (b)?(a):(b))
#define min(a, b) ((a) < (b)?(a):(b))

struct rtree_t {
	struct rtree_node_t * pool;
	intptr_t used;
	size_t size;
};

struct rtree_node_t {
/* if type == leaf, first_child = first par, if type == nonleaf first_child = first child node*/
	int type;  /*0 = leaf, 1= nonleaf, nonlastsibling, 2=lastsibling */
/* number of children or pars */
	int length;
/* bounding box */
	float bot[3];
	float top[3];
	intptr_t parent;
	intptr_t first_child;
	intptr_t lhv; /* largest hilbert key stored in the node */
};


void rtree_build(struct rtree_t * rtree, float pos[][3], float sml[], intptr_t key[], size_t npar) {
	/* build a rtree, in input, pos is sorted by the key */

	size_t leaves = (npar + 11)/ 12;
	size_t nonleaves = 0;
	size_t t;
	for(t = (leaves + 7)/ 8; t >= 1; t = (t + 7)/ 8) {
		nonleaves += t;
		if(t == 1) break;
	}

	rtree->size = leaves + nonleaves;
	rtree->pool = malloc(sizeof(struct rtree_node_t) * rtree->size);
	rtree->used = 0;
	intptr_t i;
	/* fill the leaves*/
	for(i = 0; i < leaves; i++) {
		struct rtree_node_t * node = &rtree->pool[i];
		node->first_child = i * 12;
		node->type = 0;
		node->length = min(npar - i * 12, 12);
		node->lhv = key[node->first_child + node->length - 1];
		int d;
		for(d = 0; d < 3; d++) {
			node->bot[d] = pos[node->first_child][d];
			node->top[d] = pos[node->first_child][d];
		}
		intptr_t j;
		for(j = 0; j < node->length; j++) {
			for(d = 0; d < 3; d++) {
				node->bot[d] = fmin(node->bot[d],
						pos[node->first_child + j][d] - sml[node->first_child + j]);
				node->top[d] = fmax(node->top[d], 
						pos[node->first_child + j][d] + sml[node->first_child + j]);
			}
		}
		rtree->used ++;
	}
	intptr_t child_cur = 0;
	intptr_t parent_cur = leaves;
	for(t = (leaves + 7)/ 8; t >= 1; t = (t + 7)/ 8) {
		for(i = 0; i < t; i++) {
			struct rtree_node_t * node = &rtree->pool[parent_cur + i];
			node->type = 1;
			node->first_child = child_cur + i * 8;
			node->length = min(parent_cur - child_cur - i * 8, 8);
			intptr_t j;
			node->lhv = 0;
			int d;
			struct rtree_node_t * child = &rtree->pool[node->first_child];
			for(d = 0; d < 3; d++) {
				node->bot[d] = child[0].bot[d];
				node->top[d] = child[0].top[d];
			}
			for(j = 0; j < node->length; j++) {
				child[j].parent = parent_cur + i;
				node->lhv = max(node->lhv, child[j].lhv);
				for(d = 0; d < 3; d++) {
					node->bot[d] = fmin(node->bot[d], child[j].bot[d]);
					node->top[d] = fmax(node->top[d], child[j].top[d]);
				}
			}
			child[j-1].type += 2;
		rtree->used ++;
		}
		child_cur = parent_cur;
		parent_cur += t;
		if(t == 1) break;
	}

}

int intersect(float top1[], float bot1[], float top2[], float bot2[]) {
	int d;
/*
	printf("%g %g %g - %g %g %g x %g %g %g - %g %g %g",
		top1[0], top1[1], top1[2],
		bot1[0], bot1[1], bot1[2],
		top2[0], top2[1], top2[2],
		bot2[0], bot2[1], bot2[2]);
*/
	for(d = 0 ;d < 3; d++) {
		if(bot1[d] > top2[d]) {
/*
			printf("false\n");
*/
			return 0;
		}
		if(bot2[d] > top1[d]) {
/*
			printf("false\n");
*/
			return 0;
		}
	}
//	printf("true\n");
	return 1;
}
void add_to_list(intptr_t ipar, float pos[][3], float p[3], intptr_t neighbours[], float dist[], int n, int * nused) {
	int d;
	float d2 = 0.0;
	for(d = 0; d < 3; d++) {
		float s = pos[ipar][d] - p[d];
		d2 += s * s;
	}
	float dd = sqrt(d2);
	if(*nused == n && dd > dist[n - 1]) return;
	int i;

	int ip = 0;
	while(ip < *nused && dist[ip] < dd) ip++;

	for(i = min(*nused, n - 1); i > ip; i--) {
		dist[i] = dist[i - 1];
		neighbours[i] = neighbours[i - 1];
	}
	dist[ip] = dd;
	neighbours[ip] = ipar;
	if(*nused < n) (*nused)++;
}
void rtree_neighbours(struct rtree_t * rtree, float p[3], float pos[][3], size_t npar,
		intptr_t neighbours[], float dist[], int n, int *nused, float * hhint) {
	float hh = 1.0;
	if(hhint == NULL) { hhint = &hh; }
	if(n > npar) n = npar;
	float bot[3];
	float top[3];
	int d;
	struct rtree_node_t * N = rtree->pool;
	do {
		for(d = 0; d < 3; d++) {
			bot[d] = p[d] - *hhint;
			top[d] = p[d] + *hhint;
		}
		intptr_t cur = rtree->used - 1;
		*nused = 0;
		while(1) {
			if(intersect(N[cur].top, N[cur].bot, top, bot)) {
				if(N[cur].type == 1 || N[cur].type == 3) {
					cur = N[cur].first_child;
					continue;
				}
				if(N[cur].type == 0 || N[cur].type == 2) {
					int j;
					for(j = 0; j < N[cur].length; j++) {
						add_to_list(N[cur].first_child + j, pos, p, neighbours, dist, n, nused);
					}
				}
			}
			if(cur == rtree->used - 1) break;
			if(N[cur].type == 0 || N[cur].type == 1) {
				cur++;
				continue;
			}
			while((N[cur].type == 2 || N[cur].type == 3)) {
				cur = N[cur].parent;
				if(cur == rtree->used - 1) break;
			}
			if(cur == rtree->used - 1) break;
			cur++;
		}
		(*hhint) *= 2.0;
	} while(*nused < n);
}

#ifdef __MAIN__

intptr_t peano_hilbert_key(int x, int y, int z, int bits);


#include <gsl/gsl_permutation.h>
#include <gsl/gsl_heapsort.h>
static void * permute (const size_t * p, void * data, const size_t ele_bytes, const size_t stride, const size_t n, const size_t n_max)
{
	void * out = malloc(stride * n_max);
	intptr_t i;
	#pragma omp parallel for private(i)
	for(i = 0; i < n; i++) {
		if(p[i] > n ) {
			printf("p[i] out of range check peano key\n");
		}
		memcpy((char*) out+ i * stride, (char*) data + p[i] * stride, ele_bytes);
	}
	return out;
}

static int intptr_t_compare(const intptr_t * p1, const intptr_t * p2) {
	if(*p1 > *p2) return 1;
	if(*p1 < *p2) return -1;
	if(*p1 == *p2) return 0;
}

void main() {
	size_t * perm = malloc(sizeof(size_t) * 1000);
	struct rtree_t t = {0};
	int i, j, k;
	float pos[1000][3];
	float sml[1000];
	intptr_t key[1000];
	int l = 0;
	for(i = 0; i < 10; i++) {
		for(j = 0; j < 10; j++) {
			for(k = 0; k < 10; k++) {
				pos[l][0] = i * 0.1 * (1<<20);
				pos[l][1] = j * 0.1 * (1<<20);
				pos[l][2] = k * 0.1 * (1<<20);
				sml[l] = 0;
				key[l] = peano_hilbert_key(pos[l][0], pos[l][1], pos[l][2], 20);
				l++;
			}
		}
	}

	gsl_heapsort_index(perm, key, 1000, sizeof(intptr_t), (void*)intptr_t_compare);

	float (*__pos)[3] = permute(perm, pos, 3 * sizeof(float), 3 * sizeof(float), 1000, 1000);
	float * __sml = permute(perm, sml, sizeof(float), sizeof(float), 1000, 1000);
	intptr_t * __key = permute(perm, key, sizeof(intptr_t), sizeof(intptr_t), 1000, 1000);
	rtree_build(&t, __pos, __sml, __key, 1000);
	float dist[10];
	intptr_t nei[10];
	int nused = 0;
	float hhint = 100000;
	float p[3] = {499999, 499999,499999};

	rtree_neighbours(&t, p, __pos, 1000, nei, dist, 10, &nused, NULL);
	printf("hello\n");
}
static const unsigned char rottable3[48][8] = {
  {36, 28, 25, 27, 10, 10, 25, 27},
  {29, 11, 24, 24, 37, 11, 26, 26},
  {8, 8, 25, 27, 30, 38, 25, 27},
  {9, 39, 24, 24, 9, 31, 26, 26},
  {40, 24, 44, 32, 40, 6, 44, 6},
  {25, 7, 33, 7, 41, 41, 45, 45},
  {4, 42, 4, 46, 26, 42, 34, 46},
  {43, 43, 47, 47, 5, 27, 5, 35},
  {33, 35, 36, 28, 33, 35, 2, 2},
  {32, 32, 29, 3, 34, 34, 37, 3},
  {33, 35, 0, 0, 33, 35, 30, 38},
  {32, 32, 1, 39, 34, 34, 1, 31},
  {24, 42, 32, 46, 14, 42, 14, 46},
  {43, 43, 47, 47, 25, 15, 33, 15},
  {40, 12, 44, 12, 40, 26, 44, 34},
  {13, 27, 13, 35, 41, 41, 45, 45},
  {28, 41, 28, 22, 38, 43, 38, 22},
  {42, 40, 23, 23, 29, 39, 29, 39},
  {41, 36, 20, 36, 43, 30, 20, 30},
  {37, 31, 37, 31, 42, 40, 21, 21},
  {28, 18, 28, 45, 38, 18, 38, 47},
  {19, 19, 46, 44, 29, 39, 29, 39},
  {16, 36, 45, 36, 16, 30, 47, 30},
  {37, 31, 37, 31, 17, 17, 46, 44},
  {12, 4, 1, 3, 34, 34, 1, 3},
  {5, 35, 0, 0, 13, 35, 2, 2},
  {32, 32, 1, 3, 6, 14, 1, 3},
  {33, 15, 0, 0, 33, 7, 2, 2},
  {16, 0, 20, 8, 16, 30, 20, 30},
  {1, 31, 9, 31, 17, 17, 21, 21},
  {28, 18, 28, 22, 2, 18, 10, 22},
  {19, 19, 23, 23, 29, 3, 29, 11},
  {9, 11, 12, 4, 9, 11, 26, 26},
  {8, 8, 5, 27, 10, 10, 13, 27},
  {9, 11, 24, 24, 9, 11, 6, 14},
  {8, 8, 25, 15, 10, 10, 25, 7},
  {0, 18, 8, 22, 38, 18, 38, 22},
  {19, 19, 23, 23, 1, 39, 9, 39},
  {16, 36, 20, 36, 16, 2, 20, 10},
  {37, 3, 37, 11, 17, 17, 21, 21},
  {4, 17, 4, 46, 14, 19, 14, 46},
  {18, 16, 47, 47, 5, 15, 5, 15},
  {17, 12, 44, 12, 19, 6, 44, 6},
  {13, 7, 13, 7, 18, 16, 45, 45},
  {4, 42, 4, 21, 14, 42, 14, 23},
  {43, 43, 22, 20, 5, 15, 5, 15},
  {40, 12, 21, 12, 40, 6, 23, 6},
  {13, 7, 13, 7, 41, 41, 22, 20}
};

static const unsigned char subpix3[48][8] = {
  {0, 7, 1, 6, 3, 4, 2, 5},
  {7, 4, 6, 5, 0, 3, 1, 2},
  {4, 3, 5, 2, 7, 0, 6, 1},
  {3, 0, 2, 1, 4, 7, 5, 6},
  {1, 0, 6, 7, 2, 3, 5, 4},
  {0, 3, 7, 4, 1, 2, 6, 5},
  {3, 2, 4, 5, 0, 1, 7, 6},
  {2, 1, 5, 6, 3, 0, 4, 7},
  {6, 1, 7, 0, 5, 2, 4, 3},
  {1, 2, 0, 3, 6, 5, 7, 4},
  {2, 5, 3, 4, 1, 6, 0, 7},
  {5, 6, 4, 7, 2, 1, 3, 0},
  {7, 6, 0, 1, 4, 5, 3, 2},
  {6, 5, 1, 2, 7, 4, 0, 3},
  {5, 4, 2, 3, 6, 7, 1, 0},
  {4, 7, 3, 0, 5, 6, 2, 1},
  {6, 7, 5, 4, 1, 0, 2, 3},
  {7, 0, 4, 3, 6, 1, 5, 2},
  {0, 1, 3, 2, 7, 6, 4, 5},
  {1, 6, 2, 5, 0, 7, 3, 4},
  {2, 3, 1, 0, 5, 4, 6, 7},
  {3, 4, 0, 7, 2, 5, 1, 6},
  {4, 5, 7, 6, 3, 2, 0, 1},
  {5, 2, 6, 1, 4, 3, 7, 0},
  {7, 0, 6, 1, 4, 3, 5, 2},
  {0, 3, 1, 2, 7, 4, 6, 5},
  {3, 4, 2, 5, 0, 7, 1, 6},
  {4, 7, 5, 6, 3, 0, 2, 1},
  {6, 7, 1, 0, 5, 4, 2, 3},
  {7, 4, 0, 3, 6, 5, 1, 2},
  {4, 5, 3, 2, 7, 6, 0, 1},
  {5, 6, 2, 1, 4, 7, 3, 0},
  {1, 6, 0, 7, 2, 5, 3, 4},
  {6, 5, 7, 4, 1, 2, 0, 3},
  {5, 2, 4, 3, 6, 1, 7, 0},
  {2, 1, 3, 0, 5, 6, 4, 7},
  {0, 1, 7, 6, 3, 2, 4, 5},
  {1, 2, 6, 5, 0, 3, 7, 4},
  {2, 3, 5, 4, 1, 0, 6, 7},
  {3, 0, 4, 7, 2, 1, 5, 6},
  {1, 0, 2, 3, 6, 7, 5, 4},
  {0, 7, 3, 4, 1, 6, 2, 5},
  {7, 6, 4, 5, 0, 1, 3, 2},
  {6, 1, 5, 2, 7, 0, 4, 3},
  {5, 4, 6, 7, 2, 3, 1, 0},
  {4, 3, 7, 0, 5, 2, 6, 1},
  {3, 2, 0, 1, 4, 5, 7, 6},
  {2, 5, 1, 6, 3, 4, 0, 7}
};

/*! This function computes a Peano-Hilbert key for an integer triplet (x,y,z),
  *  with x,y,z in the range between 0 and 2^bits-1.
  */
intptr_t peano_hilbert_key(int x, int y, int z, int bits)
{
  int mask;
  unsigned char rotation = 0;
  intptr_t key = 0;

  for(mask = 1 << (bits - 1); mask > 0; mask >>= 1)
    {
      unsigned char pix = ((x & mask) ? 4 : 0) | ((y & mask) ? 2 : 0) | ((z & mask) ? 1 : 0);

      key <<= 3;
      key |= subpix3[rotation][pix];
      rotation = rottable3[rotation][pix];
    }

  return key;
}

#endif

