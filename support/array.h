#ifndef __ARRAY_H__
#include <stdlib.h>
#include <string.h>
#define __ARRAY_H__

#define ARRAY_DEFINE(array, type) \
type * array = NULL; \
size_t array ## _length = 0; \
size_t array ## _size = 0;

#define ARRAY_DEFINE_S(array, type) \
type * array; \
size_t array ## _length; \
size_t array ## _size;

#define ARRAY_CLEAR(array) (array ## _length = 0, array)

#define ARRAY_ENSURE(array, type, size) \
	((array == NULL) && \
		(array ## _size = 128, \
		array = calloc(sizeof(type), array ## _size)), \
	(array ## _size < (size)) &&  \
		(array ## _size = __array_roundup__(array ## _size, (size)), \
		array = realloc(array, sizeof(type) * array ## _size)), \
	array)

#define ARRAY_ENSURE0(array, type, size) \
	((array == NULL) && \
		(array ## _size = 128, \
		array = calloc(sizeof(type), array ## _size)), \
	(array ## _size < (size)) &&  \
		(array = __array_realloc0__(array, array ## _size, (size), sizeof(type)), array ## _size = __array_roundup__(array ## _size, (size))), \
	array)

#define ARRAY_RESIZE(array, type, length) \
	(ARRAY_ENSURE0(array, type, (length)), array ## _length = (length), array)

#define ARRAY_APPEND(array, type) \
	(array ## _length = array ## _length + 1, \
	&(ARRAY_ENSURE(array, type, array ## _length)[array ## _length - 1]))

/*append element, if there is free element in the chain return it untainted
 * otherwise return a zeroed new element */
#define ARRAY_APPEND_REUSE(array, type) \
	(array ## _length = array ## _length + 1, \
	&(ARRAY_ENSURE0(array, type, array ## _length)[array ## _length - 1]))

#define ARRAY_FREE(array) \
	(free(array), array ## _length = 0, array ## _size = 0, array = NULL, array)
static inline size_t __array_roundup__(size_t old, size_t new_) {
	while(old < new_) old <<= 1;
	return old;
}
static inline void * __array_realloc0__(void * pt, size_t old, size_t new_, size_t unit) {
	size_t new__ = __array_roundup__(old, new_);
	pt = realloc(pt, new__ * unit);
	memset((char*)pt + old * unit, 0, (new__ - old) * unit);
	return pt;
}


#define SEARCH_DIR_LEFT -1 
#define SEARCH_DIR_RIGHT 1

static const void * _hash(void (*hash)(const void * p, void * k), const void * p, void * k) {
	if(hash) { hash(p, k); return k; }
	return p;
}

static intptr_t searchsorted(const void * key, const void * array, 
		const size_t nmemb, const size_t bsize, const int flag,
		int (*keycmp)(const void * k1, const void * k2), 
		void (*hash)(const void * p, void * k), 
		const size_t keysize) {

	intptr_t mid = nmemb / 2;
	const size_t ksize = hash?keysize:bsize;
	int c = 0;
	char keybuf[ksize];
	const void *otherkey;

	otherkey = _hash(hash, array, keybuf);

	c = keycmp(key, otherkey) ;
	if(c < 0 || (flag == SEARCH_DIR_LEFT && c == 0)) return 0;
	otherkey = _hash(hash, array + (nmemb - 1) * bsize, keybuf);
	c = keycmp(key, otherkey) ;
	if(c > 0 || (flag == SEARCH_DIR_RIGHT && c == 0)) return nmemb;
	otherkey = _hash(hash, array + mid * bsize, keybuf);
	c = keycmp(key, otherkey);
	if(c > 0 || ((flag == SEARCH_DIR_RIGHT) && c == 0))
		return searchsorted(key, array + mid * bsize, nmemb - mid, bsize, flag, keycmp, hash, ksize) + mid;
	if(c < 0 || (flag == SEARCH_DIR_LEFT && c == 0))
		return searchsorted(key, array, mid, bsize, flag, keycmp, hash, ksize);
	abort();
}

static void argmaxmin(const void * array, 
		const size_t nmemb, const size_t bsize, 
		intptr_t * argmax, intptr_t * argmin, 
		int (*keycmp)(const void * k1, const void * k2), 
		void (*hash)(const void *p, void * k), 
		const size_t keysize) {
	intptr_t i;
	const size_t ksize = hash?keysize:bsize;
	char min[ksize], max[ksize];
	char keybuf[ksize];
	const void * key;
	if(nmemb == 0) {
		if(argmin) *argmin = -1;
		if(argmax) *argmax = -1;
		return;
	}
	/* else */
	if(argmin) {
		*argmin = 0;
		key = _hash(hash, array, keybuf);
		memcpy(min, key, ksize);
	}
	if(argmax) {
		*argmax = 0;
		key = _hash(hash, array, keybuf);
		memcpy(max, key, ksize);
	}
	for(i = 1; i < nmemb; i ++ ) {
		key = _hash(hash, array + i * bsize, keybuf);
		if(argmin && keycmp(key, min) < 0) {
			*argmin = i;
			memcpy(min, key, ksize);
		}
		if(argmax && keycmp(key, max) > 0) {
			*argmax = i;
			memcpy(max, key, ksize);
		}
	}
}

static void argsort(const void * array,
		const size_t nmemb, const size_t bsize,
		intptr_t * arg,
		int (*keycmp)(const void * p1, const void * p2),
		void (*hash)(const void * p, void * k),
		const size_t keysize) {
	intptr_t * argb = malloc(sizeof(intptr_t) * nmemb);
	intptr_t * A = arg;
	intptr_t * B = argb;
	intptr_t * tmp;

	const size_t ksize = hash?keysize:bsize;
	char k0buf[ksize], k1buf[ksize];
	intptr_t width;
	for(width = 1; width < nmemb; width = 2 * width) {
		/* decide which arg buf to use so that when exiting
         * the correct one, arg contains the result A. We need
         * an even total number of swaps. So we preswap the
         * same numbers first */
		tmp = A;
		A = B;
		B = tmp;
	}
	/* and fill in A with the original order */
	intptr_t i;
	for(i = 0; i < nmemb; i++) {
		A[i] = i;
	}
	for(width = 1; width < nmemb; width = 2 * width) {
		intptr_t i;
		for(i = 0; i < nmemb; i += 2 * width) {
			intptr_t pivot = (i+width<nmemb)?(i+width):nmemb;
			intptr_t end = (i+2*width<nmemb)?(i+2*width):nmemb;
			intptr_t i0 = i;
			intptr_t i1 = pivot;
			intptr_t j;
			for(j = i; j < end; j++) {
				const void * k0, * k1;
				if(i0 < pivot && (i1 >= end || (
					k0 = _hash(hash, array + bsize * A[i0], k0buf),
					k1 = _hash(hash, array + bsize * A[i1], k1buf),
					keycmp(k0, k1) <= 0)
				)) {
					B[j] = A[i0++];
				} else {
					B[j] = A[i1++];
				}
			}
		}
		/* swap A and B */
		tmp = A;
		A = B;
		B = tmp;
	}
	free(argb);
}
static void maxmin(const void * array, 
		const size_t nmemb, const size_t bsize, 
		void * max, void * min, 
		int (*keycmp)(const void * p1, const void * p2), 
		void (*hash)(const void * p, void * k), 
		const size_t keysize) {
	intptr_t argmax, argmin;
	if(nmemb == 0) return;
	argmaxmin(array, nmemb, bsize, &argmax, &argmin, keycmp, hash, keysize);

	const size_t ksize = hash?keysize:bsize;
	char keybuf[ksize];
	const void * key;
	if(max) {
		key = _hash(hash, array + argmax * bsize, keybuf);
		memcpy(max, key, ksize);
	}
	if(min) {
		key = _hash(hash, array + argmin * bsize, keybuf);
		memcpy(min, key, ksize);
	}
}

static void argpermute (void * array, void * target, size_t nmemb, size_t bsize, intptr_t * p) {
	intptr_t i, k, pk;

	if(target != NULL && target != array) {
		for(i = 0; i < nmemb; i++) {
			memcpy(target + i * bsize, array + p[i] * bsize, bsize);
		}
	}
	for (i = 0; i < nmemb; i++)
	{
		k = p[i];

		while (k > i) 
			k = p[k];

		if (k < i)
			continue ;

		/* Now have k == i, i.e the least in its cycle */

		pk = p[k];

		if (pk == i)
			continue ;

		/* shuffle the elements of the cycle */

		{
			unsigned int a;

			char t[bsize];

			memcpy(t, array + i * bsize, bsize);

			while (pk != i)
			{
				memcpy(array + k * bsize, array + pk * bsize, bsize);
				k = pk;
				pk = p[k];
			};

			memcpy(array + k * bsize, t, bsize);
		}
	}

}

#endif
