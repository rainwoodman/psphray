#define BLOCK_TYPE unsigned int

#define BLOCK (sizeof(BLOCK_TYPE) * 8)
static int __bitmask_log2[] = {
// 0 , 1, 2, 3,  4,  5,  6,  7, 8, 9, 10, ...
   -1, 3, 4, -1, 5, -1, -1, -1, 6, -1, -1, 
};

#define LOG2_BLOCK  5
//(__bitmask_log2[sizeof(BLOCK_TYPE)])
#define BLOCK_1 (BLOCK - 1)

typedef struct bitmask {
	size_t bits;
	size_t bits_alloc;
	size_t bytes;
	size_t length;
	BLOCK_TYPE data[];
} bitmask_t;

static inline bitmask_t * bitmask_alloc(size_t bits) {
	/* ensure the last 2 bytes so that we can shuffle around the bits for
       cacheline optimization of continueous access */
	size_t bits_alloc = bits;
	if(bits_alloc & 0xFFFF) bits_alloc = (bits_alloc + 0x10000) & ~ ((size_t) 0xFFFF);
	size_t length = (bits_alloc + BLOCK_1) >> LOG2_BLOCK;
	bits_alloc = length * BLOCK;
	bitmask_t * rt = malloc(length * sizeof(BLOCK_TYPE) + sizeof(bitmask_t));
	rt->bits = bits;
	rt->bits_alloc = bits_alloc;
	rt->length = length;
	rt->bytes = length * sizeof(BLOCK_TYPE);
	return rt;
}
static inline void bitmask_set_all(bitmask_t * mask) {
	memset(mask->data, ~0x00, mask->bytes);
}

static inline void bitmask_clear_all(bitmask_t * mask) {
	memset(mask->data, 0x00, mask->bytes);
}

static const inline intptr_t bitmask_shuffle(const intptr_t idx) {
#ifdef BITMASK_SHUFFLE
	uint16_t v = idx;
	v = ((v >> 1) & 0x5555) | ((v & 0x5555) << 1);
	v = ((v >> 2) & 0x3333) | ((v & 0x3333) << 2);
	v = ((v >> 4) & 0x0F0F) | ((v & 0x0F0F) << 4);
	v = ((v >> 8) & 0x00FF) | ((v & 0x00FF) << 8);
	
	return (idx & ~0xffff) | v;
#else
	return idx;
#endif
}

static inline void bitmask_set(bitmask_t * mask, intptr_t idx) {
	idx = bitmask_shuffle(idx);
	int offset = idx & BLOCK_1;
	BLOCK_TYPE bit = ((BLOCK_TYPE)1) << offset;
	__sync_or_and_fetch(&mask->data[idx >> LOG2_BLOCK], bit);
}
static inline void bitmask_set_unsafe(bitmask_t * mask, intptr_t idx) {
	idx = bitmask_shuffle(idx);
	int offset = idx & BLOCK_1;
	BLOCK_TYPE bit = ((BLOCK_TYPE)1) << offset;
	mask->data[idx >> LOG2_BLOCK] |= bit;

}
static inline void bitmask_clear(bitmask_t * mask, intptr_t idx) {
	idx = bitmask_shuffle(idx);
	int offset = idx & BLOCK_1;
	BLOCK_TYPE bit = ((BLOCK_TYPE)1) << offset;
	__sync_and_and_fetch(&mask->data[idx >> LOG2_BLOCK], ~bit);
}

static inline void bitmask_clear_unsafe(bitmask_t * mask, intptr_t idx) {
	idx = bitmask_shuffle(idx);
	int offset = idx & BLOCK_1;
	BLOCK_TYPE bit = ((BLOCK_TYPE)1) << offset;
	mask->data[idx >> LOG2_BLOCK] &= ~bit;
}

static inline int bitmask_test_and_clear(bitmask_t * mask, intptr_t idx) {
	idx = bitmask_shuffle(idx);
	int offset = idx & BLOCK_1;
	BLOCK_TYPE bit = ((BLOCK_TYPE)1) << offset;
	return (__sync_fetch_and_and(&mask->data[idx >> LOG2_BLOCK], ~bit) & bit) != 0;
}

static inline int bitmask_test_and_set(bitmask_t * mask, intptr_t idx) {
	idx = bitmask_shuffle(idx);
	int offset = idx & BLOCK_1;
	BLOCK_TYPE bit = ((BLOCK_TYPE)1) << offset;
	return (__sync_fetch_and_or(&mask->data[idx >> LOG2_BLOCK], bit) & bit) != 0;
}

static inline int bitmask_aquire(bitmask_t * mask, intptr_t idx, intptr_t timeout) {
	idx = bitmask_shuffle(idx);
	int offset = idx & BLOCK_1;
	BLOCK_TYPE bit = ((BLOCK_TYPE)1) << offset;
	intptr_t c = 0;
	volatile BLOCK_TYPE * p = mask->data + (idx >> LOG2_BLOCK);

	while( *p & bit) {
		c++;
		if(c >= timeout) return 0;
	}
	return !(__sync_fetch_and_or(p, bit) & bit);
}

static inline int bitmask_test(bitmask_t * mask, intptr_t idx) {
	idx = bitmask_shuffle(idx);
	int offset = idx & BLOCK_1;
	BLOCK_TYPE bit = ((BLOCK_TYPE)1) << offset;
	return (mask->data[idx >> LOG2_BLOCK] & bit) != 0;
}

static inline size_t bitmask_sum(bitmask_t * mask) {
	size_t size = *(size_t *) mask;
	size_t sum = 0;
	intptr_t i;
	for(i = 0; i < mask->bytes; i++) {
		sum += __builtin_popcount((BLOCK_TYPE) mask->data[i]);
	}
	return sum;
}
#undef BLOCK_TYPE
#undef LOG2_BLOCK
#undef BLOCK
#undef BLOCK_1

