#define BLOCK_TYPE unsigned int

#define BLOCK (sizeof(BLOCK_TYPE) * 8)
static int __bitmask_log2[] = {
// 0 , 1, 2, 3,  4,  5,  6,  7, 8, 9, 10, ...
   -1, 3, 4, -1, 5, -1, -1, -1, 6, -1, -1, 
};

#define LOG2_BLOCK  5
//(__bitmask_log2[sizeof(BLOCK_TYPE)])
#define BLOCK_1 (BLOCK - 1)

static inline void * bitmask_alloc(size_t size) {
	size_t total = (size + BLOCK_1 ) >> LOG2_BLOCK;
	BLOCK_TYPE * rt = malloc(total * sizeof(BLOCK_TYPE) + sizeof(size_t));
	*(size_t *) rt = size;
	return rt;
}
static inline void bitmask_set_all(void * mask) {
	size_t size = *(size_t *) mask;
	size_t total = (size + BLOCK_1 ) >> LOG2_BLOCK;
	BLOCK_TYPE * buf = (BLOCK_TYPE*)(((size_t *) mask) + 1);
	memset(buf, -1, total * sizeof(BLOCK_TYPE));
	if((size & BLOCK_1) != 0) {
		BLOCK_TYPE m = ((((BLOCK_TYPE)1) << ((size & BLOCK_1))) -1);
		buf[total - 1] &= m;
	}
}

static inline void bitmask_clear_all(void * mask) {
	size_t size = *(size_t *) mask;
	size_t total = (size + BLOCK_1 ) >> LOG2_BLOCK;
	BLOCK_TYPE * buf = (BLOCK_TYPE*)(((size_t *) mask) + 1);
	memset(buf, 0, total * sizeof(BLOCK_TYPE));
}

static inline void bitmask_set(void * mask, intptr_t idx) {
	BLOCK_TYPE * buf = (BLOCK_TYPE*)(((size_t *) mask) + 1);
	int offset = idx & BLOCK_1;
	BLOCK_TYPE bit = ((BLOCK_TYPE)1) << offset;
	//buf[idx >> LOG2_BLOCK] |= bit;
	__sync_or_and_fetch(&buf[idx >> LOG2_BLOCK], bit);
}

static inline void bitmask_clear(void * mask, intptr_t idx) {
	BLOCK_TYPE * buf = (BLOCK_TYPE*)(((size_t *) mask) + 1);
	int offset = idx & BLOCK_1;
	BLOCK_TYPE bit = ((BLOCK_TYPE)1) << offset;
	//buf[idx >> LOG2_BLOCK] &= ~bit;
	__sync_and_and_fetch(&buf[idx >> LOG2_BLOCK], ~bit);
}

static inline int bitmask_test_and_clear(void * mask, intptr_t idx) {
	BLOCK_TYPE * buf = (BLOCK_TYPE*)(((size_t *) mask) + 1);
	int offset = idx & BLOCK_1;
	BLOCK_TYPE bit = ((BLOCK_TYPE)1) << offset;
//	return (buf[idx >> LOG2_BLOCK] & bit) != 0;
	return __sync_fetch_and_and(&buf[idx >> LOG2_BLOCK], ~bit) & bit;
}

static inline int bitmask_test_and_set(void * mask, intptr_t idx) {
	BLOCK_TYPE * buf = (BLOCK_TYPE*)(((size_t *) mask) + 1);
	int offset = idx & BLOCK_1;
	BLOCK_TYPE bit = ((BLOCK_TYPE)1) << offset;
//	return (buf[idx >> LOG2_BLOCK] & bit) != 0;
	return __sync_fetch_and_or(&buf[idx >> LOG2_BLOCK], bit) & bit;
}

static inline int bitmask_test(void * mask, intptr_t idx) {
	BLOCK_TYPE * buf = (BLOCK_TYPE*)(((size_t *) mask) + 1);
	int offset = idx & BLOCK_1;
	BLOCK_TYPE bit = ((BLOCK_TYPE)1) << offset;
	return (buf[idx >> LOG2_BLOCK] & bit) != 0;
}

static inline size_t bitmask_sum(void * mask) {
//int __builtin_popcount(unsigned int x);
	BLOCK_TYPE * buf = (BLOCK_TYPE*)(((size_t *) mask) + 1);
	size_t size = *(size_t *) mask;
	size_t sum = 0;
	intptr_t i;
	for(i = 0; i < (size + BLOCK_1)>> LOG2_BLOCK; i++) {
		sum += __builtin_popcount((BLOCK_TYPE) buf[i]);
	}
	return sum;
}
#undef BLOCK_TYPE
#undef LOG2_BLOCK
#undef BLOCK
#undef BLOCK_1

