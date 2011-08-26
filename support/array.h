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
	(array ## _size < size) &&  \
		(array ## _size = __array_roundup__(array ## _size, size), \
		array = realloc(array, sizeof(type) * array ## _size)), \
	array)

#define ARRAY_ENSURE0(array, type, size) \
	((array == NULL) && \
		(array ## _size = 128, \
		array = calloc(sizeof(type), array ## _size)), \
	(array ## _size < size) &&  \
		(array = __array_realloc0__(array, array ## _size, size, sizeof(type)), array ## _size = __array_roundup__(array ## _size, size)), \
	array)

#define ARRAY_RESIZE(array, type, length) \
	(ARRAY_ENSURE0(array, type, length), array ## _length = length, array)

#define ARRAY_APPEND(array, type) \
	(array ## _length = array ## _length + 1, \
	&(ARRAY_ENSURE(array, type, array ## _length)[array ## _length - 1]))

#define ARRAY_FREE(array) \
	(free(array), array ## _length = 0, array ## _size = 0, array = NULL, array)
static inline size_t __array_roundup__(size_t old, size_t new_) {
	while(old < new_) old <<= 1;
	return old;
}
static inline void * __array_realloc0__(void * pt, size_t old, size_t new_, size_t unit) {
	size_t new__ = __array_roundup__(old, new_);
	void * rt = realloc(pt, new__ * unit);
	memset((char*)rt + old * unit, 0, (new__ - old) * unit);
	return rt;
}


