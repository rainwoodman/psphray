struct header_t;

typedef struct {
	int (*blockid)(const char * blk, char ** error);
	const char * (*blockname)(int blockid, char ** error);
	size_t (*itemsize)(const char * blk, char ** error);
	size_t (*pstart)(const struct header_t * h, const char * blk, int ptype, char ** error);
	size_t (*length)(const struct header_t * h, const char * blk, char ** error);
	void (*get_constants)(const struct header_t * h, ReaderConstants * c);
	void (*set_constants)(struct header_t * h, const ReaderConstants * c);
	void (*def_header)(struct header_t * h);
	size_t (*offset)(const struct header_t * h, const char * blk, char ** error);
	void (*read)(struct header_t * h, const char * blk, void * buffer, int start, int length, FILE * fp, char ** error);
	void (*write)(const struct header_t * h, const char * blk, const void * buffer, int start, int length, FILE * fp, char ** error);
} ReaderFuncs;

