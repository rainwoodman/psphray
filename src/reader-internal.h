struct header_t;
typedef struct {
	int (*blockid)(char * blk, char ** error);
	char * (*blockname)(int blockid, char ** error);
	size_t (*itemsize)(char * blk, char ** error);
	size_t (*pstart)(struct header_t * h, char * blk, int ptype, char ** error);
	size_t (*length)(struct header_t * h, char * blk, char ** error);
	size_t (*npar)(struct header_t * h, int ptype, char ** error);
	size_t (*offset)(struct header_t * h, char * blk, char ** error);
	void (*read)(struct header_t * h, char * blk, void * buffer, int start, int length, FILE * fp, char ** error);
	void (*write)(struct header_t * h, char * blk, void * buffer, int start, int length, FILE * fp, char ** error);
} ReaderFuncs;
