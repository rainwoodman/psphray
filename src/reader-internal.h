struct header_t;

typedef struct {
	int (*blockid)(char * blk, char ** error);
	char * (*blockname)(int blockid, char ** error);
	size_t (*itemsize)(char * blk, char ** error);
	size_t (*pstart)(struct header_t * h, char * blk, int ptype, char ** error);
	size_t (*length)(struct header_t * h, char * blk, char ** error);
	void (*get_constants)(struct header_t * h, ReaderConstants * c);
	void (*set_constants)(struct header_t * h, ReaderConstants * c);
	size_t (*npar_total)(struct header_t * h, int ptype);
	double (*time)(struct header_t * h);
	double (*boxsize)(struct header_t * h);
	void (*def_header)(struct header_t * h);
	size_t (*offset)(struct header_t * h, char * blk, char ** error);
	void (*read)(struct header_t * h, char * blk, void * buffer, int start, int length, FILE * fp, char ** error);
	void (*write)(struct header_t * h, char * blk, void * buffer, int start, int length, FILE * fp, char ** error);
} ReaderFuncs;

