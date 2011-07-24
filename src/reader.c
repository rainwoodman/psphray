#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <messages.h>
#include "reader.h"
#include "reader-internal.h"

struct header_t;
struct _Reader {
	char * filename;
/* private */
	struct header_t * header;
	FILE * fp;
	int read_only;
	ReaderFuncs funcs;
	ReaderConstants constants;
};

#define _blockid (reader->funcs.blockid)
#define _blockname (reader->funcs.blockname)
#define _itemsize (reader->funcs.itemsize)
#define _pstart (reader->funcs.pstart)
#define _length (reader->funcs.length)
#define _get_constants (reader->funcs.get_constants)
#define _set_constants (reader->funcs.set_constants)
#define _time (reader->funcs.time)
#define _boxsize (reader->funcs.boxsize)
#define _def_header (reader->funcs.def_header)
#define _offset (reader->funcs.offset)
#define _read (reader->funcs.read)
#define _write (reader->funcs.write)

void massiveblack_get_reader_funcs(ReaderFuncs * funcs);
void psphray_get_reader_funcs(ReaderFuncs * funcs);

Reader * reader_new(const char * format) {
	Reader * reader = calloc(sizeof(Reader), 1);
	if(!strcmp(format, "massiveblack")) {
		massiveblack_get_reader_funcs(&reader->funcs);
	} else if(!strcmp(format, "e5")) {
		massiveblack_get_reader_funcs(&reader->funcs);
	} else if(!strcmp(format, "psphray")) {
		psphray_get_reader_funcs(&reader->funcs);
	} else {
		ERROR("format %s unknown", format);
	}
	return reader;
}
void reader_destroy(Reader * reader) {
	if(reader->fp) reader_close(reader);
	free(reader->filename);
	free(reader);
}
char * reader_make_filename(const char * snapshot, int id) {
	if(strchr(snapshot, '%')) {
		char * t;
		asprintf(&t, snapshot, id);
		return t;
	}
	return strdup(snapshot);
}
void reader_open(Reader * reader, const char * filename) {
	reader->filename = strdup(filename);
	reader->fp = fopen(filename, "r");
	if(reader->fp == NULL) {
		ERROR("can't open %s", filename);
	}
	reader->read_only = 1;
	char * error = NULL;
	reader->header = reader_alloc(reader, "header", -1);
	_read(reader->header, "header", NULL, 0, 1, reader->fp, &error);
	if(error) {
		ERROR("%s", error);
		free(error);
	}
	_get_constants(reader->header, &reader->constants);
}
void reader_create(Reader * reader, const char * filename) {
	reader->filename = strdup(filename);
	reader->fp = fopen(filename, "w+");
	if(reader->fp == NULL) {
		ERROR("can't open %s", filename);
	}
	reader->header = reader_alloc(reader, "header", -1);
	_def_header(reader->header);
	reader->read_only = 0;
}

void reader_update_header(Reader * reader) {
	_set_constants(reader->header, &reader->constants);
}
void reader_update_constants(Reader * reader) {
	_get_constants(reader->header, &reader->constants);
}

void reader_close(Reader * reader) {
	char * error = NULL;
	if(!reader->read_only) {
		_write(reader->header, "header", NULL, 0, 1, reader->fp, &error);
	}
	free(reader->header);
	fclose(reader->fp);
	reader->fp = NULL;
}

size_t reader_length(Reader * reader, const char * blk) {
	char * error = NULL;
	size_t l = _length(reader->header, blk, &error);
	if(error) {
		ERROR("%s", error);
	}
	return l;
}

size_t reader_npar(Reader * reader, int ptype) {
	return reader->constants.N[ptype];
}

ReaderConstants * reader_constants(Reader * reader) {
	return &reader->constants;
}

void * reader_alloc(Reader * reader, const char * blk, int ptype) {
	char * error = NULL;
	size_t npar = 0;
	if(ptype != -1) {
		npar = reader_npar(reader, ptype);
	} else {
		npar = reader_length(reader, blk);
	}
	size_t b = _itemsize(blk, &error);
	if(error) {
		ERROR("%s", error);
		free(error);
	}
	return calloc(b, npar);
}

size_t reader_itemsize(Reader * reader, const char * blk) {
	char * error = NULL;
	size_t b = _itemsize(blk, &error);
	if(error) {
		ERROR("%s", error);
	}
	return b;
}

void reader_read(Reader * reader, const char * blk, int ptype, void *buf) {
	char * error = NULL;
	size_t pstart = 0;
	size_t npar = 0;
	if(ptype != -1) {
		pstart = _pstart(reader->header, blk, ptype, &error);
		if(error) {
			ERROR("%s", error);
			free(error);
		}
		npar = reader_npar(reader, ptype);
	} else {
		npar = reader_length(reader, blk);
	}
	_read(reader->header, blk, buf, pstart, npar, reader->fp, &error);
	if(error) {
		ERROR("%s", error);
		free(error);
	}
}

void reader_write(Reader * reader, const char * blk, int ptype, void *buf) {
	char * error = NULL;
	size_t pstart = 0;
	size_t npar = 0;
	if(ptype != -1) {
		pstart = _pstart(reader->header, blk, ptype, &error);
		npar = reader_npar(reader, ptype);
		if(error) {
			ERROR("%s", error);
			free(error);
		}
	} else {
		npar = reader_length(reader, blk);
	}
	_write(reader->header, blk, buf, pstart, npar, reader->fp, &error);
	if(error) {
		ERROR("%s", error);
		free(error);
	}
}

int reader_has(Reader * reader, const char * blk) {
	char * error = NULL;
	int id = _blockid(blk, &error);
	if(error) {
		free(error);
		return 0;
	}
	return 1;
}
