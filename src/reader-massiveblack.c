
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "reader.h"
#include "reader-internal.h"
#define constants_t _ReaderConstants
#include "reader-massiveblack.inc"

void massiveblack_get_reader_funcs(ReaderFuncs * funcs) {
	funcs->blockid = _blockid;
	funcs->blockname = _blockname;
	funcs->itemsize = _itemsize;
	funcs->pstart = _pstart;
	funcs->length = _length;
	funcs->get_constants = _get_constants;
	funcs->set_constants = _set_constants;
	funcs->def_header = _def_header;
	funcs->offset = _offset;
	funcs->read = _read;
	funcs->write = _write;
};
