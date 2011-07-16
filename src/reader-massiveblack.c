
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "reader-internal.h"
#include "reader-massiveblack.inc"

void massiveblack_get_reader_funcs(ReaderFuncs * funcs) {
	funcs->blockid = _blockid;
	funcs->blockname = _blockname;
	funcs->itemsize = _itemsize;
	funcs->pstart = _pstart;
	funcs->length = _length;
	funcs->npar = _npar;
	funcs->offset = _offset;
	funcs->read = _read;
	funcs->write = _write;
};
