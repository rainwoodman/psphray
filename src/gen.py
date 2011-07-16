from sys import argv

_temp = __import__('gaepsi.readers.%s' % argv[1], globals(), locals(),
	  ['Reader'],  -1)
reader = _temp.Reader()

from gaepsi.tools.cgen import gen_reader

f = file("reader-%s.inc" % argv[1], 'w')
f.write(gen_reader(reader))
f.close()
f = file("reader-%s.c" % argv[1], 'w')
f.write(
"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "reader-internal.h"
#define constants_t _ReaderConstants
#include "reader-%s.inc"

void %s_get_reader_funcs(ReaderFuncs * funcs) {
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
""" % (argv[1], argv[1]))
f.close()
