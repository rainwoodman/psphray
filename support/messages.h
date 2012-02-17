#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <stdint.h>
#include <printf.h>
#define PRINTF_BINARY(c) register_printf_function((c), printf_binary, printf_binary_arginfo)
#define MESSAGE(fmt,...) (fprintf(stderr, "MM:" fmt " @(%s:%d)\n", \
						## __VA_ARGS__, __FILE__, __LINE__), fflush(stderr) )

#define WARNING(fmt,...) (fprintf(stderr, "WW:" fmt " @(%s:%d)\n", \
						## __VA_ARGS__, __FILE__, __LINE__), fflush(stderr) )

#define LOG_ERROR(fmt,...) (fprintf(stderr, "EE:" fmt " @(%s:%d)\n",\
				 ## __VA_ARGS__, __FILE__, __LINE__), fflush(stderr))

#ifndef MPI_INCLUDED
#define ERROR(fmt,...) ( LOG_ERROR(fmt, ## __VA_ARGS__ ), raise(SIGTRAP))
#else
#define ERROR(fmt,...) ( LOG_ERROR(fmt, ## __VA_ARGS__ ), MPI_Abort(MPI_COMM_WORLD, -1))
#endif

static int printf_binary(FILE * stream, const struct printf_info * info, const void * const *args);
static int printf_binary_arginfo (const struct printf_info *info, size_t n, int *argtypes);

static int printf_binary(FILE * stream, const struct printf_info * info, const void * const *args) {
	const int bits = info->width;
	int buflen = bits + 1;
	if(info->prec != -1) {
		buflen == bits / info->prec + 1;
	}
	char buf[buflen];
	int p = 0;
	unsigned int mask = (1 << (bits - 1));
	intptr_t l;
	if(info->is_long) {
		l = *(long * )(args[0]);
	} else if(info->is_short) {
		l = *(short *)(args[0]);
	} else if(info->is_char) {
		l = *(char *)(args[0]);
	} else {
		l = *(int *)(args[0]);
	}
	int j = 0;
	while(p < bits) {
		buf[j] = (l & mask)?'1':'0';
		mask >>= 1;
		p++;
		j++;
		if((info->prec != -1) && ((bits - p) % info->prec == 0)) {
			buf[j] = ' ';
			j++;
		}
	}
	buf[j] = 0;
	fputs(buf, stream);
	return bits;
}
static int printf_binary_arginfo (const struct printf_info *info, size_t n, int *argtypes) {
	if (n > 0)
		argtypes[0] = PA_POINTER;
	return 1;
}

