#include <stdio.h>
#include <stdlib.h>
#include <signal.h>

#define MESSAGE(fmt,...) ( fprintf(stderr, "MM:" fmt " @(%s:%d)\n", \
						## __VA_ARGS__, __FILE__, __LINE__), fflush(stderr) )

#define WARNING(fmt,...) ( fprintf(stderr, "WW:" fmt " @(%s:%d)\n", \
						## __VA_ARGS__, __FILE__, __LINE__), fflush(stderr) )

#define LOG_ERROR(fmt,...) ( fprintf(stderr, "EE:" fmt " @(%s:%d)\n",\
				 ## __VA_ARGS__, __FILE__, __LINE__), fflush(stderr))

#ifndef MPI_INCLUDED
#define ERROR(fmt,...) ( LOG_ERROR(fmt, ## __VA_ARGS__ ), raise(SIGTRAP))
#else
#define ERROR(fmt,...) ( LOG_ERROR(fmt, ## __VA_ARGS__ ), MPI_Abort(MPI_COMM_WORLD, -1))
#endif

