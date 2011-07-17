#define MESSAGE(fmt,...) { fprintf(stderr, "MM:" fmt " %s:%d\n", \
						## __VA_ARGS__, __FILE__, __LINE__); fflush(stderr); }

#define WARNING(fmt,...) { fprintf(stderr, "WW:" fmt " %s:%d\n", \
						## __VA_ARGS__, __FILE__, __LINE__); fflush(stderr); }
#ifndef MPI_INCLUDED
#define ERROR(fmt,...) { fprintf(stderr, "EE:" fmt " %s:%d\n",\
					 ## __VA_ARGS__, __FILE__, __LINE__); fflush(stderr); abort();}
#else
#define ERROR(fmt,...) { fprintf(stderr, "EE:" fmt " %s:%d\n",\
					 ## __VA_ARGS__, __FILE__, __LINE__); fflush(stderr); MPI_Abort(MPI_COMM_WORLD, -1);}
#endif
#define ERROR_R(code, fmt,...) { fprintf(stderr, "EE:" fmt " %s:%d\n",\
					 ## __VA_ARGS__, __FILE__, __LINE__); fflush(stderr); return(code);}

