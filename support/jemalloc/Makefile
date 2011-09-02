
SOURCES=\
arena.c \
atomic.c \
base.c \
bitmap.c \
chunk.c \
chunk_dss.c \
chunk_mmap.c \
chunk_swap.c \
ckh.c \
ctl.c \
extent.c \
hash.c \
huge.c \
jemalloc.c \
mb.c \
mutex.c \
prof.c \
rtree.c \
stats.c \
tcache.c

all: libjemalloc.a
.c.o:
	$(CC) -D_GNU_SOURCE -I.. -I. -I./internal -O3 -c -o $@ $<

libjemalloc.a: $(SOURCES:.c=.o)
	ar cru $@ $^
clean:
	$(RM) *.[oa]
