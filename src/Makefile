CC = gcc
CFLAGS = -std=gnu99 -D_GNU_SOURCE -g -ffast-math -I ../support 
OPT= -O3 -fopenmp
F90FLAGS = -ffast-math -O3 
FC = ifort
LDFLAGS = -L ../support/libconfig/ -L ../support/jemalloc 
LIBS = -lconfig -lgfortran -lgsl -lgslcblas -lguide -liomp5 -limf -lpthread
LIBS += -ljemalloc
SOURCES = \
config.c \
units.c \
epochs.c \
reader.c \
psystem.c \
raytrace.c \
spectra.c \
solve.c \
peano.c \
lte.c \
atomic-rate.c \
loop.c \
sph-depth.c \
reader-massiveblack.c \
reader-psphray.c

OBJS = $(SOURCES:.c=.o) pluecker-sse.o verner.o
ALL = psphray sane test-solve test-peano

.SUFFIXES: .inc .F90 .s

all: $(ALL)
asm: $(SOURCES:.c=.s) pluecker.s pluecker-sse.s verner.s

.c.o:
	$(CC) $(CFLAGS) $(OPT) -c -o $@ $^
.c.s:
	$(CC) $(CFLAGS) $(OPT) -fverbose-asm -S $^
.F90.s:
	$(FC) $(F90FLAGS) -fverbose-asm -S $^

sane.o: sane.c
	$(CC) $(CFLAGS) -c -o $@ $^

.F90.o:
	$(FC) $(F90FLAGS) -c -o $@ $<

psphray: $(OBJS) main.o
	$(CC) $(LDFLAGS) -o $@  $^ $(LIBS)

sane: $(OBJS) sane.o
	$(CC) $(LDFLAGS) -o $@  $^ $(LIBS)

test-peano: $(OBJS) test-peano.o
	$(CC) $(LDFLAGS) -o $@  $^ $(LIBS)

test-solve: $(OBJS) test-solve.o
	$(CC) $(LDFLAGS) -o $@  $^ $(LIBS)

test-bitmask: $(OBJS) test-bitmask.o
	$(CC) -g -fopenmp -O3 $(LDFLAGS) -o $@  $^ $(LIBS)

clean:
	rm $(ALL) *.o

.inc.c: reader-internal.h
	touch $@

