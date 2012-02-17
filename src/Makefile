CC = mpigcc
LD = mpigcc
CFLAGS = -std=gnu99 -D_GNU_SOURCE -g -ffast-math -I ../support -DUSE_LSODA
#OPT= -O3 -fopenmp
F90FLAGS = -ffast-math -O3 
FC = ifort
LDFLAGS = -L ../support/libconfig/ -L ../support/jemalloc -L ../support/lsoda
#LIBS += -ljemalloc
LIBS = -lconfig -lgfortran -lgsl -lgslcblas -llsoda -limf -liomp5 -lpthread
include Makefile.local

SOURCES = \
config.c \
tabfun.c \
units.c \
stat.c \
epochs.c \
reader.c \
psystem.c \
sources.c \
raytrace.c \
spectra.c \
solve.c \
peano.c \
secion.c \
lte.c \
eos.c \
atomic-rate.c \
loop.c \
sph-depth.c \
discrete.c \
ompdiscrete.c \
reader-massiveblack.c \
reader-psphray.c

OBJS = $(SOURCES:.c=.o) pluecker-sse.o verner.o
ALL = psphray sane test-solve test-peano

.SUFFIXES: .inc .F90 .s

all: $(ALL)
asm: $(SOURCES:.c=.s) pluecker.s pluecker-sse.s verner.s

Makefile.local:
	touch $@

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
	$(LD) $(LDFLAGS) -o $@  $^ $(LIBS)

balance: balance.o mortonkey.o psystem.o
	$(LD) $(LDFLAGS) -o $@  $^ $(LIBS)

sane: $(OBJS) sane.o
	$(LD) $(LDFLAGS) -o $@  $^ $(LIBS)

test-peano: $(OBJS) test-peano.o
	$(LD) $(LDFLAGS) -o $@  $^ $(LIBS)

test-solve: $(OBJS) test-solve.o
	$(LD) $(LDFLAGS) -o $@  $^ $(LIBS)

test-bitmask: $(OBJS) test-bitmask.o
	$(LD) -g -fopenmp -O3 $(LDFLAGS) -o $@  $^ $(LIBS)

test-lte: $(OBJS) test-lte.o
	$(LD) -g -fopenmp -O3 $(LDFLAGS) -o $@  $^ $(LIBS)

test-spectra: $(OBJS) test-spectra.o
	$(LD) -g -fopenmp -O3 $(LDFLAGS) -o $@  $^ $(LIBS)

clean:
	rm $(ALL) *.o

.inc.c: reader-internal.h
	touch $@

