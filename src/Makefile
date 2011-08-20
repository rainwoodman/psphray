CC = gcc
CFLAGS = -std=gnu99 -D_GNU_SOURCE -g -mfpmath=387,sse -ffast-math -fopenmp -I ../support
OPT=-O3
F90FLAGS = -g -ffast-math -O3 
FC = gfortran
LDFLAGS = -pg -L ../support/libconfig/
LIBS = -lconfig -lgfortran -lgsl -lgslcblas -lgomp

SOURCES = \
config.c \
units.c \
epochs.c \
reader.c \
psystem.c \
raytrace.c \
spectra.c \
solve.c \
atomic-rate.c \
loop.c \
sph-depth.c \
reader-massiveblack.c \
reader-psphray.c

OBJS = $(SOURCES:.c=.o) pluecker.o verner.o
ALL = psphray sane test-solve

.SUFFIXES: .inc .F90

all: $(ALL)

.c.o:
	$(CC) $(CFLAGS) $(OPT) -c -o $@ $^
#	$(CC) $(CFLAGS) $(OPT) -fverbose-asm -S $^

sane.o: sane.c
	$(CC) $(CFLAGS) -c -o $@ $^

.F90.o:
	$(FC) $(F90FLAGS) -c -o $@ $<

psphray: $(OBJS) main.o
	$(CC) $(LDFLAGS) -o $@  $^ $(LIBS)

sane: $(OBJS) sane.o
	$(CC) $(LDFLAGS) -o $@  $^ $(LIBS)

test-solve: $(OBJS) test-solve.o
	$(CC) $(LDFLAGS) -o $@  $^ $(LIBS)

clean:
	rm $(ALL) *.o

.inc.c: reader-internal.h
	touch $@

