CFLAGS = -g -ffast-math -fopenmp -I ../support
OPT=-O3
F90FLAGS = -g -ffast-math -O3 
FC = gfortran
LDFLAGS = -pg -L ../support/libconfig/
LIBS = -lconfig -lm -lgfortran -lgsl -lgslcblas -lgomp

SOURCES = \
config.c \
units.c \
epochs.c \
reader.c \
psystem.c \
raytrace.c \
solve.c \
atomic-rate.c \
loop.c \
sph-depth.c \
reader-massiveblack.c \
reader-psphray.c

OBJS = $(SOURCES:.c=.o) pluecker.o verner.o
ALL = sane test-solve

.SUFFIXES: .inc .F90

all: $(ALL)

.c.o:
	$(CC) $(OPT) $(CFLAGS) -c -o $@ $^

sane.o: sane.c
	$(CC) -O3 $(CFLAGS) -c -o $@ $^

loop.o: loop.c
	$(CC) -O3 $(CFLAGS) -c -o $@ $^

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

