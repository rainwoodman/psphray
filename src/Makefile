CFLAGS = -g -I ../support
FC = gfortran
LDFLAGS = -L ../support/libconfig/ 
LIBS = -lconfig -lm -lgfortran

SOURCES = \
config.c \
units.c \
epochs.c \
reader.c \
psystem.c \
raytrace.c \
reader-massiveblack.c \

OBJS = $(SOURCES:.c=.o) pluecker.o
ALL = sane

.SUFFIXES: .inc

all: $(ALL)

.c.o:
	$(CC) $(CFLAGS) -c -o $@ $^

pluecker.o: pluecker.F90
	$(FC) -g -c -o $@ $<

psphray: $(OBJS) main.o
	$(CC) $(LDFLAGS) -o $@  $^ $(LIBS)

sane: $(OBJS) sane.o
	$(CC) $(LDFLAGS) -o $@  $^ $(LIBS)

clean:
	rm $(ALL) *.o

.inc.c: reader-internal.h
	touch $@

