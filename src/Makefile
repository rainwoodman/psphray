CFLAGS = -g -I ../support

LDFLAGS = -L ../support/libconfig/ 
LIBS = -lconfig

SOURCES = \
config.c \
units.c \
reader.c \
reader-massiveblack.c

ALL = sane

.SUFFIXES: .inc

all: $(ALL)

.c.o:
	$(CC) $(CFLAGS) -c -o $@ $^

psphray: $(SOURCES:.c=.o) main.c
	$(CC) $(LDFLAGS) -o $@  $^ $(LIBS)

sane: $(SOURCES:.c=.o) sane.o
	$(CC) $(LDFLAGS) -o $@  $^ $(LIBS)

clean:
	rm $(ALL) *.o

.inc.c: reader-internal.h
	touch $@

