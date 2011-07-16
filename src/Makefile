CFLAGS = -g -I ../support

LDFLAGS = -L ../support/libconfig/ 
LIBS = -lconfig

SOURCES = \
config.c \
main.c \
reader.c \
reader-massiveblack.c

ALL = psphray

.SUFFIXES: .inc

all: $(ALL)

.c.o:
	$(CC) $(CFLAGS) -c -o $@ $^

psphray: $(SOURCES:.c=.o)
	$(CC) $(LDFLAGS) -o $@  $^ $(LIBS)

clean:
	rm $(ALL) *.o

.inc.c: reader-internal.h
	touch $@

