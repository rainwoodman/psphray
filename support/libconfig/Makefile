SOURCES= grammar.c \
libconfig.c \
scanctx.c \
scanner.c \
strbuf.c

all: libconfig.a

.c.o:
	$(CC) -g -c -o $@ $<
libconfig.a: $(SOURCES:.c=.o)
	$(AR) cru $@ $^

clean:
	$(RM) *.[oa]
