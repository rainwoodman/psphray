#include <stdio.h>
#include <stdint.h>
#include <messages.h>
#include "reader.h"

extern void cfg_init(char * filename);
extern void cfg_dump(char * filename);

int main(int argc, char* argv[]) {
	MESSAGE("hello, world, %s\n", argv[1]);
	cfg_init(argv[1]);
	cfg_dump("used-config");
	Reader * r = reader_new("massiveblack");
	reader_open(r, "../test/snapshot_036.0");
	float (*pos)[3] = reader_alloc(r, "pos", -1);
	size_t l = reader_length(r, "pos");

	reader_read(r, "pos", -1, pos);
	intptr_t i;
	float x, y, z;
	for(i = 0; i < l; i++) {
		x += pos[i][0];
		y += pos[i][1];
		z += pos[i][2];
	}
	printf("%f %f %f\n", x/l, y/l, z/l);
	return 0;
}
