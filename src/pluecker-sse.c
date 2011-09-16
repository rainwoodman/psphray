
typedef float v4sf __attribute__ ((vector_size(16)));
typedef union {
	float a[4];
	v4sf v;
} ssef;
typedef union {
	short int a[2];
	int b;
} cmpr;

static inline ssef prodsub(const float f1, const float f2, const float f3, const float f4,
				const float g1, const float g2, const float g3, const float g4) {
	register ssef f;
	f.a[0] = f1;
	f.a[1] = f2;
	f.a[2] = f3;
	f.a[3] = f4;
	register ssef g;
	g.a[0] = g1;
	g.a[1] = g2;
	g.a[2] = g3;
	g.a[3] = g4;
	f.v = f.v * g.v;
	return f;
}

int pluecker_(const float dir[3], const float * dist, const float const s2b[3], const float const s2t[3]) {
	const int class = ((dir[0] > 0)<< 0) + ((dir[1] > 0)<< 1) + ((dir[2] > 0) << 2);

	float e2t[3];
	float e2b[3];
	int i;
	for(i = 0; i < 3; i++) {
		e2t[i] = s2t[i] - *dist * dir[i];
	}
	for(i = 0; i < 3; i++) {
		e2b[i] = s2b[i] - *dist * dir[i];
	}
	ssef r;
	switch(class) {
		case 0:
			if(s2b[0] > 0 || s2b[1] > 0 || s2b[2] > 0) return 0;
			if(e2t[0] < 0 || e2t[1] < 0 || e2t[2] < 0) return 0;

			r = prodsub(
				dir[0], dir[1], dir[0], dir[1],
				s2b[1], s2t[0], s2t[1], s2b[0]);
			if(r.a[0] < r.a[1] || r.a[2] > r.a[3]) return 0;

			r = prodsub(
				dir[0], dir[2], dir[0], dir[2],
				s2b[2], s2t[0], s2t[2], s2b[0]);
			if(r.a[0] < r.a[1] || r.a[2] > r.a[3]) return 0;

			r = prodsub(
				dir[1], dir[2], dir[1], dir[2],
				s2b[2], s2t[1], s2t[2], s2b[1]);
			if(r.a[0] < r.a[1] || r.a[2] > r.a[3]) return 0;

		break;
		case 1:
			if(s2b[0] < 0 || s2b[1] > 0 || s2b[2] > 0) return 0;
			if(e2t[0] > 0 || e2t[1] < 0 || e2t[2] < 0) return 0;

			r = prodsub(
				dir[0], dir[1], dir[0], dir[1],
				s2t[1], s2t[0], s2b[1], s2b[0]);
			if(r.a[0] < r.a[1] || r.a[2] > r.a[3]) return 0;

			r = prodsub(
				dir[0], dir[2], dir[0], dir[2],
				s2t[2], s2t[0], s2b[2], s2b[0]);
			if(r.a[0] < r.a[1] || r.a[2] > r.a[3]) return 0;

			r = prodsub(
				dir[1], dir[2], dir[1], dir[2],
				s2b[2], s2t[1], s2t[2], s2b[1]);
			if(r.a[0] < r.a[1] || r.a[2] > r.a[3]) return 0;
		break;
		case 2:
			if(s2b[0] < 0 || s2t[1] < 0 || s2b[2] > 0) return 0;
			if(e2t[0] > 0 || e2b[1] > 0 || e2t[2] < 0) return 0;

			r = prodsub(
				dir[0], dir[1], dir[0], dir[1],
				s2b[1], s2b[0], s2t[1], s2t[0]);
			if(r.a[0] < r.a[1] || r.a[2] > r.a[3]) return 0;

			r = prodsub(
				dir[0], dir[2], dir[0], dir[2],
				s2b[2], s2t[0], s2t[2], s2b[0]);
			if(r.a[0] < r.a[1] || r.a[2] > r.a[3]) return 0;

			r = prodsub(
				dir[1], dir[2], dir[1], dir[2],
				s2t[2], s2t[1], s2b[2], s2b[1]);
			if(r.a[0] < r.a[1] || r.a[2] > r.a[3]) return 0;
		break;
		case 3:
			if(s2t[0] < 0 || s2t[1] < 0 || s2b[2] > 0) return 0;
			if(e2b[0] > 0 || e2b[1] > 0 || e2t[2] < 0) return 0;

			r = prodsub(
				dir[0], dir[1], dir[0], dir[1],
				s2b[1], s2b[0], s2t[1], s2t[0]);
			if(r.a[0] < r.a[1] || r.a[2] > r.a[3]) return 0;

			r = prodsub(
				dir[0], dir[2], dir[0], dir[2],
				s2b[2], s2t[0], s2t[2], s2b[0]);
			if(r.a[0] < r.a[1] || r.a[2] > r.a[3]) return 0;

			r = prodsub(
				dir[1], dir[2], dir[1], dir[2],
				s2t[2], s2t[1], s2b[2], s2b[1]);
			if(r.a[0] < r.a[1] || r.a[2] > r.a[3]) return 0;
		break;
		case 4:
			if(s2b[0] > 0 || s2b[1] > 0 || s2t[2] < 0) return 0;
			if(e2t[0] < 0 || e2t[1] < 0 || e2b[2] > 0) return 0;

			r = prodsub(
				dir[0], dir[1], dir[0], dir[1],
				s2b[1], s2t[0], s2t[1], s2b[0]);
			if(r.a[0] < r.a[1] || r.a[2] > r.a[3]) return 0;

			r = prodsub(
				dir[0], dir[2], dir[0], dir[2],
				s2b[2], s2b[0], s2t[2], s2t[0]);
			if(r.a[0] < r.a[1] || r.a[2] > r.a[3]) return 0;

			r = prodsub(
				dir[1], dir[2], dir[1], dir[2],
				s2b[2], s2b[1], s2t[2], s2t[1]);
			if(r.a[0] < r.a[1] || r.a[2] > r.a[3]) return 0;
		break;
		case 5:
			if(s2t[0] < 0 || s2b[1] > 0 || s2t[2] < 0) return 0;
			if(e2b[0] > 0 || e2t[1] < 0 || e2b[2] > 0) return 0;

			r = prodsub(
				dir[0], dir[1], dir[0], dir[1],
				s2t[1], s2t[0], s2b[1], s2b[0]);
			if(r.a[0] < r.a[1] || r.a[2] > r.a[3]) return 0;

			r = prodsub(
				dir[0], dir[2], dir[0], dir[2],
				s2t[2], s2b[0], s2b[2], s2t[0]);
			if(r.a[0] < r.a[1] || r.a[2] > r.a[3]) return 0;

			r = prodsub(
				dir[1], dir[2], dir[1], dir[2],
				s2b[2], s2b[1], s2t[2], s2t[1]);
			if(r.a[0] < r.a[1] || r.a[2] > r.a[3]) return 0;
		break;
		case 6:
			if(s2b[0] > 0 || s2t[1] < 0 || s2t[2] < 0) return 0;
			if(e2t[0] < 0 || e2b[1] > 0 || e2b[2] > 0) return 0;

			r = prodsub(
				dir[0], dir[1], dir[0], dir[1],
				s2b[1], s2b[0], s2t[1], s2t[0]);
			if(r.a[0] < r.a[1] || r.a[2] > r.a[3]) return 0;

			r = prodsub(
				dir[0], dir[2], dir[0], dir[2],
				s2b[2], s2b[0], s2t[2], s2t[0]);
			if(r.a[0] < r.a[1] || r.a[2] > r.a[3]) return 0;

			r = prodsub(
				dir[1], dir[2], dir[1], dir[2],
				s2t[2], s2b[1], s2b[2], s2t[1]);
			if(r.a[0] < r.a[1] || r.a[2] > r.a[3]) return 0;
		break;
		case 7:
			if(s2t[0] < 0 || s2t[1] < 0 || s2t[2] < 0) return 0;
			if(e2b[0] > 0 || e2b[1] > 0 || e2b[2] > 0) return 0;

			r = prodsub(
				dir[0], dir[1], dir[0], dir[1],
				s2t[1], s2b[0], s2b[1], s2t[0]);
			if(r.a[0] < r.a[1] || r.a[2] > r.a[3]) return 0;

			r = prodsub(
				dir[0], dir[2], dir[0], dir[2],
				s2t[2], s2b[0], s2b[2], s2t[0]);
			if(r.a[0] < r.a[1] || r.a[2] > r.a[3]) return 0;

			r = prodsub(
				dir[1], dir[2], dir[1], dir[2],
				s2t[2], s2b[1], s2b[2], s2t[1]);
			if(r.a[0] < r.a[1] || r.a[2] > r.a[3]) return 0;
		break;

	}
	return 1;
}
