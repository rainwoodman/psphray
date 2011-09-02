
	float pos[3] = {0., 0., 0.};
	float dir[3] = {0.688749, 0.722979, 0.0540979};
	float dist = 4.632;
	intptr_t * ipars = NULL;
	size_t size =0;
	size_t length = rt_trace(pos, dir, dist, &ipars, &size);
	for(i = 0; i < length; i++) {
		intptr_t ipar = ipars[i];
		printf("%ld %g %g %g\n", ipar, psys.pos[ipar][0], psys.pos[ipar][1], psys.pos[ipar][2]);
	}
	free(ipars);
