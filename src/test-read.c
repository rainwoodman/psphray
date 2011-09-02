
	Reader * r = reader_new("massiveblack");
	reader_open(r, "../test/snapdir_680/snapshot_680.0");
	ReaderConstants * c = reader_constants(r);
	MESSAGE("-------reader sanity check--------");
	MESSAGE("boxsize = %f", c->boxsize);
	MESSAGE("redshift = %f", c->redshift);
	MESSAGE("Ntot[0] = %lu", c->Ntot[0]);
	MESSAGE("N[0] = %lu", c->N[0]);
	MESSAGE("Nfiles = %lu", c->Nfiles);
	MESSAGE("h = %f", c->h);
	MESSAGE("-------reader sanity check--------");

	float *ie = reader_alloc(r, "ie", 0);
	float *ye = reader_alloc(r, "ye", 0);
	reader_read(r, "ie", 0, ie);
	reader_read(r, "ye", 0, ye);

	float Tmin = 1000000000;
	float Tmax = -1;
	for(i = 0; i < reader_npar(r, 0); i++) {
		float T = ieye2T(ie[i], ye[i]);
		if(T > Tmax) Tmax = T;
		if(T < Tmin) Tmin = T;
	}
	free(ie);
	free(ye);
	reader_destroy(r);
	printf("%f %f\n", Tmax, Tmin);

