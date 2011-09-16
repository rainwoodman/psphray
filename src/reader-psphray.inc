struct header_t{
unsigned int N[6];
double mass[6];
double time;
double redshift;
int flag_sfr;
int flag_feedback;
unsigned int Nparticle_total_low[6];
int flag_cool;
int Nfiles;
double boxsize;
double OmegaM;
double OmegaL;
double h;
int flag_sft;
int flag_met;
unsigned int Nparticle_total_high[6];
int flag_entropy;
double OmegaB;
unsigned long long nrays;
int flag_Hmf;
int flag_Hemf;
int flag_Helium;
int flag_gammaHI;
int flag_cloudy;
int flag_eos;
int unused[5];
};

static int _blockid(const char * blk, char ** error) {
	if(!strcmp(blk, "pos")) return 1 ;
	if(!strcmp(blk, "id")) return 2 ;
	if(!strcmp(blk, "mass")) return 3 ;
	if(!strcmp(blk, "ie")) return 4 ;
	if(!strcmp(blk, "rho")) return 5 ;
	if(!strcmp(blk, "ye")) return 6 ;
	if(!strcmp(blk, "xHI")) return 7 ;
	if(!strcmp(blk, "sml")) return 8 ;
	if(!strcmp(blk, "Hmf")) return 9 ;
	if(!strcmp(blk, "Hemf")) return 10 ;
	if(!strcmp(blk, "HeIa")) return 11 ;
	if(!strcmp(blk, "HeIIa")) return 12 ;
	if(!strcmp(blk, "gammaHI")) return 13 ;
	if(!strcmp(blk, "HIa_cloudy")) return 14 ;
	if(!strcmp(blk, "eos")) return 15 ;
	if(!strcmp(blk, "lasthit")) return 16 ;
	if(!strcmp(blk, "hits")) return 17 ;
	if(!strcmp(blk, "sft")) return 18 ;
	if(!strcmp(blk, "ngammas")) return 19 ;
	if(!strcmp(blk, "spec")) return 20 ;
	if(!strcmp(blk, "header")) return 0;
	asprintf(error, "block %s is unknown", blk);
	return -1;}

static const char * _blockname(const int blockid, char ** error) {
	static char * names[] = {"header", 
"pos","id","mass","ie","rho","ye","xHI","sml","Hmf","Hemf","HeIa","HeIIa","gammaHI","HIa_cloudy","eos","lasthit","hits","sft","ngammas","spec"
};

	if(blockid < 0 || blockid >= sizeof(names) / sizeof(char*)) {
		asprintf(error, "blockid %d is unknown", blockid);
		return NULL;
	}
	return names[blockid];
     
}

static size_t _itemsize(const char * blk, char ** error) {
	if(!strcmp(blk, "pos")) return 12 ;
	if(!strcmp(blk, "id")) return 8 ;
	if(!strcmp(blk, "mass")) return 4 ;
	if(!strcmp(blk, "ie")) return 4 ;
	if(!strcmp(blk, "rho")) return 4 ;
	if(!strcmp(blk, "ye")) return 4 ;
	if(!strcmp(blk, "xHI")) return 4 ;
	if(!strcmp(blk, "sml")) return 4 ;
	if(!strcmp(blk, "Hmf")) return 4 ;
	if(!strcmp(blk, "Hemf")) return 4 ;
	if(!strcmp(blk, "HeIa")) return 4 ;
	if(!strcmp(blk, "HeIIa")) return 4 ;
	if(!strcmp(blk, "gammaHI")) return 4 ;
	if(!strcmp(blk, "HIa_cloudy")) return 4 ;
	if(!strcmp(blk, "eos")) return 4 ;
	if(!strcmp(blk, "lasthit")) return 8 ;
	if(!strcmp(blk, "hits")) return 4 ;
	if(!strcmp(blk, "sft")) return 4 ;
	if(!strcmp(blk, "ngammas")) return 8 ;
	if(!strcmp(blk, "spec")) return 8 ;
	if(!strcmp(blk, "header")) return 256;
	asprintf(error, "block %s is unknown", blk);
	return (size_t)-1;}

static size_t _pstart(const struct header_t * h, const char * blk, int ptype, char ** error) {
	if(!strcmp(blk, "pos")) return (0 < ptype)?(size_t)(h->N[0]): 0+(4 < ptype)?(size_t)(h->N[4]): 0+(5 < ptype)?(size_t)(h->N[5]): 0 ;
	if(!strcmp(blk, "id")) return (0 < ptype)?(size_t)(h->N[0]): 0+(4 < ptype)?(size_t)(h->N[4]): 0+(5 < ptype)?(size_t)(h->N[5]): 0 ;
	if(!strcmp(blk, "mass")) return (0 < ptype)?(size_t)(h->N[0]): 0+(4 < ptype)?(size_t)(h->N[4]): 0+(5 < ptype)?(size_t)(h->N[5]): 0 ;
	if(!strcmp(blk, "ie")) return (0 < ptype)?(size_t)(h->N[0]): 0 ;
	if(!strcmp(blk, "rho")) return (0 < ptype)?(size_t)(h->N[0]): 0 ;
	if(!strcmp(blk, "ye")) return (0 < ptype)?(size_t)(h->N[0]): 0 ;
	if(!strcmp(blk, "xHI")) return (0 < ptype)?(size_t)(h->N[0]): 0 ;
	if(!strcmp(blk, "sml")) return (0 < ptype)?(size_t)(h->N[0]): 0 ;
	if(!strcmp(blk, "Hmf")) return (0 < ptype)?(size_t)(h->N[0]): 0 ;
	if(!strcmp(blk, "Hemf")) return (0 < ptype)?(size_t)(h->N[0]): 0 ;
	if(!strcmp(blk, "HeIa")) return (0 < ptype)?(size_t)(h->N[0]): 0 ;
	if(!strcmp(blk, "HeIIa")) return (0 < ptype)?(size_t)(h->N[0]): 0 ;
	if(!strcmp(blk, "gammaHI")) return (0 < ptype)?(size_t)(h->N[0]): 0 ;
	if(!strcmp(blk, "HIa_cloudy")) return (0 < ptype)?(size_t)(h->N[0]): 0 ;
	if(!strcmp(blk, "eos")) return (0 < ptype)?(size_t)(h->N[0]): 0 ;
	if(!strcmp(blk, "lasthit")) return (0 < ptype)?(size_t)(h->N[0]): 0 ;
	if(!strcmp(blk, "hits")) return (0 < ptype)?(size_t)(h->N[0]): 0 ;
	if(!strcmp(blk, "sft")) return (4 < ptype)?(size_t)(h->N[4]): 0 ;
	if(!strcmp(blk, "ngammas")) return (5 < ptype)?(size_t)(h->N[5]): 0 ;
	if(!strcmp(blk, "spec")) return (5 < ptype)?(size_t)(h->N[5]): 0 ;
	asprintf(error, "block %s is unknown/unstable for _pstart", blk);
	return (size_t)-1;}


static size_t _length(const struct header_t * h, const char * blk, char ** error) {


	if(!strcmp(blk, "pos")) {


		if(1) {


			return (size_t)(h->N[0])+(size_t)(h->N[4])+(size_t)(h->N[5]);


		} else {
			return 0;
		}
	}


	if(!strcmp(blk, "id")) {


		if(1) {


			return (size_t)(h->N[0])+(size_t)(h->N[4])+(size_t)(h->N[5]);


		} else {
			return 0;
		}
	}


	if(!strcmp(blk, "mass")) {


		if(1) {


			return (size_t)(h->N[0])+(size_t)(h->N[4])+(size_t)(h->N[5]);


		} else {
			return 0;
		}
	}


	if(!strcmp(blk, "ie")) {


		if(1) {


			return (size_t)(h->N[0]);


		} else {
			return 0;
		}
	}


	if(!strcmp(blk, "rho")) {


		if(1) {


			return (size_t)(h->N[0]);


		} else {
			return 0;
		}
	}


	if(!strcmp(blk, "ye")) {


		if(1) {


			return (size_t)(h->N[0]);


		} else {
			return 0;
		}
	}


	if(!strcmp(blk, "xHI")) {


		if(1) {


			return (size_t)(h->N[0]);


		} else {
			return 0;
		}
	}


	if(!strcmp(blk, "sml")) {


		if(1) {


			return (size_t)(h->N[0]);


		} else {
			return 0;
		}
	}


	if(!strcmp(blk, "Hmf")) {


		if(1&&h->flag_Hmf) {


			return (size_t)(h->N[0]);


		} else {
			return 0;
		}
	}


	if(!strcmp(blk, "Hemf")) {


		if(1&&h->flag_Hemf) {


			return (size_t)(h->N[0]);


		} else {
			return 0;
		}
	}


	if(!strcmp(blk, "HeIa")) {


		if(1&&h->flag_Helium) {


			return (size_t)(h->N[0]);


		} else {
			return 0;
		}
	}


	if(!strcmp(blk, "HeIIa")) {


		if(1&&h->flag_Helium) {


			return (size_t)(h->N[0]);


		} else {
			return 0;
		}
	}


	if(!strcmp(blk, "gammaHI")) {


		if(1&&h->flag_gammaHI) {


			return (size_t)(h->N[0]);


		} else {
			return 0;
		}
	}


	if(!strcmp(blk, "HIa_cloudy")) {


		if(1&&h->flag_cloudy) {


			return (size_t)(h->N[0]);


		} else {
			return 0;
		}
	}


	if(!strcmp(blk, "eos")) {


		if(1&&h->flag_eos) {


			return (size_t)(h->N[0]);


		} else {
			return 0;
		}
	}


	if(!strcmp(blk, "lasthit")) {


		if(1) {


			return (size_t)(h->N[0]);


		} else {
			return 0;
		}
	}


	if(!strcmp(blk, "hits")) {


		if(1) {


			return (size_t)(h->N[0]);


		} else {
			return 0;
		}
	}


	if(!strcmp(blk, "sft")) {


		if(1) {


			return (size_t)(h->N[4]);


		} else {
			return 0;
		}
	}


	if(!strcmp(blk, "ngammas")) {


		if(1) {


			return (size_t)(h->N[5]);


		} else {
			return 0;
		}
	}


	if(!strcmp(blk, "spec")) {


		if(1) {


			return (size_t)(h->N[5]);


		} else {
			return 0;
		}
	}


	if(!strcmp(blk, "header")) return 1;
	asprintf(error, "block %s is unknown", blk);
	return (size_t)-1;
}



static size_t _offset(const struct header_t * h, const char * blk, char ** error) {
	int id = _blockid(blk, error);
	int i;
	if(*error) return (size_t) -1;
	size_t offset = 0;
	for(i = 0; i < id; i++) {
/* we are sure there won't be errors here*/
		const char * blknm = _blockname(i, error);
		size_t bsize = _length(h, blknm, error) * _itemsize(blknm, error);
		if(bsize > 0) offset += (bsize + 8);
	}
	return offset;
}
  

static void _def_header(struct header_t * h) {
h->flag_entropy = 0;
h->flag_cloudy = 0;
h->flag_cool = 1;
h->flag_eos = 0;
h->flag_met = 1;
h->flag_Hemf = 0;
h->flag_sft = 1;
h->flag_sfr = 1;
h->flag_Helium = 0;
h->flag_gammaHI = 0;
h->flag_feedback = 1;
h->flag_Hmf = 0;
}


static void _get_constants(const struct header_t * h, struct constants_t * c) {
	int i;
	for(i = 0; i< 6; i++) {
		c->N[i] = h->N[i];
		c->Ntot[i] = h->Nparticle_total_low[i]
			+ (((size_t) h->Nparticle_total_high[i]) << 32);
		c->mass[i] = h->mass[i];
	}

c->OmegaM = h->OmegaM;
c->OmegaL = h->OmegaL;
c->time = h->time;
c->redshift = h->redshift;
c->Nfiles = h->Nfiles;
c->h = h->h;
c->boxsize = h->boxsize;

c->OmegaB = h->OmegaB;

}



static void _set_constants(struct header_t * h, const struct constants_t * c) {
	int i;
	for(i = 0; i< 6; i++) {
		h->N[i] = c->N[i];
		h->Nparticle_total_low[i] = c->Ntot[i];
		h->Nparticle_total_high[i] = c->Ntot[i] >> 32;
		h->mass[i] = c->mass[i];
	}

h->OmegaM = c->OmegaM;
h->OmegaL = c->OmegaL;
h->time = c->time;
h->redshift = c->redshift;
h->Nfiles = c->Nfiles;
h->h = c->h;
h->boxsize = c->boxsize;

h->OmegaB = c->OmegaB;

}



static void _read(struct header_t * h, const char * blk, void * buffer, int start, int length, FILE * fp, char ** error) {
	size_t off = _offset(h, blk, error);
	if(*error) return;
	size_t l = _length(h, blk, error);
	size_t b = _itemsize(blk, error);
	if(l == 0) return;
	size_t bsize = l * b;
	fseek(fp, off, SEEK_SET);
	int blksize = 0;
	fread(&blksize, sizeof(int), 1, fp);
	if(blksize != (int)bsize) {
		asprintf(error, "block start size of %s mismatch; file says %u, reader says %lu\n",
				blk, blksize, bsize);
		return;
	}
	fseek(fp, start * b, SEEK_CUR);
	if(!strcmp(blk, "header")) buffer = h;
	fread(buffer, b, length, fp);
	fseek(fp, (l - length - start) * b, SEEK_CUR);
	blksize = 0;
	fread(&blksize, sizeof(int), 1, fp);
	if(blksize != (int)bsize) {
		asprintf(error, "block end size of %s mismatch; file says %u, reader says %lu\n",
				blk, blksize, bsize);
		return;
	}
}
  


static void _write(const struct header_t * h, const char * blk, const void * buffer, int start, int length, FILE * fp, char ** error) {
	size_t off = _offset(h, blk, error);
	if(*error) return;
	size_t l = _length(h, blk, error);
	if(l == 0) return;
	size_t b = _itemsize(blk, error);
    size_t bsize = l * b;
	fseek(fp, off, SEEK_SET);
	int blksize = bsize;
	if(!strcmp(blk, "header")) buffer = h;
	fwrite(&blksize, sizeof(int), 1, fp);
	fseek(fp, start * b, SEEK_CUR);
	fwrite(buffer, b, length, fp);
	fseek(fp, (l - length - start) * b, SEEK_CUR);
	blksize = bsize;
	fwrite(&blksize, sizeof(int), 1, fp);
}
  


