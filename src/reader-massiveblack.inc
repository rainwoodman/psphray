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
int flag_double;
int flag_ic_info;
int flag_lpt_scalingfactor;
int unused[12];
};

static int _blockid(const char * blk, char ** error) {
	if(!strcmp(blk, "pos")) return 1 ;
	if(!strcmp(blk, "vel")) return 2 ;
	if(!strcmp(blk, "id")) return 3 ;
	if(!strcmp(blk, "mass")) return 4 ;
	if(!strcmp(blk, "ie")) return 5 ;
	if(!strcmp(blk, "rho")) return 6 ;
	if(!strcmp(blk, "ye")) return 7 ;
	if(!strcmp(blk, "xHI")) return 8 ;
	if(!strcmp(blk, "sml")) return 9 ;
	if(!strcmp(blk, "sfr")) return 10 ;
	if(!strcmp(blk, "sft")) return 11 ;
	if(!strcmp(blk, "met")) return 12 ;
	if(!strcmp(blk, "bhmass")) return 13 ;
	if(!strcmp(blk, "bhmdot")) return 14 ;
	if(!strcmp(blk, "header")) return 0;
	asprintf(error, "block %s is unknown", blk);
	return -1;}

static const char * _blockname(const int blockid, char ** error) {
	static char * names[] = {"header", 
"pos","vel","id","mass","ie","rho","ye","xHI","sml","sfr","sft","met","bhmass","bhmdot"
};

	if(blockid < 0 || blockid >= sizeof(names) / sizeof(char*)) {
		asprintf(error, "blockid %d is unknown", blockid);
		return NULL;
	}
	return names[blockid];
     
}

static size_t _itemsize(const char * blk, char ** error) {
	if(!strcmp(blk, "pos")) return 12 ;
	if(!strcmp(blk, "vel")) return 12 ;
	if(!strcmp(blk, "id")) return 8 ;
	if(!strcmp(blk, "mass")) return 4 ;
	if(!strcmp(blk, "ie")) return 4 ;
	if(!strcmp(blk, "rho")) return 4 ;
	if(!strcmp(blk, "ye")) return 4 ;
	if(!strcmp(blk, "xHI")) return 4 ;
	if(!strcmp(blk, "sml")) return 4 ;
	if(!strcmp(blk, "sfr")) return 4 ;
	if(!strcmp(blk, "sft")) return 4 ;
	if(!strcmp(blk, "met")) return 4 ;
	if(!strcmp(blk, "bhmass")) return 4 ;
	if(!strcmp(blk, "bhmdot")) return 4 ;
	if(!strcmp(blk, "header")) return 256;
	asprintf(error, "block %s is unknown", blk);
	return (size_t)-1;}

static size_t _pstart(const struct header_t * h, const char * blk, int ptype, char ** error) {
	if(!strcmp(blk, "pos")) return (0 < ptype)?(size_t)(h->N[0]): 0+(1 < ptype)?(size_t)(h->N[1]): 0+(4 < ptype)?(size_t)(h->N[4]): 0+(5 < ptype)?(size_t)(h->N[5]): 0 ;
	if(!strcmp(blk, "vel")) return (0 < ptype)?(size_t)(h->N[0]): 0+(1 < ptype)?(size_t)(h->N[1]): 0+(4 < ptype)?(size_t)(h->N[4]): 0+(5 < ptype)?(size_t)(h->N[5]): 0 ;
	if(!strcmp(blk, "id")) return (0 < ptype)?(size_t)(h->N[0]): 0+(1 < ptype)?(size_t)(h->N[1]): 0+(4 < ptype)?(size_t)(h->N[4]): 0+(5 < ptype)?(size_t)(h->N[5]): 0 ;
	if(!strcmp(blk, "mass")) return (0 < ptype)?(size_t)(h->N[0]): 0+(4 < ptype)?(size_t)(h->N[4]): 0+(5 < ptype)?(size_t)(h->N[5]): 0 ;
	if(!strcmp(blk, "ie")) return (0 < ptype)?(size_t)(h->N[0]): 0 ;
	if(!strcmp(blk, "rho")) return (0 < ptype)?(size_t)(h->N[0]): 0 ;
	if(!strcmp(blk, "ye")) return (0 < ptype)?(size_t)(h->N[0]): 0 ;
	if(!strcmp(blk, "xHI")) return (0 < ptype)?(size_t)(h->N[0]): 0 ;
	if(!strcmp(blk, "sml")) return (0 < ptype)?(size_t)(h->N[0]): 0 ;
	if(!strcmp(blk, "sfr")) return (0 < ptype)?(size_t)(h->N[0]): 0 ;
	if(!strcmp(blk, "sft")) return (4 < ptype)?(size_t)(h->N[4]): 0 ;
	if(!strcmp(blk, "met")) return (0 < ptype)?(size_t)(h->N[0]): 0+(4 < ptype)?(size_t)(h->N[4]): 0 ;
	if(!strcmp(blk, "bhmass")) return (5 < ptype)?(size_t)(h->N[5]): 0 ;
	if(!strcmp(blk, "bhmdot")) return (5 < ptype)?(size_t)(h->N[5]): 0 ;
	asprintf(error, "block %s is unknown/unstable for _pstart", blk);
	return (size_t)-1;}

static size_t _length(const struct header_t * h, const char * blk, char ** error) {
	if(!strcmp(blk, "pos")) return (size_t)(h->N[0])+(size_t)(h->N[1])+(size_t)(h->N[4])+(size_t)(h->N[5]) ;
	if(!strcmp(blk, "vel")) return (size_t)(h->N[0])+(size_t)(h->N[1])+(size_t)(h->N[4])+(size_t)(h->N[5]) ;
	if(!strcmp(blk, "id")) return (size_t)(h->N[0])+(size_t)(h->N[1])+(size_t)(h->N[4])+(size_t)(h->N[5]) ;
	if(!strcmp(blk, "mass")) return (size_t)(h->N[0])+(size_t)(h->N[4])+(size_t)(h->N[5]) ;
	if(!strcmp(blk, "ie")) return (size_t)(h->N[0]) ;
	if(!strcmp(blk, "rho")) return (size_t)(h->N[0]) ;
	if(!strcmp(blk, "ye")) return (size_t)(h->N[0]) ;
	if(!strcmp(blk, "xHI")) return (size_t)(h->N[0]) ;
	if(!strcmp(blk, "sml")) return (size_t)(h->N[0]) ;
	if(!strcmp(blk, "sfr")) return (size_t)(h->N[0]) ;
	if(!strcmp(blk, "sft")) return (size_t)(h->N[4]) ;
	if(!strcmp(blk, "met")) return (size_t)(h->N[0])+(size_t)(h->N[4]) ;
	if(!strcmp(blk, "bhmass")) return (size_t)(h->N[5]) ;
	if(!strcmp(blk, "bhmdot")) return (size_t)(h->N[5]) ;
	if(!strcmp(blk, "header")) return 1;
	asprintf(error, "block %s is unknown", blk);
	return (size_t)-1;}


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
h->flag_sft = 1;
h->flag_sfr = 1;
h->flag_double = 0;
h->flag_cool = 1;
h->flag_feedback = 1;
h->flag_ic_info = 0;
h->flag_met = 1;
}


static void _get_constants(const struct header_t * h, struct constants_t * c) {
	int i;
	for(i = 0; i< 6; i++) {
		c->N[i] = h->N[i];
		c->Ntot[i] = h->Nparticle_total_low[i]
			+ (((size_t) h->Nparticle_total_high[i]) << 32);
		c->mass[i] = h->mass[i];
	}
	c->boxsize = h->boxsize;
	c->OmegaM = h->OmegaM;
	c->OmegaL = h->OmegaL;
	c->redshift = h->redshift;
	c->time = h->time;
	c->Nfiles = h->Nfiles;
	c->h = h->h;
}
static void _set_constants(struct header_t * h, const struct constants_t * c) {
	int i;
	for(i = 0; i< 6; i++) {
		h->N[i] = c->N[i];
		h->Nparticle_total_low[i] = c->Ntot[i];
		h->Nparticle_total_high[i] = c->Ntot[i] >> 32;
		h->mass[i] = c->mass[i];
	}
	h->boxsize = c->boxsize;
	h->OmegaM = c->OmegaM;
	h->OmegaL = c->OmegaL;
	h->redshift = c->redshift;
	h->time = c->time;
	h->Nfiles = c->Nfiles;
	h->h = c->h;
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
  


