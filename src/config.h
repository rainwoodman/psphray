#include <libconfig/libconfig.h>
#include <gsl/gsl_rng.h>

typedef struct {
	double redshift;
	double age;
	double duration;
	size_t ngas;
	size_t nticks;
	size_t nray;
	const char * source;
	const char * snapshot;
	const char * format;
	struct {
		size_t * steps;
		size_t nsteps;
		size_t nfiles;
		const char * filename;
	} output;
} Epoch;

extern config_t CFG[];
extern int CFG_WRITE_INIT;
extern int CFG_ADIABATIC;
extern int CFG_ISOTHERMAL;
extern int CFG_DISABLE_2ND_GEN_PHOTONS;
extern int CFG_COMOVING;

extern gsl_rng * RNG;
extern Epoch * EPOCHS;
extern int N_EPOCHS;


void cfg_init(char * filename);
void cfg_dump(char * filename);
void cfg_dump_stream(char * FILE);

double config_setting_parse_units(config_setting_t * e);
config_setting_t * config_setting_ensure_member(config_setting_t * e, char * member, int type);
config_setting_t * config_ensure(config_t * config, char * path, int type);

double units_parse(char * units);
double units_format(double value, char * units);

extern double U_CM;
extern double U_GRAM;
extern double U_SEC;
extern double U_M;
extern double U_KM;
extern double U_KPC;
extern double U_YR;
extern double U_MYR;
extern double U_KG;
extern double U_MSUN;
extern double U_MPROTON;

extern double C_HMF;
extern double C_OMEGA_L;
extern double C_OMEGA_M;
extern double C_OMEGA_B;
extern double C_H;
extern double C_BOLTZMANN;
extern double U_MPROTON_OVER_C_BOLTZMANN;
extern double C_BOLTZMANN_OVER_U_MPROTON;

double z2t(double z);
double t2z(double z);
static inline float ieye2T(const float ie, const float ye) {
	return U_MPROTON_OVER_C_BOLTZMANN * ie / ( (1.0 - C_HMF) * 0.25 + 
			C_HMF + C_HMF * ye) * 2.0 / 3.0 ;

}
static inline float Tye2ie(const float T, const float ye) {
	return T * C_BOLTZMANN_OVER_U_MPROTON * ( (1.0 - C_HMF) * 0.25 + 
			C_HMF + C_HMF * ye) / 2.0 * 3.0 ;
}

extern int AR_LOG_T;
extern int AR_HI_CI;
extern int AR_HII_RC_A;
extern int AR_HII_RCC_A;
extern int AR_HII_RC_B;
extern int AR_HII_RCC_B;

double ar_get(int id, double logT);
double ar_verner(double Ry);

const int spec_get(const char * name);
double spec_gen_freq(const int id);

