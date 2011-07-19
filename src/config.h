#include <libconfig/libconfig.h>
typedef struct {
	double redshift;
	double age;
	double duration;
	size_t ngas;
	size_t nticks;
	const char * source;
	const char * snapshot;
	const char * format;
	const char * output;
	size_t output_nfiles;
	size_t output_nsteps;
} Epoch;

extern config_t CFG[];
extern int CFG_WRITE_INIT;

extern Epoch * EPOCHS;
extern int N_EPOCHS;


void cfg_init(char * filename);
void cfg_dump(char * filename);
double config_setting_parse_units(config_setting_t * e);
config_setting_t * config_setting_ensure_member(config_setting_t * e, char * member, int type);
config_setting_t * config_ensure(config_t * config, char * path, int type);

double units_parse(char * units);
double units_format(double value, char * units);
double z2t(double z);
double t2z(double z);
float ieye2T(const float ie, const float ye);

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


extern int AR_LOG_T;
extern int AR_HI_CI;
extern int AR_HII_RC_A;
extern int AR_HII_RCC_A;
extern int AR_HII_RC_B;
extern int AR_HII_RCC_B;

double ar_get(int id, double logT);
double ar_verner(double Ry);

