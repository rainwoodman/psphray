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
} Epoch;

extern config_t CFG[];
extern Epoch * EPOCHS;
extern int N_EPOCHS;

void cfg_init(char * filename);
void cfg_dump(char * filename);
double config_setting_parse_units(config_setting_t * e);

double units_parse(char * units);
double units_format(double value, char * units);
double z2t(double z);
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

extern double C_HMF;
extern double C_OMEGA_L;
extern double C_OMEGA_M;
extern double C_H;
extern double C_BOLTZMANN;


