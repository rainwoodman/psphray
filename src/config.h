#include <libconfig/libconfig.h>
extern config_t CFG[];
void cfg_init(char * filename);
void cfg_dump(char * filename);
double units_parse(char * units);

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

typedef struct {
	char * snapshot;
	char * source;
	char * format;
	char * redshift;
} Epoch;
