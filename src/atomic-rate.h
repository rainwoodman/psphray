
extern int AR_LOG_T;
extern int AR_HI_CI;
extern int AR_HII_RC_A;
extern int AR_HII_RCC_A;

double ar_get(int id, double logT);
double ar_verner(double Ry);

void ar_init(char * filename);
void ar_dump(char * filename);

