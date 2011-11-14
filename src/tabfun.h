
typedef struct {
	ARRAY_DEFINE_S(data, double *);
	ARRAY_DEFINE_S(headers, char *);
	int nrows;
	int ncols;
	double min; 
	double max; 
	double step;
	double step_inv;
} TabFun;

const double tabfun_get(const TabFun * tabfun, const int id, const double value);
int tabfun_col_id(const TabFun * tabfun, const char * col);
void tabfun_init(TabFun * tabfun, const char * filename);
void tabfun_dump(const TabFun * tabfun, const char * filename);
int tabfun_setup_col(TabFun * tabfun, char * col, double (*func)(double), double unit);

typedef struct {
	ARRAY_DEFINE_S(data, double *);
	ARRAY_DEFINE_S(headers, char *);
	int nrows1;
	int nrows2;
	int ncols;
	double min1; 
	double min2; 
	double max1; 
	double max2; 
	double step1;
	double step2;
	double step1_inv;
	double step2_inv;
	ARRAY_DEFINE_S(x, double);
	ARRAY_DEFINE_S(y, double);
} TabFun2;

int tabfun2_ensure_col(TabFun2 * tabfun, char * col, double (*func)(double , double));
double tabfun2_get(TabFun2 * tabfun, int id, double x, double y);

