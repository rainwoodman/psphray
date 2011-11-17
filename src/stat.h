typedef struct {
	struct {
		size_t total;
		size_t subtotal;
	} src_ray_count;
	struct {
		size_t total;
		size_t subtotal;
	} rec_ray_count;
	struct {
		double total;
		double subtotal;
	} src_photon_count;
	struct {
		double total;
		double subtotal;
		double subtotalHII;
		double subtotalHeII;
		double subtotalHeIII;
	} rec_photon_count;

	struct {
		double HI;
		double HeI;
	} secondary_ionization;
	struct {
		double HI;
		double HeI;
		double HeII;
	} first_ionization;

/* total number of photon travel out of box/ remain in pretruncated ray*/
	double lost_photon_count_sum;
/* total number of recombination photons */
	double rec_photon_count_sum;

	size_t saturated_deposit_count;
	size_t total_deposit_count;
	size_t disordered_count;
	size_t gsl_error_count;
	size_t evolve_count;
	size_t tick_subtotal;
	size_t fast_recombination_count;

	FILE * parlogfile;
	FILE * hitlogfile;

	double spinlock_time;
	double deposit_time;
	double raytrace_time;
	double update_time;
	double emit_time;
	double merge_time;
	double total_time;
} stat_t;

extern stat_t stat;
void stat_restart();
void stat_subtotal();
void stat_stop();
