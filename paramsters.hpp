/*
 * parameters and global variables
 */

/* E edge type, V vertex type */
typedef double E;
typedef int    V;

/* parameters of the data for ./example */
int         numberofhits     = 15;
int         max              = 20;
E           trackenergy      = 35;
E           deltaE           = 1;
std::string format           = "png";

/* energy calibration for SDHCALSim */ 
double alpha0 = 0.0540243;
double alpha1 = -1.52143e-05;
double alpha2 = -1.58951e-09;
double beta0  = 0.0322513;
double beta1  = 0.000108338;
double beta2  = -6.84856e-08;
double gamma0 = 0.0986593;
double gamma1 = 0.000222751;
double gamma2 = 8.51319e-13;

/* energy calibration for SDHCAL prototype */
/*
double alpha0 = 0.022004;
double alpha1 = 3.13715e-05;
double alpha2 = -1.87027e-08;
double beta0  = 0.0953221;
double beta1  = 2.40243e-06;
double beta2  = -1.1425e-08;
double gamma0 = 0.182883;
double gamma1 = 2.00882e-05;
double gamma2 = 4.48185e-08;
*/

/* variables for time measurements */
auto   t0 = std::chrono::high_resolution_clock::now();
auto   t1 = std::chrono::high_resolution_clock::now();
auto   t2 = std::chrono::high_resolution_clock::now();
auto   t3 = std::chrono::high_resolution_clock::now();
auto   t4 = std::chrono::high_resolution_clock::now();
auto   t5 = std::chrono::high_resolution_clock::now();
auto   t6 = std::chrono::high_resolution_clock::now();
double t[6];
