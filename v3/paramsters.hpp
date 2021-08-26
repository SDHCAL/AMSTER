/*
 * parameters and global variables
 */

/* E edge type, V vertex type */
typedef double E;
typedef int    V;

/* parameters for Edmonds, Kruskal, Prim */
V     eventnumber     = 0;
float leakedenergy    = 0.0;
float depositedenergy = 0.0;
 
/* parameters for the example */
V           hits             = 15;
V           max              = 20;
E           trackenergy      = 35;
E           deltaE           = 1;
std::string format           = "png";

/* variables for time measurements and means */
auto t0 = std::chrono::high_resolution_clock::now();
auto t1 = std::chrono::high_resolution_clock::now();
auto t2 = std::chrono::high_resolution_clock::now();
auto t3 = std::chrono::high_resolution_clock::now();
auto t4 = std::chrono::high_resolution_clock::now();
auto t5 = std::chrono::high_resolution_clock::now();
auto t6 = std::chrono::high_resolution_clock::now();
std::vector<E> t(6,0.0);
std::vector<E> m(12,0.0);

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

/* improved calibration for SDHCALSim using hit density */
/*
double alpha0 = 0.0775376;
double alpha1 = 0.0486906;
double alpha2 = 0.062525;
double alpha3 = 0.0430229;
double alpha4 = 0.0700889;
double alpha5 = 0.0493535;
double alpha6 = 0.0322972;
double alpha7 = 0.0860822;
double alpha8 = 0.054635;
double beta0  = 0.111031;
double beta1  = 0.090856;
double beta2  = 0.0744228;
double beta3  = 0.0656031;
double beta4  = 0.0744405;
double beta5  = 0.0153327;
double beta6  = 0.0209479;
double beta7  = 0.00296485;
double beta8  = 0.14694;
double gamma0 = 0.224636;
double gamma1 = 0.172162;
double gamma2 = 0.0784331;
double gamma3 = 0.0416039;
double gamma4 = 0.00168717;
double gamma5 = 1.73649e-10;
double gamma6 = 3.96967e-10;
double gamma7 = 1.36294e-10;
double gamma8 = 0.33539;
*/
