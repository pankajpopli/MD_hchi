#ifndef GLOBAL_H_
#define GLOBAL_H_
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include<string.h>
#include<stdbool.h>
#include<complex.h>
#include <gsl/gsl_randist.h>
#include<sys/stat.h>
#include <unistd.h>
#include<time.h>
#include<omp.h>
#ifndef PI
	#define PI (M_PI)
#endif

// color codes
#define RED   "\x1B[31m"      //|Errors
#define GRN   "\x1B[32m"      //|Messages
#define YEL   "\x1B[33m"      //|Warnings
#define BLU   "\x1B[34m"      //|
#define MAG   "\x1B[35m"      //|
#define CYN   "\x1B[36m"      //|
#define WHT   "\x1B[37m"      //|Information
#define RESET "\x1B[0m"

//**************************************************************************************************************//
//*******************************************structure used****************************************************//
//************************************************************************************************************//
// [1] structure for particles
typedef struct particle{
	double chi;
   double R[2];   //these cords are used to calculate referenc for X
   double ri[2];
   double ri0[2];
   double vel[2];
   // double vel_vec[2];
   // double velVck[2];
   // double Ymat[4];
   // double chi_temp;
   double h_xy;
   double forceChi[2];
   double forceC[2];
   int 	 nn_static[10];
   int	 pid;
   int	 cid;
}Particle;
extern  Particle p[];
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// [2] stores nbr cord. of centre particle in generate_lattice2.c traingular or square used only once-------//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//[3] used to find cid for any given pid (pid = index of array map)-------------------//
typedef struct Map_p2c{
	int pid;
	int cid;
}map_p2c;
extern map_p2c map[];
//#######################################################################################################################





//************************************************************************************************************************//
//****************************************** Some default parameters ****************************************************//
//**********************************************************************************************************************//
//[1] constants and array for monte carlo and lattice system
extern double  sigma;
extern double  red_rho;
extern double  cutoff;
extern double  cutoff_pot;
extern double  hchi_0;
extern double  hchi_mean;
extern double  hchi_width;
extern double  h_step;
extern double  k1;
extern double  k2;
extern double  epsilon;
extern double  K_B;
extern double  red_T;
extern int		nsteps;
extern int		printFreq;
extern int 		part_id;
extern int 		totalN; /// ((2*part_id)^2 + 2) 
extern int     h_step_intrvl;
extern char 	print_config[50];
extern char    pot_type[50];
extern double  rho,a,lbox_x,lbox_y,Ymat[4];
//[2] constants and array used in cell list
extern double 	lcell_x, lcell_y;
extern int		ncell_x, ncell_y, ncell, c_list[], *clen, *cbgn, c_updt;
//[3] constants for pointer arrays
extern int pot_slctr,clist_ntwrk_pot_slctr,latt_slctr,h_chi_force_slctr;
extern int n_list_slctr,updtChi_slctr,dchi_slctr,initYmat_slctr,Xfull_slctr;
extern char glatt[500],gpot[500],latt_type[50];
//[4] constants for integrator
extern double deltat,mass,damp,eta,stdDev;
//[5] constants for vchk model
// extern double Vck0,Vnoise,cutoffV;
//[6] constants for turbulent vel field
//#########################################################################################################################################################

extern unsigned int thread_num;



//******************************************************************************************************************************************************//
//****************************************************** functions Used in code ***********************************************************************//
//****************************************************************************************************************************************************//


//~~~~~~~~~~~~~~~intitailizers~~~~~~~~~~~~~~~~~~~~~~~~~~~
void initializer();
//[1] Generate lattice traing/square but index not ordered
void generate_latt();
//[2] Creates a nbr list which does not change with time.
void n_list_static();
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


//~~~~~~~~~~~~~~~~potentials & PBC ~~~~~~~~~~~~~~~~~~~~~
//---PBC functions--
///[3] Applies pbc on x-cordinate and output the new dx
double pbcx(double dx);
//[4] Applies pbc on y-cordinate and output the new dy
double pbcy(double dy);
//[5] Applies pbc on radial distance and output the new dr
double PBCr(double rx, double ry);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


//~~~~~~~~~~~~~~~~cell list functions~~~~~~~~~~~~~~~~
//[6] Creates cell_list
void create_c_list();
//[7] Updates the cell list full...works like[10]
void update_c_list_full();
double c_list_totalE(char *string);
void c_list_forceC();
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


//~~~~~~~~~~~~~~~~chi and Y matrix calculations~~~~~~~~~
//[8] Initalize Y matrix for the calc of Chi
void init_Y_mat();
//[9] Calculate total Chi for a given configuration
double X();
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



//~~~~~~~~~~~~~~~~~force calc~~~~~~~~~~~~~~~~~~~~~~~~~~~
//[11] calclte forces because of non-affine parameter
void h_chi_force();
double calc_Temp();
double potential(double distance, int mvdpartIndex);
void force_intr_atmc(double distance, int idx, double *fi_x, double *fi_y);
double h_field(double distance);
void spatial_random_h(const gsl_rng *ran3);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



//#############################################################################################################################


#endif /* GLOBAL_H_ */








