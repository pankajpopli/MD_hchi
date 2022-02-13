#include "global.h"


//--------------Type of potential and lattice-------------------------------//
//|put yes to print on terminal screen
char	print_config[50]	=  "yes";

//|possible use {no_pot, triang_ntwrk,triang_ntwrk_bending, sq_ntwrk, honeycomb_ntwrk,honeycomb_ntwrk_bending,kagome_ntwrk, r12, lj, GCM};
char	pot_type[50]	   =  "GCM";

//| possible use {triang, square,honeycomb,kagome,rectangle};
char  latt_type[50]     = "rectangle";

//----------------------Constants for the problem(system profile)----------------------------------//
// units scaled to 1---------------
double   K_B         = 1.0;
double 	sigma		  	= 1.0;
double   mass        = 1.0;
double   epsilon     = 1.0;
// parameters, star quantities-----
double   red_T             =  0.001;
double 	red_rho	  	      =  1.0;
double   damp              =  5.0;
double 	cutoff	  	      =  3.2;  //|in terms of sigma cutoff cell list. for potentials 
double 	cutoff_pot	      =  2.5;  //|potential cuttoff in terms of sigma cuttoff_pot > 3(cutoff)/2
double   hchi_0            =  0.00;
double   hchi_mean         =  -0.0015;
double   hchi_width        =  0.31;
int      h_step_intrvl     =  5000;
double   h_step            =  0.1; //|final hchi_0 = (nteps/h_step_intrvl)*h_ramp_step, plot change time = h_step_intrvl/print_freq
double   deltat            =  0.001;
double   k1                =  1.0;  //| not being used if potential is not ntwrk
double   k2                =  0.5;  //| not being used if potential is not ntwrk
//-----integration parameters-------
int 		nsteps	  	= 600000;      //|no of verlet iterations
int	 	printFreq  	= 100;         //|print freq to write data in files
int 	 	part_id	   = 10;			   //|	2*part_id + 2;
int 		totalN 	  	=  484;        //|   totalN = 2*part_id + 2;
#define lenarray (484)
int 	c_list[lenarray]	= {0};		//|initializing array to zero
Particle p[lenarray]	  	= {0.0};		//|initializing struct to zero
map_p2c	map[lenarray]	= {0}; 	   //|initializing struct to zero
double Ymat[4]={0.0,0.0,0.0,0.0};
//-----------------------------------------------------------------------------------------------//

//-----------------cell list and other constant parmaters---------------------------------------//
int 		*clen,*cbgn,ncell,ncell_x,ncell_y,c_updt;
double	lcell_x,lcell_y;
double 	rho,a,lbox_x,lbox_y;

//--------for pointer array functions---------------------------------------------------------//
int pot_slctr,clist_ntwrk_pot_slctr,latt_slctr,h_chi_force_slctr;
int n_list_slctr,updtChi_slctr,dchi_slctr,initYmat_slctr,Xfull_slctr;
char glatt[50],gpot[50];
//-------------------for LD dyanmics---------------------------------------------------------//
double stdDev;
//-------------------------------------------------------------------------------------------//

unsigned int thread_num=0;
































//--------------------------------------------------------------------------------------------------------------------------------------//

