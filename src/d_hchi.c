#include"global.h"

////////////////// local function declearations ///////////////////////// add here
double d_hchi_warning(int moved_particle,double *jstChi);
double d_hchi_triang(int moved_particle,double *jstChi);
double d_hchi_sq(int moved_particle,double *jstChi);
double d_hchi_readin(int moved_particle,double *jstChi);
double d_hchi_honeycomb(int moved_particle,double *jstChi);
double d_hchi_kagome(int moved_particle,double *jstChi);
double d_hchi_rectangle(int moved_particle,double *jstChi);
//
void update_chi_warning(int mvd_part);
void update_chi_triang(int mvd_part);
void update_chi_sq(int mvd_part);
void update_chi_readin(int mvd_part);
void update_chi_honeycomb(int mvd_part);
void update_chi_kagome(int mvd_part);
void update_chi_rectangle(int mvd_part);
//

////////////////////////// pointer array function decleration ///////////////////////////////////////////////////////////////////////////// add here
void(*updtChi_ptr_arr[])()= {update_chi_warning,update_chi_triang,update_chi_sq,update_chi_readin,update_chi_honeycomb,update_chi_kagome,update_chi_rectangle};
double(*dchi_ptr_arr[])(int moved_particle, double *jstChi) = {d_hchi_warning,d_hchi_triang,d_hchi_sq,d_hchi_readin,d_hchi_honeycomb,d_hchi_kagome,d_hchi_rectangle};


////////////////////////////////////// functions to be called from main.c ///////////////////////////////
// [M1]
double d_hchi(int moved_particle, double *jstChi)
{
	double 	jstChi_dummy=0.0,return_val=0;
	return_val =  (*dchi_ptr_arr[dchi_slctr])(moved_particle, &jstChi_dummy);
	*jstChi = jstChi_dummy;
	return return_val;
}

// [M2]
void update_chi(int mvd_part)
{
	(*updtChi_ptr_arr[updtChi_slctr])(mvd_part);
}


//[M2R*] [M1R*]
// All these functions are defined in seperate respective files


//////////////////////////////////////////////// warnings and readin functions definitions /////////////////////////////////
void update_chi_warning(int mvd_part)
{
	printf(RED"\n chi matrix not updated, guess-> wrong lattice type in Global.c or wrong 'initialiser'\n"RESET);
	exit(0);
}
void update_chi_readin(int mvd_part)
{
	printf(RED"\n Readin is not supported yet\n"RESET);
	exit(0);
}
double d_hchi_warning(int moved_particle,double *jstChi)
{
	printf(RED"\n dchi not calculated, guess-> wrong lattice type in Global.c or wrong 'initialiser'\n"RESET);
	exit(0);
}
double d_hchi_readin(int moved_particle,double *jstChi)
{
	printf(RED"\n Readin is not supported yet\n"RESET);
	exit(0);
}