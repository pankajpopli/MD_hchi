#include"global.h"




/////////////// functions for these files are defined in /////////////////
// 1) dchi_square.c
// 2) dchi_triang.c
// 3) x_sq.c
// 4) x_triang.c and others
///////////////////////////////////////////////////////////////////////////



////////////////// local function declearations ///////////////////////// add here
double d_hchi_warning(int moved_particle,double *jstChi);
double d_hchi_triang(int moved_particle,double *jstChi);
double d_hchi_sq(int moved_particle,double *jstChi);
double d_hchi_readin(int moved_particle,double *jstChi);
double d_hchi_honeycomb(int moved_particle,double *jstChi);
double d_hchi_kagome(int moved_particle,double *jstChi);
double d_hchi_ring(int moved_particle,double *jstChi);
double d_hchi_rectangle(int moved_particle,double *jstChi);
//
double X_warning();
double X_triang();
double X_sq();
double X_readin();
double X_honeycomb();
double X_kagome();
double X_ring();
double X_rectangle();
//

//
void update_chi_warning(int mvd_part);
void update_chi_triang(int mvd_part);
void update_chi_sq(int mvd_part);
void update_chi_readin(int mvd_part);
void update_chi_honeycomb(int mvd_part);
void update_chi_kagome(int mvd_part);
void update_chi_ring(int mvd_part);
void update_chi_rectangle(int mvd_part);
//
void init_Y_mat_warning(int mvd_part);
void init_Y_mat_triang(int mvd_part);	 // check if the function definition is correct according to dchi_traing.c
void init_Y_mat_sq();							//  check if the function definition is correct according to dchi_sq.c
void init_Y_mat_readin();
void init_Y_mat_honeycomb();
void init_Y_mat_kagome();
void init_Y_mat_ring();
void init_Y_mat_rectangle();
//

/*
void update_X_warning(int mvd_part);
void update_X_triang(int mvd_part);
void update_X_sq(int mvd_part);
void update_X_readin(int mvd_part);
void update_X_honeycomb(int mvd_part);


void update_chi_X_full_warning();
void update_chi_X_full_triang();
void update_chi_X_full_sq();
void update_chi_X_full_readin();
void update_chi_X_full_honeycomb();
*/
////////////////////////////////////////////////////////////////////////////////





////////////////////////// pointer array function decleration ///////////////////////////////////////////////////////////////////////////// add here
double(*dchi_ptr_arr[])(int moved_particle, double *jstChi) = {d_hchi_warning,d_hchi_triang,d_hchi_sq,d_hchi_readin,d_hchi_honeycomb,d_hchi_kagome,d_hchi_ring,d_hchi_rectangle};
double(*Xfull_ptr_arr[])() = {X_warning,X_triang,X_sq,X_readin,X_honeycomb,X_kagome,X_ring,X_rectangle};
//void(*updtX_ptr_arr[])() = {update_X_warning,update_X_triang,update_X_sq,update_X_readin,update_X_honeycomb};
void(*updtChi_ptr_arr[])()= {update_chi_warning,update_chi_triang,update_chi_sq,update_chi_readin,update_chi_honeycomb,update_chi_kagome,update_chi_ring,update_chi_rectangle};
void(*initYmat_ptr_arr[])() = {init_Y_mat_warning,init_Y_mat_triang,init_Y_mat_sq,init_Y_mat_readin,init_Y_mat_honeycomb,init_Y_mat_kagome,init_Y_mat_ring,init_Y_mat_rectangle};
//void(*updtChiXfull_pt_arr[])() = {update_chi_X_full_warning,update_chi_X_full_triang,update_chi_X_full_sq,update_chi_X_full_readin,update_chi_X_full_honeycomb};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




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
double X()
{
	return (*Xfull_ptr_arr[Xfull_slctr])();
}

// [M3]
void update_chi(int mvd_part)
{
	(*updtChi_ptr_arr[updtChi_slctr])(mvd_part);
}
// {M4]
void init_Y_mat()
{
	(*initYmat_ptr_arr[initYmat_slctr])();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////





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
void init_Y_mat_warning(int mvd_part)
{
	printf(RED"\n X Y matrix not initialised, guess-> wrong lattice type in Global.c or wrong 'initialiser'\n"RESET);
	exit(0);
}
void init_Y_mat_readin()
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
double X_warning()
{
	printf(RED"\n  X not calculated, guess-> wrong lattice type in Global.c or wrong 'initialiser'\n"RESET);
	exit(0);
}
double X_readin()
{
	printf(RED"\n Readin is not supported yet\n"RESET);
	exit(0);
}