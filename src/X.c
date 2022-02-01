#include"global.h"

////////////////// local function declearations ///////////////////////// add here
double X_warning();
double X_triang();
double X_sq();
double X_readin();
double X_honeycomb();
double X_kagome();
double X_rectangle();

////////////////////////// pointer array function decleration ///////////////////////////////////////////////////////////////////////////// add here
double(*Xfull_ptr_arr[])() = {X_warning,X_triang,X_sq,X_readin,X_honeycomb,X_kagome,X_rectangle};



////////////////////////////////////// functions to be called from main.c ///////////////////////////////
// [M1]
double X()
{
	return (*Xfull_ptr_arr[Xfull_slctr])();
}

//[M1R*]
// All these functions are defined in seperate respective files
//////////////////////////////////////////////// warnings and readin functions definitions /////////////////////////////////
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