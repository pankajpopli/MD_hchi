//------------------force polar form -------------------//
#include "global.h"


/////////// local function for pot_array_pointer ///////////////////////////////////////////////////////////////////////////
//|Force_x and Force_y are only used while calculating spring interactions in sq_network potential
//|These variable are not used in other potential and thus assigned the same value which is the radial derivative
//|this radial derivative is then multiplied by cos and sin in cell_list to calculate cartesian forces.
void force_GCM(double dis, int mvdpart, double *Force_x, double *Force_y);
void force_lj(double dis, int mvdpart, double *Force_x, double *Force_y);
void force_r12(double dis,int mvdpart, double *Force_x, double *Force_y);
void force_sq_ntwrk(double dis, int mvdpart, double *Force_x, double *Force_y);
void force_no_pot(double dis, int mvdpart, double *Force_x, double *Force_y);
void force_warning(double dis, int mvdpart, double *Force_x, double *Force_y);
void force_triang_ntwrk(double dis, int mvdpart, double *Force_x, double *Force_y);
void force_triang_ntwrk_plus_bending(double dis, int mvdpart, double *Force_x, double *Force_y);
void force_honeycomb_ntwrk(double dis, int mvdpart, double *Force_x, double *Force_y);
void force_honeycomb_ntwrk_plus_bending(double dis, int mvdpart, double *Force_x, double *Force_y);
void force_kagome_ntwrk(double dis, int mvdpart, double *Force_x, double *Force_y);

////////////////////////////////////////// pointer array definition ///////////////////////////////////////////////////////
void (*force_ptr_arr[])(double dis, int mvdpart, double *Force_x, double *Force_y)  = {force_warning,force_no_pot,force_triang_ntwrk,force_triang_ntwrk_plus_bending, force_sq_ntwrk, force_honeycomb_ntwrk,force_honeycomb_ntwrk_plus_bending,force_kagome_ntwrk, force_r12,force_lj,force_GCM};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//#######################################################################################################################
//####################### function to be called from main ###############################################################
//#######################################################################################################################
//[M1] functions to be called from main @@@@
void force_intr_atmc(double distance, int idx, double *fi_x, double *fi_y)
{
	double force_x=0.0,force_y=0.0;
   (*force_ptr_arr[pot_slctr])(distance,idx, &force_x, &force_y );
   *fi_x = force_x; *fi_y = force_y;
}
//#######################################################################################################################




//#######################################################################################################################
//################################### define forces below will be used by M1 ############################################
//################################### through pointer array functions ###################################################
//#######################################################################################################################

void force_GCM(double dis, int mvdpart, double *Force_x, double *Force_y)
{
	// first Dimensionaless polar derivative of pot. below are cartsn component of forces
   double dVdr=0.0;
   dVdr = ((2.0*dis*exp(-(pow((dis/sigma),2))))/sigma)*(epsilon/sigma);   // supplyning back -( derivative of potential)
   *Force_x = dVdr;
   *Force_y = dVdr;
}
void force_lj(double dis, int mvdpart, double *Force_x, double *Force_y)
{
	// first Dimensionaless polar derivative of pot. below are cartsn component of forces
   double dVdr=0.0;
   dVdr = (4.0*( (12.0*pow((sigma/dis),13))-(6.0*pow((sigma/dis),7))))*(epsilon/sigma);   // supplyning back -( derivative of potential)
   *Force_x = dVdr;
   *Force_y = dVdr;
}
void force_r12(double dis,int mvdpart, double *Force_x, double *Force_y)
{
   double dVdr=0.0;
   dVdr = (12.0*pow((sigma/dis),13));  // supplyning back -( derivative of potential)
   *Force_x = dVdr;
   *Force_y = dVdr;
}
void force_no_pot(double dis,int mvdpart, double *Force_x, double *Force_y)
{
   double dummy;
   dummy=2;
}


void force_sq_ntwrk(double dis, int mvdpart, double *Force_x, double *Force_y)
{
	int bondno=0,Dbondno=0,dirctn=0,i=0;
	double sumx=0.0,sumy=0.0;
   double Rij_x=0.0,Rij_y=0.0,rij_x=0.0,rij_y=0.0,Rij=0.0,rij=0.0;
   
      for(i=0;i<totalN;i++){  
         sumx=0.0;sumy=0.0;
         Rij_x=0.0,Rij_y=0.0,rij_x=0.0,rij_y=0.0,Rij=0.0,rij=0.0;
         
         for(dirctn = 1;dirctn<5;dirctn++){
            
            bondno =	p[i].nn_static[dirctn];
            rij_x = pbcx(p[i].ri[0] - p[bondno].ri[0]);
            rij_y = pbcy(p[i].ri[1] - p[bondno].ri[1]);
            Rij = a;
            rij = PBCr( rij_x, rij_y);
            sumx -= k1*(1.0- (Rij)/rij)*rij_x; 
            sumy -= k1*(1.0- (Rij)/rij)*rij_y;
         }
         for(dirctn = 5;dirctn<9;dirctn++){
            
            Dbondno =	p[i].nn_static[dirctn];
            rij_x = pbcx(p[i].ri[0] - p[Dbondno].ri[0]);
            rij_y = pbcy(p[i].ri[1] - p[Dbondno].ri[1]);
            Rij = sqrt(2.0)*a;
            rij = PBCr( rij_x, rij_y);
            sumx -= k2*(1.0- (Rij)/rij)*rij_x; 
            sumy -= k2*(1.0- (Rij)/rij)*rij_y;
         }
         p[i].forceC[0] = sumx; p[i].forceC[1] = sumy;
      }
}

void force_warning(double dis, int mvdpart, double *Force_x, double *Force_y)
{
	printf(RED"\n Error ! Wrong choice of potential in Global.c\n"RESET);
	exit(0);
}

void force_triang_ntwrk(double dis, int mvdpart, double *Force_x, double *Force_y)
{
	int bondno=0,dirctn=0,i=0;
	double sumx=0.0,sumy=0.0;
   double Rij_x=0.0,Rij_y=0.0,rij_x=0.0,rij_y=0.0,Rij=0.0,rij=0.0;
	
   for(i=0;i<totalN;i++){ 
      
      sumx=0.0;sumy=0.0;
      Rij_x=0.0,Rij_y=0.0,rij_x=0.0,rij_y=0.0,Rij=0.0,rij=0.0;
      
      for(dirctn = 1;dirctn<7;dirctn++){
         bondno =	p[i].nn_static[dirctn];
         rij_x = pbcx(p[i].ri[0] - p[bondno].ri[0]);
         rij_y = pbcy(p[i].ri[1] - p[bondno].ri[1]);
         Rij = a;
         rij = PBCr( rij_x, rij_y);
         sumx -= k1*(1.0- (Rij)/rij)*rij_x; 
         sumy -= k1*(1.0- (Rij)/rij)*rij_y;
      }
      p[i].forceC[0] = sumx; p[i].forceC[1] = sumy;
   }
}
void force_triang_ntwrk_plus_bending(double dis, int mvdpart, double *Force_x, double *Force_y){
   printf(RED" Triangle bending force not implemented yet. Exiting the code now."RESET);exit(0);
}
void force_honeycomb_ntwrk(double dis, int mvdpart, double *Force_x, double *Force_y)
{	
	int bondno=0,Dbondno=0,dirctn=0,i=0;
	double sumx=0.0,sumy=0.0;
   double Rij_x=0.0,Rij_y=0.0,rij_x=0.0,rij_y=0.0,Rij=0.0,rij=0.0;
   
      for(i=0;i<totalN;i++){
         
         sumx=0.0;sumy=0.0;
         Rij_x=0.0,Rij_y=0.0,rij_x=0.0,rij_y=0.0,Rij=0.0,rij=0.0;
         
         for(dirctn = 1;dirctn<4;dirctn++){
            
            bondno =	p[i].nn_static[dirctn];
            rij_x = pbcx(p[i].ri[0] - p[bondno].ri[0]);
            rij_y = pbcy(p[i].ri[1] - p[bondno].ri[1]);
            Rij = a;
            rij = PBCr( rij_x, rij_y);
            sumx -= k1*(1.0- (Rij)/rij)*rij_x; 
            sumy -= k1*(1.0- (Rij)/rij)*rij_y;
         }
         for(dirctn = 4;dirctn<10;dirctn++){

            Dbondno =	p[i].nn_static[dirctn];
            rij_x = pbcx(p[i].ri[0] - p[Dbondno].ri[0]);
            rij_y = pbcy(p[i].ri[1] - p[Dbondno].ri[1]);
            Rij = sqrt(3.0)*a;
            rij = PBCr( rij_x, rij_y);
            sumx -= k2*(1.0- (Rij)/rij)*rij_x; 
            sumy -= k2*(1.0- (Rij)/rij)*rij_y;
            
         }
         p[i].forceC[0] = sumx; p[i].forceC[1] = sumy;
      }
}
void force_honeycomb_ntwrk_plus_bending(double dis, int mvdpart, double *Force_x, double *Force_y){
   printf(RED" Triangle bending force not implemented yet. Exiting the code now."RESET);exit(0);
}

void force_kagome_ntwrk(double dis, int mvdpart, double *Force_x, double *Force_y)
{	
	int bondno=0,Dbondno=0,dirctn=0,i=0;
	double sumx=0.0,sumy=0.0;
   double Rij_x=0.0,Rij_y=0.0,rij_x=0.0,rij_y=0.0,Rij=0.0,rij=0.0;
   
   for(i=0;i<totalN;i++){
      
      sumx=0.0;sumy=0.0;
      Rij_x=0.0,Rij_y=0.0,rij_x=0.0,rij_y=0.0,Rij=0.0,rij=0.0;
      
      for(dirctn = 1;dirctn<5;dirctn++){
         
         bondno =	p[i].nn_static[dirctn];
         rij_x = pbcx(p[i].ri[0] - p[bondno].ri[0]);
         rij_y = pbcy(p[i].ri[1] - p[bondno].ri[1]);
         Rij = PBCr( (p[i].R[0]-p[bondno].R[0]), (p[i].R[1]-p[bondno].R[1]));
         rij = PBCr( rij_x, rij_y);
         sumx -= k1*(1.0- (Rij)/rij)*rij_x; 
         sumy -= k1*(1.0- (Rij)/rij)*rij_y;
      }
      for(dirctn = 5;dirctn<9;dirctn++){
         
         Dbondno =	p[i].nn_static[dirctn];
         rij_x = pbcx(p[i].ri[0] - p[Dbondno].ri[0]);
         rij_y = pbcy(p[i].ri[1] - p[Dbondno].ri[1]);
         Rij = PBCr( (p[i].R[0]-p[Dbondno].R[0]), (p[i].R[1]-p[Dbondno].R[1]));
         rij = PBCr( rij_x, rij_y);
         sumx -= k2*(1.0- (Rij)/rij)*rij_x; 
         sumy -= k2*(1.0- (Rij)/rij)*rij_y;
      }
      p[i].forceC[0] = sumx; p[i].forceC[1] = sumy;
   }
}
//#######################################################################################################################