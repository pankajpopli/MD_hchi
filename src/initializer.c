#include"global.h"

void initializer()
{

/*:::::::::::::::::::::::::::::::::::::::::::initialize to zero ::::::::::::::::::::::::::::::::::::::::::::::::::*/
   int ii=0,jj=0;
   for(ii=0;ii<totalN;ii++){
      p[ii].chi=0;
      p[ii].R[0]=p[ii].R[1]=0;   //these cords are used to calculate referenc for X
      p[ii].ri[0]=p[ii].ri[1]=0;
      p[ii].vel[0]=p[ii].vel[1]=0;
      p[ii].forceChi[0]=p[ii].forceChi[1]=0;
      for( jj=0;jj<10;jj++){p[ii].nn_static[jj]=0;}
      p[ii].pid=0;
      p[ii].cid=0;
		p[ii].h_xy=0;
   }
/*:::::::::::::::::::::::::::::::::::::::::::initialize square lattice const ::::::::::::::::::::::::::::::::::::::::::::::*/
      sprintf(glatt,"latt_%s",latt_type); //add lattices here
		if(strcmp(glatt,"latt_triang")==0)
		{
			//------------------constants for traing ----------------------------------------------------------
				printf(GRN"\n constans have been chosen for traingular latice in Global.c\n"RESET);
				rho		=	sqrt(((sqrt(3.0))*red_rho)/2.0); 
				a			=	sigma/rho;
			//----------------------------pointers for pointer array functions  ----------------------------
				latt_slctr				= 1;
				dchi_slctr 				= 1;
				Xfull_slctr 			= 1;
				n_list_slctr			= 1;
				updtChi_slctr 		   = 1;
				initYmat_slctr		   = 1;
				h_chi_force_slctr		= 1;
			//-------------------------------------------------------------------------------------------------	
		}else if(strcmp(glatt,"latt_square")==0)
		{
			//------------------constants for square ----------------------------------------------------------
				printf(GRN"\n constans have been chosen for square latice in Global.c\n"RESET);
				rho		=	sqrt(red_rho); 
				a			=	sigma/rho;
			//----------------------------pointers for pointer array functions  ----------------------------
				latt_slctr 				= 2;
				dchi_slctr 				= 2;
				Xfull_slctr 			= 2;
				n_list_slctr 			= 2;
				updtChi_slctr 		   = 2;
				initYmat_slctr 		= 2;
				h_chi_force_slctr		= 2;
			//-------------------------------------------------------------------------------------------------
		}else if(strcmp(glatt,"latt_readin")==0)
		{
			//------------------constants for readin ----------------------------------------------------------
				// write code here //
			//----------------------------pointers for pointer array functions  ----------------------------
				latt_slctr 				= 3;
				dchi_slctr 				= 3;
				Xfull_slctr 			= 3;
				n_list_slctr 			= 3;
				updtChi_slctr 	   	= 3;
				initYmat_slctr 		= 3;
				h_chi_force_slctr		= 0; //not implemented yet
			//-------------------------------------------------------------------------------------------------
		}else if(strcmp(glatt,"latt_honeycomb")==0)
		{
			//------------------constants for honeycomb ----------------------------------------------------------
				printf(GRN"\n constans have been chosen for honeycomb latice in Global.c\n"RESET);
				rho		=	sqrt(((3.0*sqrt(3.0)*red_rho)/(4.0))); 
				a			=	sigma/rho;
			//----------------------------pointers for pointer array functions  ----------------------------
				latt_slctr 				= 4;
				dchi_slctr 				= 4;
				Xfull_slctr 			= 4;
				n_list_slctr 			= 4;
				updtChi_slctr 		   = 4;
				initYmat_slctr 		= 4;
				h_chi_force_slctr		= 0; //not implemented yet
			//-------------------------------------------------------------------------------------------------
		}else if(strcmp(glatt,"latt_kagome")==0)
      {
         //------------------constants for honeycomb ----------------------------------------------------------
				printf(GRN"\n constans have been chosen for Kagome latice in Global.c\n"RESET);
				rho		=	sqrt((red_rho/(sqrt(12.0))));
				a			=	sigma/rho;
			//----------------------------pointers for pointer array functions  ----------------------------
				latt_slctr 				= 5;
				dchi_slctr 				= 5;
				Xfull_slctr 			= 5;
				n_list_slctr 			= 5;
				updtChi_slctr 		   = 5;
				initYmat_slctr 		= 5;
				h_chi_force_slctr		= 0; //not implemented yet
			//------------------------------------------------------------------------------------------------- 
      }else if(strcmp(glatt,"latt_rectangle")==0){
			//------------------constants for square ----------------------------------------------------------
				printf(GRN"\n constans have been chosen for rectangular latice in Global.c\n"RESET);
				rho		=	sqrt(((sqrt(3.0))*red_rho)/2.0); 
				a			=	sigma/rho;
			//----------------------------pointers for pointer array functions  ----------------------------
				latt_slctr 				= 6;
				dchi_slctr 				= 6;
				Xfull_slctr 			= 6;
				n_list_slctr 			= 6;
				updtChi_slctr 		   = 6;
				initYmat_slctr 		= 6;
				h_chi_force_slctr		= 3;
			//-------------------------------------------------------------------------------------------------
		}else
		{
            latt_slctr 				= 0;
            dchi_slctr 				= 0;
            Xfull_slctr 			= 0;
            n_list_slctr 			= 0;
            updtChi_slctr 		   = 0;
            initYmat_slctr 		= 0;
				h_chi_force_slctr		= 0;
		}
/*::::::::::::::::::::::::::::::::::::::::::: Select interatomic potential :::::::::::::::::::::::::::::::::::::::::::::::*/
      sprintf(gpot,"pot_%s",pot_type);
		if(strcmp(gpot,"pot_no_pot")==0){
         pot_slctr               = 1;
         clist_ntwrk_pot_slctr 	= 1;
      }
		else if(strcmp(gpot,"pot_triang_ntwrk")==0){
			pot_slctr 					= 2;
			clist_ntwrk_pot_slctr 	= 1;
		}else if(strcmp(gpot,"pot_triang_ntwrk_bend")==0){
			pot_slctr 					= 3;
			clist_ntwrk_pot_slctr 	= 1;
		}else if(strcmp(gpot,"pot_sq_ntwrk")==0){
			pot_slctr 					= 4;
			clist_ntwrk_pot_slctr 	= 1; //clist_ntwrk_pot_slctr = 1 if it is ntwrk otherwise 2
		}else if(strcmp(gpot,"pot_honeycomb_ntwrk")==0){
			pot_slctr 					= 5;
			clist_ntwrk_pot_slctr 	= 1;
		}else if(strcmp(gpot,"pot_honeycomb_ntwrk_bend")==0){
			pot_slctr 					= 6;
			clist_ntwrk_pot_slctr 	= 1;
		}
		else if(strcmp(gpot,"pot_kagome_ntwrk")==0){
         pot_slctr 					= 7;
			clist_ntwrk_pot_slctr 	= 1;
      }else if(strcmp(gpot,"pot_r12")==0){
			pot_slctr					= 8;
			clist_ntwrk_pot_slctr 	= 2;
		}else if(strcmp(gpot,"pot_lj")==0){
			pot_slctr 					= 9;
			clist_ntwrk_pot_slctr 	= 2;
		}
		else if(strcmp(gpot,"pot_GCM")==0){
			pot_slctr 					= 10;
			clist_ntwrk_pot_slctr 	= 2;
		}
      else{
			pot_slctr					= 0;
			clist_ntwrk_pot_slctr 	= 0;
		}
/*:::::::::::::::::::::::::::::::::::::::::::save parameter to a file ::::::::::::::::::::::::::::::::::::::::::::::::::*/
FILE *f;
f = fopen("system_profile.dat","w");
fprintf(f,"-------------------System profile below------------------\n");
fprintf(f,"K_B = %lf\nepsilon = %lf\nmass = %lf\ncell_list cutoff = %lf\ndelta_t = %lf\ndamp = %lf\n",K_B,epsilon,mass,cutoff,deltat,damp);
fprintf(f,"totalN = %d\nsigma = %lf\nlatt const(a) = %lf\nred_rho = %lf\n",totalN,sigma,a,red_rho);
fprintf(f,"hchi_width = %lf\nhchi_mean = %lf\nhchi_0 = %lf\n",hchi_width,hchi_mean);
fprintf(f,"cell list cutoff = %lf\nTemp = %lf\n",cutoff,red_T);
fprintf(f,"total MD steps = %d\nPrint freq = %d\ntotal configs generated = %d\n",nsteps,printFreq,nsteps/printFreq);
// fprintf(f,"k1 = %lf\nk2 = %lf\n",k1,k2); // not used
fclose(f);
}
