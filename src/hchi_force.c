#include"global.h"


/////////// local function for pot_array_pointer ///////////////////////////////////////////////////////////////////////////
void h_chi_force_triangle();
void h_chi_force_square();
void h_chi_force_rectangle();
void h_chi_force_warning();

////////////////////////////////////////// pointer array definition ///////////////////////////////////////////////////////
void (*h_chi_force_ptr_arr[])()={h_chi_force_warning,h_chi_force_triangle,h_chi_force_square,h_chi_force_rectangle};


//#######################################################################################################################
//####################### function to be called from main ###############################################################
//#######################################################################################################################
//[M1] functions to be called from main @@@@
void h_chi_force()
{
	double force_x=0.0,force_y=0.0;
   (*h_chi_force_ptr_arr[h_chi_force_slctr])();
}



//#######################################################################################################################
// The below routine for hchi forces SHOULD work for square, triangle and square lattice. 
//since the analytical expression same, but do check if hchi forces here work for traingular lattice also.
// hchi forces for the ring MAYBE a more generalized expression of hchi forces.
void h_chi_force_triangle()
{
	
   double Deform[2][2]={{0.0,0.0},{0.0,0.0}};
   double sqRsum = (1.0/(3.0*a*a));
   double rij_x=0.0,rij_y=0.0,Rij_x=0.0,Rij_y=0.0,sumx=0.0,sumy=0.0;
   int i=0,n=0,m=0,nbr=7,nbr_of_i=0,nbr_of_n=0;
   
   //omp_set_dynamic(0);     // Explicitly disable dynamic teams
   //omp_set_num_threads(thread_num);
   //#pragma omp parallel for shared(p,sqRsum,hchi,nbr) private(sumx,sumy,nbr_of_i,Deform,rij_x,rij_y,Rij_x,Rij_y,n,m,nbr_of_n,i) 
   for(i=0;i<totalN;i++)
   {
      sumx=0.0;sumy=0.0;
      for(n=1;n<nbr;n++){
         
         nbr_of_i = p[i].nn_static[n]; //(N is nbr of I)
         Deform[0][0]=0.0; Deform[0][1]=0.0;Deform[1][0]=0.0;Deform[1][1]=0.0;
         for(m=1;m<nbr;m++){
            
            nbr_of_n = p[nbr_of_i].nn_static[m]; //  (M is nbr of nbr of I)
            rij_x = pbcx( p[nbr_of_n].ri[0] - p[nbr_of_i].ri[0] );
            rij_y = pbcy( p[nbr_of_n].ri[1] - p[nbr_of_i].ri[1] );
            Rij_x = pbcx( p[nbr_of_n].R[0]  - p[nbr_of_i].R[0]  );
            Rij_y = pbcy( p[nbr_of_n].R[1]  - p[nbr_of_i].R[1]  );
         
            Deform[0][0] += rij_x*Rij_x;
            Deform[0][1] += rij_x*Rij_y;
            Deform[1][0] += rij_y*Rij_x;
            Deform[1][1] += rij_y*Rij_y;
         }
      
         Deform[0][0] = (Deform[0][0]*sqRsum)  - 1.0; Deform[0][1] = (Deform[0][1]*sqRsum)  - 0.0;
         Deform[1][0] = (Deform[1][0]*sqRsum)  - 0.0; Deform[1][1] = (Deform[1][1]*sqRsum)  - 1.0;
   
         rij_x = pbcx( p[i].ri[0] - p[nbr_of_i].ri[0]  );
         rij_y = pbcy( p[i].ri[1] - p[nbr_of_i].ri[1]  );
         Rij_x = pbcx( p[i].R[0]  - p[nbr_of_i].R[0]   );
         Rij_y = pbcy( p[i].R[1]  - p[nbr_of_i].R[1]   );
         
         sumx += (4.0*rij_x) - 2.0*( ((Deform[0][0]+1.0)*Rij_x) + ((Deform[0][1])*Rij_y));
         sumy += (4.0*rij_y) - 2.0*( ((Deform[1][0])*Rij_x) + ((Deform[1][1]+1.0)*Rij_y));
      }
      // p[i].forceChi[0] = (hchi*sumx); p[i].forceChi[1] = (hchi*sumy);
      p[i].forceChi[0] = hchi_0*(p[i].h_xy*sumx); p[i].forceChi[1] = hchi_0*(p[i].h_xy*sumy);
   }
   

}
void h_chi_force_square()
{
	
   double Deform[2][2]={{0.0,0.0},{0.0,0.0}};
   double sqRsum = (1.0/(6.0*a*a));
   double rij_x=0.0,rij_y=0.0,Rij_x=0.0,Rij_y=0.0,sumx=0.0,sumy=0.0;
   int i=0,n=0,m=0,nbr=9,nbr_of_i=0,nbr_of_n=0;
   
   //omp_set_dynamic(0);     // Explicitly disable dynamic teams
   //omp_set_num_threads(thread_num);
   //#pragma omp parallel for shared(p,sqRsum,hchi,nbr) private(sumx,sumy,nbr_of_i,Deform,rij_x,rij_y,Rij_x,Rij_y,n,m,nbr_of_n,i) 
   for(i=0;i<totalN;i++)
   {
      sumx=0.0;sumy=0.0;
      for(n=1;n<nbr;n++){
         
         nbr_of_i = p[i].nn_static[n]; //(N is nbr of I)
         Deform[0][0]=0.0; Deform[0][1]=0.0;Deform[1][0]=0.0;Deform[1][1]=0.0;
         for(m=1;m<nbr;m++){
            
            nbr_of_n = p[nbr_of_i].nn_static[m]; //  (M is nbr of nbr of I)
            rij_x = pbcx( p[nbr_of_n].ri[0] - p[nbr_of_i].ri[0] );
            rij_y = pbcy( p[nbr_of_n].ri[1] - p[nbr_of_i].ri[1] );
            Rij_x = pbcx( p[nbr_of_n].R[0]  - p[nbr_of_i].R[0]  );
            Rij_y = pbcy( p[nbr_of_n].R[1]  - p[nbr_of_i].R[1]  );
         
            Deform[0][0] += rij_x*Rij_x;
            Deform[0][1] += rij_x*Rij_y;
            Deform[1][0] += rij_y*Rij_x;
            Deform[1][1] += rij_y*Rij_y;
         }
      
         Deform[0][0] = (Deform[0][0]*sqRsum)  - 1.0; Deform[0][1] = (Deform[0][1]*sqRsum)  - 0.0;
         Deform[1][0] = (Deform[1][0]*sqRsum)  - 0.0; Deform[1][1] = (Deform[1][1]*sqRsum)  - 1.0;
   
         rij_x = pbcx( p[i].ri[0] - p[nbr_of_i].ri[0]  );
         rij_y = pbcy( p[i].ri[1] - p[nbr_of_i].ri[1]  );
         Rij_x = pbcx( p[i].R[0]  - p[nbr_of_i].R[0]   );
         Rij_y = pbcy( p[i].R[1]  - p[nbr_of_i].R[1]   );
         
         sumx += (4.0*rij_x) - 2.0*( ((Deform[0][0]+1.0)*Rij_x) + ((Deform[0][1])*Rij_y));
         sumy += (4.0*rij_y) - 2.0*( ((Deform[1][0])*Rij_x) + ((Deform[1][1]+1.0)*Rij_y));
      }
      // p[i].forceChi[0] = (hchi*sumx); p[i].forceChi[1] = (hchi*sumy);
      p[i].forceChi[0] = hchi_0*(p[i].h_xy*sumx); p[i].forceChi[1] = hchi_0*(p[i].h_xy*sumy);
   }
   

}
void h_chi_force_rectangle()
{
	// printf("\nhchi value = %lf\n",hchi);
   double Deform[2][2]={{0.0,0.0},{0.0,0.0}};
   double sqRsum1 = (1.0/(6.0*a*a));
   double sqRsum2 = (1.0/(6.0*a*a* (3.0/4.0) ));
   double rij_x=0.0,rij_y=0.0,Rij_x=0.0,Rij_y=0.0,sumx=0.0,sumy=0.0;
   int i=0,n=0,m=0,nbr=9,nbr_of_i=0,nbr_of_n=0;
   
   //omp_set_dynamic(0);     // Explicitly disable dynamic teams
   //omp_set_num_threads(thread_num);
   //#pragma omp parallel for shared(p,sqRsum,hchi,nbr) private(sumx,sumy,nbr_of_i,Deform,rij_x,rij_y,Rij_x,Rij_y,n,m,nbr_of_n,i) 
   for(i=0;i<totalN;i++)
   {
      sumx=0.0;sumy=0.0;
      for(n=1;n<nbr;n++){
         
         nbr_of_i = p[i].nn_static[n]; //(N is nbr of I)
         Deform[0][0]=0.0; Deform[0][1]=0.0;Deform[1][0]=0.0;Deform[1][1]=0.0;
         for(m=1;m<nbr;m++){
            
            nbr_of_n = p[nbr_of_i].nn_static[m]; //  (M is nbr of nbr of I)
            rij_x = pbcx( p[nbr_of_n].ri[0] - p[nbr_of_i].ri[0] );
            rij_y = pbcy( p[nbr_of_n].ri[1] - p[nbr_of_i].ri[1] );
            Rij_x = pbcx( p[nbr_of_n].R[0]  - p[nbr_of_i].R[0]  );
            Rij_y = pbcy( p[nbr_of_n].R[1]  - p[nbr_of_i].R[1]  );
         
            Deform[0][0] += rij_x*Rij_x;
            Deform[0][1] += rij_x*Rij_y;
            Deform[1][0] += rij_y*Rij_x;
            Deform[1][1] += rij_y*Rij_y;
         }
      
         Deform[0][0] = (Deform[0][0]*sqRsum1)  - 1.0; Deform[0][1] = (Deform[0][1]*sqRsum2)  - 0.0;
         Deform[1][0] = (Deform[1][0]*sqRsum1)  - 0.0; Deform[1][1] = (Deform[1][1]*sqRsum2)  - 1.0;
   
         rij_x = pbcx( p[i].ri[0] - p[nbr_of_i].ri[0]  );
         rij_y = pbcy( p[i].ri[1] - p[nbr_of_i].ri[1]  );
         Rij_x = pbcx( p[i].R[0]  - p[nbr_of_i].R[0]   );
         Rij_y = pbcy( p[i].R[1]  - p[nbr_of_i].R[1]   );
         
         sumx += (4.0*rij_x) - 2.0*( ((Deform[0][0]+1.0)*Rij_x) + ((Deform[0][1])*Rij_y));
         sumy += (4.0*rij_y) - 2.0*( ((Deform[1][0])*Rij_x) + ((Deform[1][1]+1.0)*Rij_y));
      }
      // p[i].forceChi[0] = (hchi*sumx); p[i].forceChi[1] = (hchi*sumy);
      p[i].forceChi[0] = hchi_0*(p[i].h_xy*sumx); p[i].forceChi[1] = hchi_0*(p[i].h_xy*sumy);
   }
   

}
void h_chi_force_warning(){
   	printf(RED"\n Error ! Wrong choice of h_chi_force in Global.c|force.c\n"RESET);
	exit(0);
}
//#######################################################################################################################
