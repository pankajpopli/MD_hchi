#include"global.h"

void write_cord();

void generate_latt()
{
	
	// Predefined lattice vectors 
   double latt_vector1[2]={a,0.0};
	double latt_vector2[2]={0.0,a}; 
	// varbls for latt gnrtn
	int n1=0,n2=0,m=0;
	int nx,ny,pair_id;
	// varlbs for nearN printing
	double r=0.0,dummyx=0.0,dummyy=0.0,dummy2x=0.0,dummy2y=0.0;
	int i=0,count;
	
	// code
   nx = 2*part_id + 2;
   ny = nx;
   lbox_x = nx*latt_vector1[0];
   lbox_y = ny*latt_vector2[1];
   printf(GRN"\n the half length of box (Lx/2, Ly/2) = %.15lf,%.15lf and a = %.15lf\n"RESET,lbox_x*0.5,lbox_y*0.5,a);

   for(n1=0;n1<nx;n1++){ //square
      for(n2=0;n2<ny;n2++){
         pair_id = (n1 + nx*n2); //printf("%d\n",pair_id);
         p[pair_id].pid	= pair_id;
         p[pair_id].R[0] = -(lbox_x*0.5) + n1*a + a ;	
         p[pair_id].R[1] = -(lbox_y*0.5) + n2*a + a ;	
         p[pair_id].ri[0]=p[pair_id].R[0];
         p[pair_id].ri[1]=p[pair_id].R[1];
         p[pair_id].ri0[0]=p[pair_id].ri[0];
         p[pair_id].ri0[1]=p[pair_id].ri[1];
      }
   }

// Wrting co-ordinates in a file
   write_cord();  

// This part of the code checks whether we need to impose PBC again in other parts of the code- 
   for(m=0;m<totalN;m++){
      //periodicity in x
      if(abs(p[m].R[0])>lbox_x*0.5){
         printf(MAG"\n warning! need periodic boundary condition %d\n"RESET,1);
      }
      if(abs(p[m].R[1])>lbox_y*0.5){
         printf(MAG"\n warning! need periodic boundary condition %d\n"RESET,2);
      }
   }
}
   
 void write_cord()
{
	FILE *f1;
	int pair_no=0;
	
	f1 = fopen("cordinates.dat","w");
	for(pair_no=0;pair_no<totalN;pair_no++){
		fprintf(f1, "%d\t%.15lf\t%.15lf\n",pair_no,p[pair_no].R[0], p[pair_no].R[1]);
   }fclose(f1);
   
   f1 = fopen("system_profile.dat","a");
   fprintf(f1,"lbox_x = %lf and lbox_y = %lf.",lbox_x,lbox_y);
   fclose(f1);
}
