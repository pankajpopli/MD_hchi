//17-01-2022---> adding rectangular lattice
//15-03-2017---> using the standard method of generating lattices with centre particle always at zero.
/*
		HEXAGONE and Square LATTICE AND CONVENTION OF INDEX NO.



																				3
			2								1				4	@@ ####### @@ ####### @@ 2
			  @@ ############### @@ 					##         *          ##
			 ##               #  ##						##       *     *      ##
			##             #      ##					##     *         *    ##
		  ##             #        ##					##   *       0     *  ##
		 ##            #  theta    ##				5	@@ *        @@      * @@ 1
	3	 @@ ######### @@ ######### @@	6				##   *            *   ##
		##              0          ##					##     *        *     ##
		 ##                       ##					##       *    *       ##
			##                    ##					##          *         ##
			 ##                  ##					6	@@ ####### @@ ####### @@
			  @@ ############## @@										7				8,			honeycomb      kagome   rectangular
			4								5
*/

#include"global.h"


//////////////////////////// Local functions defntn ////////////////////////////// add here
void latt_warning();
void latt_triang();
void latt_square();
void latt_readin();
void latt_honeycomb();
void latt_kagome();
void latt_rectangle();
//
void patch_zigzag();
void patch_shear();
void patch_warning();
/////////////////////////////////////////////////////////////////////////////////
void write_cord();

///////////////////// function pointers ///////////////////////////////////////// add here
void (*latt_ptr_arr[])()  = {latt_warning,latt_triang, latt_square,latt_readin,latt_honeycomb,latt_kagome,latt_rectangle};
// void (*patch_ptr_arr[])() = {patch_warning,patch_zigzag, patch_shear};
////////////////////////////////////////////////////////////////////////////////


///////////////////// functions to be called from Main /////////////////////////
//[M1] to generate a lattice
void generate_latt()
{
	(*latt_ptr_arr[latt_slctr])();

}
// [M2] to create a patch in a lattice
// void create_patch()
// {
// 	(*patch_ptr_arr[patch_slctr])();
// }
/////////////////////////////////////////////////////////////////////////////////

typedef struct Nbour{
		double r[2];
}nbour;
nbour nn[];



//#########################################################################################################################################################
//############################### function code to generate lattices ######################################################################################
//######################################################################################################################################################### add here

// @@ [GL0] to generate wrong lattice warning @@@@@@@@@@
void latt_warning()
{
	printf(RED"\n Error! Lattice structure entered is not supported, enter correct string in Global.c."RESET);
	exit(0);
}

// @@ [GL1] routine to generate triangular lattice @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void latt_triang()
{
   // nbour nn[6];
	printf(MAG"\n triangle lattice is initialised\n"RESET);
   //Predefined lattice vectors
	double latt_vector1[2]={a,0.0};
	double latt_vector2[2]={(a/2.0),sqrt(3.0/4.0)*a};
	//varbls for latt gnrtn
	int n1=0,n2=0,m=0;
	int nx,ny,pair_id;
	// varbls for nearN printing
	double r=0.0,dummyx=0.0,dummyy=0.0,dummy2x=0.0,dummy2y=0.0;
	int i=0,count=0;

	// code
	nx = 2*part_id + 2;
	ny = nx;
   lbox_x = nx*latt_vector1[0];
   lbox_y = ny*latt_vector2[1];
   printf(GRN"\n the half length of box (Lx/2, Ly/2) = %.15lf,%.15lf and a = %.15lf\n"RESET,lbox_x*0.5,lbox_y*0.5,a);

   for(n1=0;n1<nx;n1++){
      for(n2=0;n2<ny;n2++){
         pair_id = (n1 + nx*n2);
         p[pair_id].pid	= pair_id;
         p[pair_id].ri0[0] = -(lbox_x*0.5) + n1*latt_vector1[0] + 0*latt_vector2[0];
         p[pair_id].ri0[1] = -(lbox_y*0.5) + n2*latt_vector2[1] + 0*latt_vector2[1];
         p[pair_id].R[0] = p[pair_id].ri0[0];
         p[pair_id].R[1] = p[pair_id].ri0[1];
         p[pair_id].ri[0]=p[pair_id].R[0];
         p[pair_id].ri[1]=p[pair_id].R[1];
         if(n2%2 != 0 ){
            p[pair_id].ri0[0] = p[pair_id].ri0[0] + 0.5*a;
            p[pair_id].R[0] = p[pair_id].ri0[0];
            p[pair_id].ri[0]=p[pair_id].R[0];
         }
      }
   }
/////////
	// print near Neighbour
///////
   // Wrting co-ordinates in a file
   write_cord();

   // This part of the code checks whether we need to impose PBC again in other parts of the code
   for(m=0;m<totalN;m++){
      //periodicity in x and y
      if(abs(p[m].ri0[0])>lbox_x*0.5){
         printf(MAG"\n warning! need periodic boundary condition %d\n"RESET,1);
      }
      if(abs(p[m].ri0[1])>lbox_y*0.5){
         printf(MAG"\n warning! need periodic boundary condition %d\n"RESET,2);
      }
   }
}

// @@ [GL2] to generate square lattice @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void latt_square()
{
   printf(MAG"\n square lattice is initialised\n"RESET);
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

   for(n1=0;n1<nx;n1++){
      for(n2=0;n2<ny;n2++){
         pair_id = (n1 + nx*n2); //printf("%d\n",pair_id);
         p[pair_id].pid	= pair_id;
         p[pair_id].ri0[0] = -(lbox_x*0.5) + n1*a + a ;
         p[pair_id].ri0[1] = -(lbox_y*0.5) + n2*a + a ;
         p[pair_id].R[0] = p[pair_id].ri0[0];
         p[pair_id].R[1] = p[pair_id].ri0[1];
         p[pair_id].ri[0]=p[pair_id].R[0];
         p[pair_id].ri[1]=p[pair_id].R[1];
      }
   }

// Wrting co-ordinates in a file
   write_cord();

// This part of the code checks whether we need to impose PBC again in other parts of the code-
   for(m=0;m<totalN;m++){
      //periodicity in x
      if(abs(p[m].ri0[0])>lbox_x*0.5){
         printf(MAG"\n warning! need periodic boundary condition %d\n"RESET,1);
      }
      if(abs(p[m].ri0[1])>lbox_y*0.5){
         printf(MAG"\n warning! need periodic boundary condition %d\n"RESET,2);
      }
   }

}

// @@ [GL4] to generate honeycomb lattice @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void latt_honeycomb()
{
   printf(MAG"\n honeycomb lattice is initialised\n"RESET);
	//Predefined lattice vectors
		// not doint it this way
	//varbls for latt gnrtn
   int n1=0,n2=0,m=0;
	int nx=0,ny=0,pair_id=0,count=0;

	// varbls for nearN printing
		/*double r=0.0,dummyx=0.0,dummyy=0.0,dummy2x=0.0,dummy2y=0.0,dummy3x=0.0,dummy3y=0.0;
		int i=0,count=0;*/
	// code
   nx = 2*part_id + 2;
   ny = nx;
   lbox_x = (nx/2.0)*(3.0*a);
   lbox_y =  ny*sqrt(3.0/4.0)*a;
   printf(GRN"\n the half length of box (Lx/2, Ly/2) and a is = %.15lf,%.15lf and a = %.15lf\n"RESET,lbox_x*0.5,lbox_y*0.5,a);

   for(n1=0;n1<nx;n1++){
      for(n2=0;n2<ny;n2++){
         pair_id = (n2 + nx*n1);
         p[pair_id].pid	= pair_id;
         p[pair_id].ri0[1]	= (-lbox_y*0.5) +	n1*sqrt(3.0/4.0)*a;
         if((n1%2)!=0){
            p[pair_id].ri0[0]	=	(-lbox_x*0.5) + (n2+(n2%2)*((n2-1)/2)+ ((n2+1)%2)*(n2/2)+0.5)*a;
         }else{
            p[pair_id].ri0[0]	=	(-lbox_x*0.5)	+ (n2+(n2%2)*((n2+1)/2)+((n2+1)%2)*(n2/2))*a;
         }
         p[pair_id].R[0] = p[pair_id].ri0[0];
         p[pair_id].R[1] = p[pair_id].ri0[1];
         p[pair_id].ri[0]=p[pair_id].R[0];
         p[pair_id].ri[1]=p[pair_id].R[1];
         //printf("%d\t% .9lf\t% .9lf\n",count,R[count].x,R[count].y );
         count = count+1;
      }
   }

// Wrting co-ordinates in a file
   write_cord();

// This part of the code checks whether we need to impose PBC again in other parts of the code
   for(m=0;m<totalN;m++){
      //periodicity in x and y
      if(abs(p[m].ri0[0])>lbox_x*0.5)
      {
         printf(MAG"\n warning! need periodic boundary condition %d\n"RESET,1);
      }
      if(abs(p[m].ri0[1])>lbox_y*0.5)
      {
         printf(MAG"\n warning! need periodic boundary condition %d\n"RESET,2);
      }
   }
}
// @@ [GL5] to generate kagome lattice @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void latt_kagome()
{
   printf(MAG"\n kagome lattice is initialised\n"RESET);
	//varbls for latt gnrtn
   int n1=0,n2=0,n3,m=0,ilim,jlim;
	int nx=0,ny=0,count=0;

	// code
   nx = sqrt(totalN);
   ny = (2*nx)/3;
   ilim = nx/2;
   jlim = (2*nx)/3;
   lbox_x = nx*a*0.5;
   lbox_y = (nx/sqrt(3.0))*a;
   printf(GRN"\n the half length of box (Lx/2, Ly/2) and a is = %.15lf,%.15lf and a = %.15lf\n"RESET,lbox_x*0.5,lbox_y*0.5,a);

   for(n1=0;n1<jlim;n1++){
      for(n2=0;n2<ilim;n2++){
         p[count].ri0[0] = (-lbox_x*0.5)+ n2*a+((n1%2)*0.5*a);
         p[count].ri0[1] = (-lbox_y*0.5)+n1*sqrt(3.0/4.0)*a;
         p[count].R[0]   = p[count].ri0[0];
         p[count].R[1]   = p[count].ri0[1];
         p[count].pid=count;
         count+=1;
         for(n3=0;n3<2;n3++){
            p[count].ri0[0] = (-lbox_x*0.5)+(n2*a+((n1%2)*0.5*a))+a*((!(n3%2)*0.5) + ((n3%2)*0.25));
            p[count].ri0[1] = (-lbox_y*0.5)+(n1*sqrt(3.0/4.0)*a) + a*((!(n3%2)*0.0) + ((n3%2)*sqrt(3.0/4.0)*0.5));
            p[count].R[0]   = p[count].ri0[0];
            p[count].R[1]   = p[count].ri0[1];
            p[count].ri[0]  = p[count].R[0];
            p[count].ri[1]  = p[count].R[1];
            p[count].pid=count;
            count+=1;
         }
      }
   }
	// Wrting co-ordinates in a file
   write_cord();

	// This part of the code checks whether we need to impose PBC again in other parts of the code
   for(m=0;m<totalN;m++){
      //periodicity in x and y
      if(abs(p[m].ri0[0])>lbox_x*0.5){
         printf(MAG"\n warning! need periodic boundary condition %d\n"RESET,1);
      }
      if(abs(p[m].ri0[1])>lbox_y*0.5){
         printf(MAG"\n warning! need periodic boundary condition %d\n"RESET,2);
      }
   }
}

// @@ [GL3] for readin lattice @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void latt_readin()
{
	printf(RED"\n Lattice readin is not supported yet, enter correct string in Global.c\n"RESET);
	exit(0);
}


// @@ [GL6] to generate rectangular lattice @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void latt_rectangle()
{
   printf(MAG"\n rectangular lattice is initialised\n"RESET);
   printf(YEL"\n The rectangular lattice dimensions are chosen ax and ay = sqrt(3/4)*ax |\n This makes lattice commensurate with triangular lattice\n"RESET);

   double ax = a;
   double ay = a*sqrt(3.0/4.0);
	// varbls for latt gnrtn
	int n1=0,n2=0,m=0;
	int nx,ny,pair_id;
	// varlbs for nearN printing
	double r=0.0,dummyx=0.0,dummyy=0.0,dummy2x=0.0,dummy2y=0.0;
	int i=0,count;

	// code
   nx = 2*part_id + 2;
   ny = nx;
   lbox_x = nx*ax;
   lbox_y = ny*ay;
   printf(GRN"\n the half length of box (Lx/2, Ly/2) = %.15lf,%.15lf and a = %.15lf\n"RESET,lbox_x*0.5,lbox_y*0.5,a);

   for(n1=0;n1<nx;n1++){
      for(n2=0;n2<ny;n2++){
         pair_id = (n1 + nx*n2); 
         p[pair_id].pid	= pair_id;
         p[pair_id].ri0[0] = -(lbox_x*0.5) + n1*ax + ax ;
         p[pair_id].ri0[1] = -(lbox_y*0.5) + n2*ay + ay ;
         p[pair_id].R[0] = p[pair_id].ri0[0];
         p[pair_id].R[1] = p[pair_id].ri0[1];         
         p[pair_id].ri[0]=p[pair_id].R[0];
         p[pair_id].ri[1]=p[pair_id].R[1];
      }
   }

// Wrting co-ordinates in a file
   write_cord();

// This part of the code checks whether we need to impose PBC again in other parts of the code-
   for(m=0;m<totalN;m++){
      //periodicity in x
      if(abs(p[m].ri0[0])>lbox_x*0.5){
         printf(MAG"\n warning! need periodic boundary condition %d\n"RESET,1);
      }
      if(abs(p[m].ri0[1])>lbox_y*0.5){
         printf(MAG"\n warning! need periodic boundary condition %d\n"RESET,2);
      }
   }
}
//#####################################################################################################################################################//









//#########################################################################################################################################################
//######################################### Create patch in the system ####################################################################################
//######################################################################################################################################################### add here
// @@ [CP0] create warning for patch @@@@@@@@@@@@@@@@@@
// void patch_warning()
// {
// 	printf(GRN"\n No patch created\n"RESET);
// }
// @@ [CP1] to create zigzag patch in the system @@@@@@@@@
// void patch_zigzag()
// {
// 	int i=0,ix=0,iy=0,Nx=0;
// 	FILE *f1;

// 	Nx = sqrt(totalN);
// 	f1 = fopen("cord_patched.dat","w");

// 	for(i=0;i<totalN;i++){
// 		if(((p[i].R[1] > -(patch_size*a)) && p[i].R[1] <= (patch_size*a)) && ((p[i].R[0] > -(patch_size*a)) && p[i].R[0] <= (patch_size*a))){
// 			{
// 				ix = i%Nx; iy = floor(i/Nx);
// 				if((iy%2)==0){
// 					p[i].ri0[0] = p[i].ri0[0] + 0.5*a;
// 					p[i].ri0[0] = pbcx(p[i].ri0[0]);
// 					p[i].ri0[1] = pbcy(p[i].ri0[1]);
// 				}
// 			}
// 		}
// 		fprintf(f1,"%d\t%lf\t%lf\t%d\n",i,p[i].ri0[0],p[i].ri0[1],p[i].pid);
// 	}
// 	fclose(f1);
// 	printf(GRN"\n Patch created zigzag \n"RESET);
// }
// @@ [CP2] to create shear patch in the system @@@@@@@@@@@@@@
// void patch_shear()
// {
// 	int i=0,ix=0,iy=0,Nx=0;
// 	FILE *f1;

// 	Nx = sqrt(totalN);
// 	f1 = fopen("cord_patched.dat","w");

// 	for(i=0;i<totalN;i++){
// 		if(((p[i].R[1] > -(patch_size*a)) && p[i].R[1] <= (patch_size*a)) && ((p[i].R[0] > -(patch_size*a)) && p[i].R[0] <= (patch_size*a))){
// 			{
// 				ix = i%Nx; iy = floor(i/Nx);
// 				{
// 					p[i].ri0[0] = p[i].ri0[0] + sheer_x*(p[i].R[1]+(patch_size*a));
// 					p[i].ri0[1] = p[i].ri0[1] + sheer_y*(p[i].R[0]+(patch_size*a));
// 					p[i].ri0[0] = pbcx(p[i].ri0[0]);
// 					p[i].ri0[1] = pbcy(p[i].ri0[1]);
// 				}
// 			}
// 		}
// 		fprintf(f1,"%d\t%lf\t%lf\t%d\n",i,p[i].ri0[0],p[i].ri0[1],p[i].pid);
// 	}
// 	fclose(f1);
// 	printf(GRN"\n Patch created shear \n"RESET);
// }
//#######################################################################################################################################################





//#########################################################################################################################################################
//######################################### Put Strain in the system ######################################################################################
//#########################################################################################################################################################
// @@ [PS1] put strain in the system
void strain(double ex, double ey)
{
	FILE *f1;
	int pair_no;

	f1 = fopen("cordinates_strained.dat","w");
	{
		for(pair_no=0;pair_no<totalN;pair_no++){
			p[pair_no].ri0[0] = (1.0+ex)*p[pair_no].ri0[0];
			p[pair_no].ri0[1] = (1.0+ey)*p[pair_no].ri0[1];
			fprintf(f1, "%d\t%.15lf\t%.15lf\n",pair_no,p[pair_no].ri0[0], p[pair_no].ri0[1]);
      }
		lbox_x = (1.0 + ex )*lbox_x;
		lbox_y = (1.0 + ey )*lbox_y;
	}
	printf(WHT"\n the new box full length after strain is -> %.10lf\t%.10lf"RESET,lbox_x,lbox_y);
	fclose(f1);
}
//#######################################################################################################################################################





//#########################################################################################################################################################
//######################### Auxillary functions required by this file #####################################################################################
//#########################################################################################################################################################
// (aux 1)  -------------writing cordinate of particles in a file.---------------------------//
void write_cord()
{
	FILE *f1;
	int pair_no=0;

	f1 = fopen("cordinates.dat","w");
	for(pair_no=0;pair_no<totalN;pair_no++){
		fprintf(f1, "%d\t%.15lf\t%.15lf\n",pair_no,p[pair_no].ri0[0], p[pair_no].ri0[1]);
   }fclose(f1);

   f1 = fopen("system_profile.dat","a");
   fprintf(f1,"lbox_x = %lf\nlbox_y = %lf",lbox_x,lbox_y);
   fclose(f1);
}

// (aux)
//############################################################################################################################################################

























//--------------------------------------------------------END of Program--------------------------------------------------------------//
