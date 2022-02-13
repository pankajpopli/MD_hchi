//######################################velcoity verlet rotor lattice cavagna########################################//
#include "global.h"



int main(){
/*::::::::::::::::::::::::::::::::::::::::::::::::clock:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
   clock_t start, end;
	double cpu_time_used;
	start = clock();
	if(system("clear")==0){};

/*::::::::::::::::::::::::::::::::::::::::::::::::variables::::::::::::::::::::::::::::::::::::::::::::::::*/
	setvbuf(stdout, NULL, _IONBF, 0);
	FILE *f1,*f2,*f3,*f;
	int i=0,j=0,k=0;
	double globalX=0.0,totalChi=0.0,energy_i=0.0;
   double smallb=0.0,smalla=0.0,noise_x=0.0,noise_y=0.0,totalFx=0.0,totalFy=0.0,pe=0.0,ke=0.0;
   double hchi_counter=0;
	char filename[500];

/*::::::::::::::::::::::::::::::::::::::::::::::::initi PRNG::::::::::::::::::::::::::::::::::::::::::::::::*/
   srand (time(NULL));
   //srand(2230);
	const gsl_rng_type *type;
	gsl_rng  *r1,*r2,*r3,*r4,*r5;
   gsl_rng_env_setup ();
	type = gsl_rng_default;
   r1 = gsl_rng_alloc (type);
   r2 = gsl_rng_alloc (type);
   r3 = gsl_rng_alloc (type);
   gsl_rng_set (r1, 3454);
   gsl_rng_set (r2, 7648);
   gsl_rng_set (r3, 8903);

/*::::::::::::::::::::::::::::::::::::::::::::::::initialization::::::::::::::::::::::::::::::::::::::::::::::::*/
	/*order matters */
	initializer();		      //|1
	generate_latt();        //|2
   n_list_static();        //|3
	// generate_velo(r4);   //|4
	create_c_list();        //|5
	init_Y_mat();           //|6
   spatial_random_h(r3);   //|7
/*:::::::::::::::::::::::::::::::::::::::::::::::pre-calculation----::::::::::::::::::::::::::::::::::::::::::::::::*/

   totalChi = X();   //functio X() also updates chi for each particle
   energy_i	  = c_list_totalE("oldCord");

/*:::::::::::::::::::::::::::::::::::::::::::::::integration of EOM:::::::::::::::::::::::::::::::::::::::::::::::*/

   //-------------------remove and create new time folder and log file------------------->>
      if(system("rm -rf time")!=0){printf("time/configs not deleted\n");exit(0);}
      remove("log");
      mkdir("time",0777);
      f2 = fopen("log","a");
      fprintf(f2,"#step\tPE\tTemp\tTE\tX/N\th_0\n");
      printf(GRN"\n doing MD simulation (N=%d, V=%lf, T=%lf) Total steps = %d with printf freq = %d \n\n"RESET,totalN,lbox_x*lbox_y,red_T,nsteps,printFreq);
      printf(YEL"final hchi_0 = %lf, plot change time = %d\n\n"RESET, (nsteps/h_step_intrvl)*h_step,h_step_intrvl/printFreq);
   //-------------------real thing start here--------------------------------------------->>
      // initializer.c file setsup everything for rect. lattice, hchi_force, nbrs_list, etc
      // we modify only the initial positions of the particles to triangular lattice below.
      int n1,n2,nx,ny,pair_id;
      double ac,rhoc;
      FILE *f_tring_eq;
         nx = 2*part_id + 2;
         ny = nx;
         rhoc		=	sqrt(((sqrt(3.0))*red_rho)/2.0); 
         ac			=	sigma/rhoc;
         for(n1=0;n1<nx;n1++){
            for(n2=0;n2<ny;n2++){
               pair_id = (n1 + nx*n2);
               p[pair_id].ri[0] = -(lbox_x*0.5) + n1*ac;
               p[pair_id].ri[1] = -(lbox_y*0.5) + n2*ac*sqrt(3.0/4.0);
               if(n2%2 != 0 ){
                  p[pair_id].ri[0] = p[pair_id].ri[0] + 0.5*ac;
               }
            }
         }
         // f_tring_eq = fopen("../lib/f_tring_eq.dat", "r");
         // for(n1=0;n1<nx;n1++){
         //    for(n2=0;n2<ny;n2++){
         //       pair_id = (n1 + nx*n2);
         //       fscanf(f_tring_eq, "%*d\t%lf\t%lf\t%lf\n",  &p[pair_id].ri[0],&p[pair_id].ri[1],&p[pair_id].chi);
         //    }
         // }fclose(f_tring_eq);

      //-------------------noise standard deiation and integrator params------------------>>
         stdDev = sqrt(2.0*damp*K_B*red_T*deltat);
         smallb = (1.0/(1.0 +  ((damp*deltat)/(2.0*mass))));
         smalla = (1.0 - ((damp*deltat)/(2.0*mass)))/(1.0 + ((damp*deltat)/(2.0*mass)));

      //--------------------calc forces and vel at time 0---------------------------------->>
         c_list_forceC();
         h_chi_force();
      //-------------------verlet integration (time steps)---------------------------------->>
         for(i=0;i<nsteps;i++){
               
            // increase hchi_0 at every 10000 sim steps
            // therefore in data files, hchi_0 will increase at each 100th  step
               if(hchi_0<(0.99999999)){
                  if(i%h_step_intrvl==0){
                     hchi_0 += h_step;
                     // if(hchi_counter==0){
                     //    hchi_counter=1;
                     // }
                  }
               }
            //-------------------set openMP for speed--------------------------------------->>
               //omp_set_dynamic(0);
               //omp_set_num_threads(thread_num);
               // #pragma omp parallel for private(totalFx,totalFy) shared(p,damp,mass,deltat)
            //-------------------erlet integration first half update------------------------->>
               for(k=0;k<totalN;k++){
                  noise_x = gsl_ran_gaussian_ziggurat(r1,stdDev);
                  noise_y = gsl_ran_gaussian_ziggurat(r2,stdDev);

                  totalFx = p[k].forceChi[0] + p[k].forceC[0] ;
                  totalFy = p[k].forceChi[1] + p[k].forceC[1];

                  //-------------------r_n+1 calculated using r_n and f_n , v_n-------------->>
                     p[k].ri[0] = pbcx((p[k].ri[0]) + (smallb*deltat*(p[k].vel[0])) + (((smallb*deltat*deltat)/(2.0*mass))*(totalFx)) + (((smallb*deltat)/(2.0*mass))*noise_x));
                     p[k].ri[1] = pbcy((p[k].ri[1]) + (smallb*deltat*(p[k].vel[1])) + (((smallb*deltat*deltat)/(2.0*mass))*(totalFy)) + (((smallb*deltat)/(2.0*mass))*noise_y));

                  //-------------------v_n+1/2 caclulated using v_n and f_n ------------------>>
                     p[k].vel[0] = (smalla*p[k].vel[0] + ((deltat/(2.0*mass))*smalla*(totalFx))  + ((smallb/mass)*noise_x));
                     p[k].vel[1] = (smalla*p[k].vel[1] + ((deltat/(2.0*mass))*smalla*(totalFy))  + ((smallb/mass)*noise_y));
                  //--------------------------------------------------------------------------------------------------------
               }
            //-------------------verlet integration second half update------------------------->>
               c_list_forceC();
               h_chi_force();
               update_c_list_full();
               //-------------------set openMP for speed--------------------------------------->>
                  //omp_set_dynamic(0);
                  //omp_set_num_threads(thread_num);
                  //#pragma omp parallel for private(totalFx,totalFy) shared(p,damp,mass,deltat)
               //-------------------f_n+1 forces calculated using r_n+1 ------------------------>>
                  for(k=0;k<totalN;k++){
                     totalFx = p[k].forceChi[0] + p[k].forceC[0] ;
                     totalFy = p[k].forceChi[1] + p[k].forceC[1];
                     //-------------------v_n+1 cacluated using v_n+1/2 and f_n+1 ------------------->>
                        p[k].vel[0] = (p[k].vel[0] + ((deltat/(2.0*mass))*(totalFx)) );
                        p[k].vel[1] = (p[k].vel[1] + ((deltat/(2.0*mass))*(totalFy)) );
                     //--------------------------------------------------------------------------------
                  }
               //-------------------------------------------------------
            //-------------------rint data in files---------------------------------------------->>
               if(  (i % printFreq) == 0  ){
                  pe = (c_list_totalE("newCord"))/totalN;
                  ke=calc_Temp();
                  globalX = X()/totalN;
                  fprintf(f2,"%d\t% .8lf\t% .8lf\t% .8lf\t% .8lf\t% .8lf\n",i/printFreq, pe,ke,pe+ke,globalX,hchi_0);
                  fflush(f2);
                  //-------------------Write config in a file-------------------------------------->>
                     sprintf(filename,"./time/new_cordinates_%d.dat",i/printFreq);
                     f1 = fopen(filename,"w");
                     for(j=0;j<totalN;j++){
                        fprintf(f1, "%d\t% .15lf\t% .15lf\t% .15lf\t% .15lf\t% .15lf\t% .15lf\t% .15lf\n",j,p[j].ri[0], p[j].ri[1],p[j].vel[0], p[j].vel[1],p[j].forceChi[0],p[j].forceChi[1], p[j].chi);
                     }
                     fclose(f1);
                  //-------------------print log on screen------------------------------------------>>
                  if(strcmp(print_config,"yes")==0){
                     printf("Itr No\t Potential E\t Kinetic E\t Total E\t  X global\t  h_chi_0\n\033[K %d\t% .8lf\t% .8lf\t% .8lf\t% .8lf\t% .6lf\033[A\r",i/printFreq, pe,ke,pe+ke,globalX,hchi_0);
                  }
               }
            //-----------------------------------------------------------------------------------------

      //-----------------------------------------------------------------------------------------------
   }
   end = clock();
   cpu_time_used = ( (double) (end - start))/CLOCKS_PER_SEC;
   printf(GRN"\n\n total run time = %.8f sec\n\n"RESET,cpu_time_used);
   return 0;
}
