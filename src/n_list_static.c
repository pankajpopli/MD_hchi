#include"global.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////
void sort_acc_r(double unSortArray[][3]); // auxialllry function



/////////////////////////////// local functions definition for pointer array /////////////////// add here
void n_list_static_warning();
void n_list_static_triang();
void n_list_static_square();
void n_list_static_readin();
void n_list_static_honeycomb();
void n_list_static_kagome();
void n_list_static_rectangle();
//////////////////////////////// definition of n_static_list pointer array of functions ///////// add here
void(*nlist_ptr_arr[])() = {n_list_static_warning,n_list_static_triang,n_list_static_square,n_list_static_readin,n_list_static_honeycomb,n_list_static_kagome,n_list_static_rectangle};
//////////////////////////////////////////////////////////////////////////////////////////////////





//##################################################################################################
//############################ function to be called from main #####################################
//##################################################################################################
// [M1]
void n_list_static()
{
	(*nlist_ptr_arr[n_list_slctr])();
}
//##################################################################################################






//##################################################################################################
//############################# routines to generate nbr static list  for ##########################
//##################################################################################################
// [R0] to generate warning @@
void n_list_static_warning()
{
	printf(" Error!, wrong square lattice chosen, n_list_static can't be created\n");
	exit(0);
}
// [R1] to create nbr static list for square @@
void n_list_static_square()
{
   printf(MAG"\n nbr list for square is initialised\n"RESET);
	int i=0,x=0,y=0;
	int N = sqrt(totalN);
	
	for(i=0;i<totalN;i++){
      x = i%N;
      y = floor(i/N);
      //NN[0]	=	i + N*j;
      p[i].nn_static[0] 	=	(x+N)%N + ((y+N)%N)*N;		// centre particle
      p[i].nn_static[1]	=	(x+N)%N + ((y+N+1)%N)*N;	//N
      p[i].nn_static[2]	=	(x+N)%N + ((y+N-1)%N)*N;	//S
      p[i].nn_static[3]	=	(x+N+1)%N + ((y+N)%N)*N;	//E
      p[i].nn_static[4]	=	(x+N-1)%N + ((y+N)%N)*N;	//W
      p[i].nn_static[5]	=	(x+N+1)%N + ((y+N+1)%N)*N;	//NE
      p[i].nn_static[6]	=	(x+N-1)%N + ((y+N+1)%N)*N;	//NW
      p[i].nn_static[7]	=	(x+N+1)%N + ((y+N-1)%N)*N;	//SE
      p[i].nn_static[8]	=	(x+N-1)%N + ((y+N-1)%N)*N;	//SW
      p[i].nn_static[9]	=	0;
		
   }

}
// [R2] to create nbr static list for triang @@
void n_list_static_triang()
{
   printf(MAG"\n nbr list for triangle is initialised\n"RESET);
  int i=0,x=0,y=0;
  int N = sqrt(totalN);
  
  for(i=0;i<totalN;i++){	
      x = i%N;
      y = floor(i/N);
      p[i].nn_static[0] =	(x+N)%N + ((y+N)%N)*N;	// center particle
      
      // generating and re-aranging the nbrs accordng to diagram in genrt_latt2.c
      if(y % 2 == 0){ 
         p[i].nn_static[1] = ((x+N) % N  + ((y+N+1)%N)*N);    //NE
         p[i].nn_static[5] = ((x+N) % N  + ((y+N-1)%N)*N);    //SE  
         p[i].nn_static[2] = ((x+N-1)%N  + ((y+N+1)%N)*N);    //NW
         p[i].nn_static[4] = ((x+N-1)%N  + ((y+N-1)%N)*N);    //SW    
      }else{
         p[i].nn_static[1] = ((x+N+1)%N + ((y+N+1)%N)*N);   //NE
         p[i].nn_static[5] = ((x+N+1)%N + ((y+N-1)%N)*N);   //SE
         p[i].nn_static[2] = ((x+N) % N + ((y+N+1)%N)*N);   //NW
         p[i].nn_static[4] = ((x+N) % N + ((y+N-1)%N)*N);   //SW
      }
      p[i].nn_static[3] = ((x+N-1)%N + ((y+N)%N)*N);  		 //W
      p[i].nn_static[6] = ((x+N+1)%N + ((y+N)%N)*N);  		 //E

      p[i].nn_static[7] = 0; //bcoz triangular lattice has 6 nbrs + 1 centre particle
      p[i].nn_static[8] = 0;
      p[i].nn_static[9] = 0;
   }	
}

// [R3] to create nbr static list for readin @@
void n_list_static_readin()
{
	printf(RED"\n Stattic nbr list, Not yet supported, contact pankajp@tifrh.res.in\n"RESET);
	exit(0);
}

// [R4] to create nbr static list for hanoeycomb @@
void n_list_static_honeycomb()
{
   printf(MAG"\n nbr list for honeycomb is initialised\n"RESET);
   int i=0,xi=0,yi=0,N=0;
   
   N = sqrt(totalN);
   for(i=0;i<totalN;i++){
      
      xi = i%N;
      yi = floor(i/N);
      p[i].nn_static[0] = (xi+N)%N + ((yi+N)%N)*N;
      
      if((i+yi)%2==0){
         p[i].nn_static[1] = (xi+N)%N + ((yi+N-1)%N)*N;
         p[i].nn_static[2] = (xi+N)%N + ((yi+N+1)%N)*N;
         p[i].nn_static[3] = (xi+N-1)%N + ((yi+N)%N)*N;
      }else if((i+yi)%2==1){
         p[i].nn_static[1] = (xi+N)%N + ((yi+N-1)%N)*N;
         p[i].nn_static[2] = (xi+N+1)%N + ((yi+N)%N)*N;
         p[i].nn_static[3] = (xi+N)%N + ((yi+N+1)%N)*N;
      }
      p[i].nn_static[4] = (xi+N+1)%N + ((yi+N+1)%N)*N;
      p[i].nn_static[5] = (xi+N+1)%N + ((yi+N-1)%N)*N;
      p[i].nn_static[6] = (xi+N-1)%N + ((yi+N+1)%N)*N;
      p[i].nn_static[7] = (xi+N-1)%N + ((yi+N-1)%N)*N;
      p[i].nn_static[8] = (xi+N)%N   + ((yi+N+2)%N)*N;
      p[i].nn_static[9] = (xi+N)%N   + ((yi+N-2)%N)*N;
   }
}
// [R1] to create nbr static list for kagome @@
void n_list_static_kagome()
{
   printf(MAG"\n nbr list for kagome is initialised\n"RESET);
   int id_prime=0,row=0,col=0,basis=0,Nx_t=0,Ny_t=0,i=0,id_part=0;
   
   Nx_t = sqrt(totalN)/2;
   Ny_t = (4*Nx_t)/3;
   
   for(i=0;i<totalN;i++){  
      
      id_part=p[i].pid;
      id_prime = floor(id_part/3);
      basis = id_part%3;
      row=id_prime%(Nx_t);
      col=floor(id_prime/(Nx_t));
      
      if(basis==0){
         p[i].nn_static[0]=3*(((row+Nx_t)%Nx_t) + (((col+Ny_t)%Ny_t)*Nx_t) ); //self
         p[i].nn_static[1]=3*(((row+Nx_t)%Nx_t) + (((col+Ny_t)%Ny_t)*Nx_t) )+1;
         p[i].nn_static[2]=3*(((row+Nx_t)%Nx_t) + (((col+Ny_t)%Ny_t)*Nx_t))+2;
         p[i].nn_static[3]=3*((row+Nx_t-1)%Nx_t + ((col+Ny_t)%Ny_t)*Nx_t)+1;
         p[i].nn_static[4]=3*((row+Nx_t-1+((col%2)*1))%Nx_t  + ((col+Ny_t-1)%Ny_t)*Nx_t)+2; 
         p[i].nn_static[5]=3*((row+Nx_t-1)%Nx_t + ((col+Ny_t)%Ny_t)*Nx_t)+2;
         p[i].nn_static[6]=3*((row+Nx_t-1+((col%2)*1))%Nx_t  + ((col+Ny_t-1)%Ny_t)*Nx_t)+1; 
         p[i].nn_static[7]=3*((row+Nx_t-1+((col%2)*1))%Nx_t  + ((col+Ny_t+1)%Ny_t)*Nx_t)+1;
         p[i].nn_static[8]=3*((row+Nx_t+((col%2)*1))%Nx_t  + ((col+Ny_t-1)%Ny_t)*Nx_t)+2;
         p[i].nn_static[9]=0;
      }else if(basis==1){
         p[i].nn_static[0]=3*(((row+Nx_t)%Nx_t) + (((col+Ny_t)%Ny_t)*Nx_t) )+1; //self 
         p[i].nn_static[1]=3*(((row+Nx_t)%Nx_t) + (((col+Ny_t)%Ny_t)*Nx_t) );
         p[i].nn_static[2]=3*(((row+Nx_t)%Nx_t) + (((col+Ny_t)%Ny_t)*Nx_t))+2;
         p[i].nn_static[3]=3*((row+Nx_t+((col%2)*1))%Nx_t  + ((col+Ny_t-1)%Ny_t)*Nx_t)+2;
         p[i].nn_static[4]=3*((row+Nx_t+1)%Nx_t + ((col+Ny_t)%Ny_t)*Nx_t); 
         p[i].nn_static[5]=3*((row+Nx_t-1+((col%2)*1))%Nx_t  + ((col+Ny_t-1)%Ny_t)*Nx_t)+2;
         p[i].nn_static[6]=3*((row+Nx_t+((col%2)*1))%Nx_t  + ((col+Ny_t+1)%Ny_t)*Nx_t);
         p[i].nn_static[7]=3*((row+Nx_t+1)%Nx_t + ((col+Ny_t)%Ny_t)*Nx_t)+2;
         p[i].nn_static[8]=3*((row+Nx_t+((col%2)*1))%Nx_t  + ((col+Ny_t-1)%Ny_t)*Nx_t);
         p[i].nn_static[9]=0;
      }else if(basis==2){
         p[i].nn_static[0]=3*(((row+Nx_t)%Nx_t) + (((col+Ny_t)%Ny_t)*Nx_t) )+2;	//self
         p[i].nn_static[1]=3*(((row+Nx_t)%Nx_t) + (((col+Ny_t)%Ny_t)*Nx_t) );
         p[i].nn_static[2]=3*(((row+Nx_t)%Nx_t) + (((col+Ny_t)%Ny_t)*Nx_t))+1;
         p[i].nn_static[3]=3*((row+Nx_t-1+((col%2)*1))%Nx_t  + ((col+Ny_t+1)%Ny_t)*Nx_t)+1;
         p[i].nn_static[4]=3*((row+Nx_t+((col%2)*1))%Nx_t  + ((col+Ny_t+1)%Ny_t)*Nx_t);
         p[i].nn_static[5]=3*((row+Nx_t-1+((col%2)*1))%Nx_t  + ((col+Ny_t+1)%Ny_t)*Nx_t); 
         p[i].nn_static[6]=3*((row+Nx_t+((col%2)*1))%Nx_t  + ((col+Ny_t+1)%Ny_t)*Nx_t)+1;
         p[i].nn_static[7]=3*((row+Nx_t-1)%Nx_t + ((col+Ny_t)%Ny_t)*Nx_t)+1;
         p[i].nn_static[8]=3*((row+Nx_t+1)%Nx_t + ((col+Ny_t)%Ny_t)*Nx_t);
         p[i].nn_static[9]=0;
      }else{
         printf(RED"some problem in finding nbr of kagome lattice\n"RESET);
         exit(0);
      }
   }
}
void n_list_static_rectangle()
{
   printf(MAG"\n nbr list for rectangle is initialised\n"RESET);
	int i=0,x=0,y=0;
	int N = sqrt(totalN);
	
	for(i=0;i<totalN;i++){
      x = i%N;
      y = floor(i/N);
      //NN[0]	=	i + N*j;
      p[i].nn_static[0] 	=	(x+N)%N + ((y+N)%N)*N;		// centre particle
      p[i].nn_static[1]	=	(x+N)%N + ((y+N+1)%N)*N;	//N
      p[i].nn_static[2]	=	(x+N)%N + ((y+N-1)%N)*N;	//S
      p[i].nn_static[3]	=	(x+N+1)%N + ((y+N)%N)*N;	//E
      p[i].nn_static[4]	=	(x+N-1)%N + ((y+N)%N)*N;	//W
      p[i].nn_static[5]	=	(x+N+1)%N + ((y+N+1)%N)*N;	//NE
      p[i].nn_static[6]	=	(x+N-1)%N + ((y+N+1)%N)*N;	//NW
      p[i].nn_static[7]	=	(x+N+1)%N + ((y+N-1)%N)*N;	//SE
      p[i].nn_static[8]	=	(x+N-1)%N + ((y+N-1)%N)*N;	//SW
      p[i].nn_static[9]	=	0;
		
   }

}
//##################################################################################################




//#########################################################################################################################################################
//######################### Auxillary functions required by this file #####################################################################################
//#########################################################################################################################################################

void sort_acc_r(double unSortArray[][3])
{
	int g=0,h=0;
   
	for( g=1;g<totalN;g++){
      for( h=0; h<(totalN-1); h++){
         if(unSortArray[h][2]>unSortArray[h+1][2]){
            double temp=unSortArray[h][2];
            unSortArray[h][2]=unSortArray[h+1][2];
            unSortArray[h+1][2]=temp;

            double temp1 = unSortArray[h][0];
            unSortArray[h][0] = unSortArray[h+1][0];
            unSortArray[h+1][0]=temp1;
         
            double temp2 = unSortArray[h][1];
            unSortArray[h][1] = unSortArray[h+1][1];
            unSortArray[h+1][1]=temp2;
		   }
      }
   }
}
//##################################################################################################