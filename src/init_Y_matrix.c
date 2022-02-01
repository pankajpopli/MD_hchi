#include"global.h"


////////////////// local function declearations ///////////////////////// add here
void init_Y_mat_warning();
void init_Y_mat_triang();	
void init_Y_mat_sq();		
void init_Y_mat_readin();
void init_Y_mat_honeycomb();
void init_Y_mat_kagome();
void init_Y_mat_rectangle();

////////////////////////// pointer array function decleration ///////////////////////////////////////////////////////////////////////////// add here
void(*initYmat_ptr_arr[])() = {init_Y_mat_warning,init_Y_mat_triang,init_Y_mat_sq,init_Y_mat_readin,init_Y_mat_honeycomb,init_Y_mat_kagome,init_Y_mat_rectangle};


////////////////////////////////////// functions to be called from main.c ///////////////////////////////
// {M1]
void init_Y_mat()
{
	(*initYmat_ptr_arr[initYmat_slctr])();
}

//[M1R*]
// All these functions are defined in seperate respective files

//////////////////////////////////////////////// warnings and readin functions definitions /////////////////////////////////
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

////////////////// functions definitions ////////////////////////////
void init_Y_mat_sq()
{
   printf(MAG"\n Y matrix for square is initialised\n"RESET);
	int nbr_sq = 8;
	double Y1=0.0,Y2=0.0,Y3=0.0,Y4=0.0,Rn0x=0.0,Rn0y=0.0;
	int N=0,h=0,k=0,nh=0;
	
		N = sqrt(totalN);
	
		for( h=0;h<totalN;h++)
		{
			Y1=0.0,Y2=0.0,Y3=0.0,Y4=0.0;
			for( k=1;k<nbr_sq+1;k++)
			{	
				nh = p[h].nn_static[k];
				Rn0x=pbcx(p[nh].R[0] - p[h].R[0]); Rn0y=pbcy(p[nh].R[1] - p[h].R[1] );
						
				Y1 += (Rn0x)*(Rn0x);
				Y2 += (Rn0x)*(Rn0y);
				Y3 += (Rn0y)*(Rn0x);
				Y4 += (Rn0y)*(Rn0y);
			}
			Ymat[0]	=	Y1;
			Ymat[1]	=	Y2;
			Ymat[2]	=	Y3;
			Ymat[3]	=	Y4;
		}


}
//////////////////
void init_Y_mat_honeycomb()
{
	printf(MAG"\n Y matrix for honeycomb is initialised\n"RESET);
   int nbr_honey=9; 
	double Y1=0.0,Y2=0.0,Y3=0.0,Y4=0.0,Rn0x=0.0,Rn0y=0.0;
	int xi=0,yi=0,xi2=0,yi2=0,N=0,h=0,k=0,l=0,sit=0,sit_dmy1=0;
	
		N = sqrt(totalN);
		for(h=0;h<totalN;h++)
		{
         Y1=0.0;Y2=0.0;Y3=0.0;Y4=0.0;
			xi = h%N;
			yi = floor(h/N);
			for(l=0;l<2;l++)
			{	
				if( ((h+yi)%2==0)  )
				{
					sit = ((l+1)%2)*h + (l%2)*((xi +N- 1)%N + ((yi + N)%N)*N);
					xi2 = sit%N; yi2 = floor(sit/N);
				}else if( ((h+yi)%2==1)  )
				{
					sit = ((l+1)%2)*h + (l%2)*((xi +N+ 1)%N + ((yi + N)%N)*N);
					xi2 = sit%N; yi2 = floor(sit/N);
				}
				for(k=1;k<nbr_honey+1;k++)
				{
					sit_dmy1 = p[sit].nn_static[k];
					if(((sit+yi2)%2==0)&&(((xi2 +N- 1)%N + ((yi2 + N)%N)*N)==sit_dmy1))
					{continue;}
					Rn0x = pbcx(p[sit_dmy1].R[0] - p[sit].R[0]); Rn0y = pbcy(p[sit_dmy1].R[1] - p[sit].R[1]);
					Y1 += (Rn0x)*(Rn0x);
					Y2 += (Rn0x)*(Rn0y);
					Y3 += (Rn0y)*(Rn0x);
					Y4 += (Rn0y)*(Rn0y);
				}
			}
			Ymat[0]	=	Y1;
			Ymat[1]	=	Y2;
			Ymat[2]	=	Y3;
			Ymat[3]	=	Y4;
}
}
//////////////////
void init_Y_mat_kagome()
{
   printf(MAG"\n Y matrix for kagome is initialised\n"RESET);
	int nbr_kagome=8;
	printf(GRN" Intitalising  Y matrix now\n"RESET);
	double Y1=0.0,Y2=0.0,Y3=0.0,Y4=0.0,Rn0x=0.0,Rn0y=0.0;
	int xi=0,yi=0,xi2=0,yi2=0,N=0,h=0,k=0,l=0,sit=0,sit_dmy1=0,sit_lcl=0,sit_id=0;
   int k_start=0,basis=0;
   
      N = sqrt(totalN);
		for(h=0;h<totalN;h++)
		{
         Y1=0.0;Y2=0.0;Y3=0.0;Y4=0.0;
         sit = p[h].pid;
         sit_id = floor(sit/3); basis = sit%3;
			for(l=0;l<3;l++)
			{	
            sit_lcl = (3*sit_id) + l;     
            k_start=(l==0?1:(l==1?2:3));
				for(k=k_start;k<nbr_kagome+1;k++)
				{
					sit_dmy1 = p[sit_lcl].nn_static[k];
					Rn0x = pbcx(p[sit_dmy1].R[0] - p[sit_lcl].R[0]); Rn0y = pbcy(p[sit_dmy1].R[1] - p[sit_lcl].R[1]);
               // printf("%d\t%d\t%lf\t%lf\n",h,k,Rn0x,Rn0x*Rn0x+Rn0y*Rn0y);
					Y1 += (Rn0x)*(Rn0x);
					Y2 += (Rn0x)*(Rn0y);
					Y3 += (Rn0y)*(Rn0x);
					Y4 += (Rn0y)*(Rn0y);
				}
			}
			Ymat[0]	=	Y1;
			Ymat[1]	=	Y2;
			Ymat[2]	=	Y3;
			Ymat[3]	=	Y4;
     }
     //printf("hello\n%lf\t%lf\n%lf\t%lf\n",Ymat[0],Ymat[1],Ymat[2],Ymat[3]);
     //exit(0);
}
//////////////////
void init_Y_mat_triang()
{
	printf(MAG"\n Y matrix for triangle is initialised\n"RESET);
	int nbr_triang =  6;
	double Y1=0.0,Y2=0.0,Y3=0.0,Y4=0.0,Rn0x=0.0,Rn0y=0.0;
	int N=0,h=0,k=0,nsit=0;
	
	N = sqrt(totalN);
	
	for( h=0;h<totalN;h++)
	{
		Y1=0.0,Y2=0.0,Y3=0.0,Y4=0.0;
		for( k=1;k<nbr_triang+1;k++)
		{	
			nsit = p[h].nn_static[k];
			Rn0x=pbcx(p[nsit].R[0] - p[h].R[0]); Rn0y=pbcy(p[nsit].R[1] - p[h].R[1] );
						
			Y1 += (Rn0x)*(Rn0x);
			Y2 += (Rn0x)*(Rn0y);
			Y3 += (Rn0y)*(Rn0x);
			Y4 += (Rn0y)*(Rn0y);
		}
		Ymat[0]	=	Y1;
		Ymat[1]	=	Y2;
		Ymat[2]	=	Y3;
		Ymat[3]	=	Y4;
		// printf("y matrix\n%.15lf\t%.15lf\n%.15lf\t%.15lf\n",Y1,Y2,Y3,Y4);
	
	}

}
//////////////////
void init_Y_mat_rectangle()
{
	int nbr_rect =8;
	double Y1=0.0,Y2=0.0,Y3=0.0,Y4=0.0,Rn0x=0.0,Rn0y=0.0;
	int N=0,h=0,k=0,nh=0;
	
		N = sqrt(totalN);
	
		for( h=0;h<totalN;h++)
		{
			Y1=0.0,Y2=0.0,Y3=0.0,Y4=0.0;
			for( k=1;k<nbr_rect+1;k++)
			{	
				nh = p[h].nn_static[k];
				Rn0x=pbcx(p[nh].R[0] - p[h].R[0]); Rn0y=pbcy(p[nh].R[1] - p[h].R[1] );
						
				Y1 += (Rn0x)*(Rn0x);
				Y2 += (Rn0x)*(Rn0y);
				Y3 += (Rn0y)*(Rn0x);
				Y4 += (Rn0y)*(Rn0y);
			}
			Ymat[0]	=	Y1;
			Ymat[1]	=	Y2;
			Ymat[2]	=	Y3;
			Ymat[3]	=	Y4;
		}


}