/*
	This code may not work properly because nbr_honey = 9 instead of 10 please be carefull and sort this discrepancy
*/
#include"global.h"


double X_honeycomb()
{

	int i=0,j=0,k=0,l=0,xi=0,yi=0,xi2=0,yi2=0,neighbours=9,N=0,sit=0,sit_dmy1=0;
	double detY=0,detY_inv=0,rn0x=0.0,rn0y=0.0,Rn0x=0.0,Rn0y=0.0;
	double chi=0.0,chi_total=0.0,X11=0.0,X12=0.0,X21=0.0,X22=0.0,Y11=0.0,Y12=0.0,Y21=0.0,Y22=0.0,eps11=0.0,eps12=0.0,eps21=0.0,eps22=0.0;
	int nbr_honey=9;
	N = sqrt(totalN);

//-------------------------------
		for(j=0;j<totalN;j++)
		{
			X11 = X21 = X12 = X22 = 0;
         Y11 = Y21 = Y12 = Y22 = 0;
         chi = 0;
         
         xi = j%N; yi = floor(j/N);
			for(l=0;l<2;l++)
			{	
				if( ((j+yi)%2==0)  )
				{
					sit = ((l+1)%2)*j + (l%2)*((xi +N- 1)%N + ((yi + N)%N)*N);
					xi2 = sit%N; yi2 = floor(sit/N);
				}else if( ((j+yi)%2==1)  )
				{
					sit = ((l+1)%2)*j + (l%2)*((xi +N+ 1)%N + ((yi + N)%N)*N);
					xi2 = sit%N; yi2 = floor(sit/N);
				}
				for(k=1;k<nbr_honey;k++)
				{
					sit_dmy1=p[sit].nn_static[k];
					if(((sit+yi2)%2==0)&(((xi2 +N- 1)%N + ((yi2 + N)%N)*N)==sit_dmy1))
					{continue;}
					rn0x = pbcx(p[sit_dmy1].ri[0] - p[sit].ri[0]); rn0y = pbcy(p[sit_dmy1].ri[1] - p[sit].ri[1]);
					Rn0x = pbcx(p[sit_dmy1].R[0] - p[sit].R[0]);   Rn0y = pbcy(p[sit_dmy1].R[1] - p[sit].R[1]);
					
					X11 += (rn0x)*(Rn0x);
					X12 += (rn0x)*(Rn0y);
					X21 += (rn0y)*(Rn0x);
					X22 += (rn0y)*(Rn0y);
				
					Y11 += (Rn0x)*(Rn0x);
					Y12 += (Rn0x)*(Rn0y);
					Y21 += (Rn0y)*(Rn0x);
					Y22 += (Rn0y)*(Rn0y);
				}
			}
			
			detY = (Y11*Y22) - (Y12*Y21) ; detY_inv = 1.0/detY;
			
			eps11	=	 ( X11*Y22 - X12*Y12)*detY_inv -1;
			eps12	=	 (-X11*Y21 + X12*Y11)*detY_inv;
			eps21	=	 ( X21*Y22 - X22*Y12)*detY_inv;
			eps22	=	 (-X21*Y21 + X22*Y11)*detY_inv -1;
			
			for(k=1;k<nbr_honey;k++)
			{
				sit_dmy1 = p[j].nn_static[k];
				rn0x = pbcx(p[sit_dmy1].ri[0]-p[j].ri[0]); rn0y = pbcy(p[sit_dmy1].ri[1]-p[j].ri[1]); 
				Rn0x = pbcx(p[sit_dmy1].R[0]-p[j].R[0]); Rn0y = pbcy(p[sit_dmy1].R[1]-p[j].R[1] );
				if(((j+yi)%2==0)&(((xi +N- 1)%N + ((yi + N)%N)*N)==sit_dmy1))
				{continue;}
				chi += pow((rn0x - (1.0 +eps11)*(Rn0x) -(eps12)*(Rn0y)),2) +  pow((rn0y-(1 +eps22)*(Rn0y) -(eps21)*(Rn0x)),2);
			}
			p[j].chi = chi;
			chi_total+= chi;
      }
return chi_total;
}
