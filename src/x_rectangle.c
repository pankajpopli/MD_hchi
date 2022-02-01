#include"global.h"


double X_rectangle()
{

	int i=0,j=0,k=0,xi=0,yi=0,neighbours=8,N=0,sit_dmy1=0;
	double detY=0,detY_inv=0,rn0x=0.0,rn0y=0.0,Rn0x=0.0,Rn0y=0.0;
	double chi=0.0,chi_total=0.0,X11=0.0,X12=0.0,X21=0.0,X22=0.0,Y11=0.0,Y12=0.0,Y21=0.0,Y22=0.0,eps11=0.0,eps12=0.0,eps21=0.0,eps22=0.0;
	N = sqrt(totalN);

//-------------------------------
	
		for(j=0;j<totalN;j++) 
		{
			X11 = X21 = X12 = X22 = 0;
         Y11 = Y21 = Y12 = Y22 = 0;
         chi = 0;
			for(k=1;k<neighbours+1;k++)
			{
				sit_dmy1 = p[j].nn_static[k];
				rn0x = pbcx(p[sit_dmy1].ri[0] - p[j].ri[0]); rn0y = pbcy(p[sit_dmy1].ri[1] - p[j].ri[1]);
				Rn0x = pbcx(p[sit_dmy1].R[0] - p[j].R[0]); Rn0y = pbcy(p[sit_dmy1].R[1] - p[j].R[1] );
				
				X11 += (rn0x)*(Rn0x);
				X12 += (rn0x)*(Rn0y);
				X21 += (rn0y)*(Rn0x);
				X22 += (rn0y)*(Rn0y);
				
				Y11 += (Rn0x)*(Rn0x);
				Y12 += (Rn0x)*(Rn0y);
				Y21 += (Rn0y)*(Rn0x);
				Y22 += (Rn0y)*(Rn0y);
			}
			
			detY = (Y11*Y22) - (Y12*Y21) ; detY_inv = 1.0/detY;
			
			eps11	=	(-X12*Y21 + X11*Y22 )*detY_inv - 1;
			eps12	=	( X12*Y11 - X11*Y12 )*detY_inv;
			eps21	=	(-X22*Y21 + X21*Y22 )*detY_inv;
			eps22	=	( X22*Y11 - X21*Y12 )*detY_inv - 1;
			
			for(k=1;k<neighbours+1;k++)
			{
				sit_dmy1= p[j].nn_static[k];
				rn0x = pbcx(p[sit_dmy1].ri[0]-p[j].ri[0]); rn0y = pbcy(p[sit_dmy1].ri[1]-p[j].ri[1]); 
				Rn0x = pbcx(p[sit_dmy1].R[0]-p[j].R[0]); Rn0y = pbcy(p[sit_dmy1].R[1]-p[j].R[1] );
				
				chi += pow((rn0x-(1.0+eps11)*(Rn0x)-(eps12)*(Rn0y)),2)+pow((rn0y-(1 +eps22)*(Rn0y)-(eps21)*(Rn0x)),2);
			}
			p[j].chi = chi;
			chi_total+= chi;
		}
return chi_total;
}