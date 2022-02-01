#include"global.h"



double d_hchi_rectangle(int moved_particle,double *jstdelChi)
{
	int nbr_sq =8;
	//-------------------------------------Define variables---------------------------------------------------------//
	int i=0,j=0,k=0,xi=0,yi=0,N=0,sit=0,sit_dmy1=0;
	double change_in_bias=0,rn0x=0.0,rn0y=0.0,Rn0x=0.0,Rn0y=0.0;

	double detY_new=0.0,detY_inv_new=0.0,X11=0.0,X12=0.0,X21=0.0,X22=0.0,Y11=0.0,Y12=0.0,Y21=0.0,Y22=0.0,eps11=0.0,eps12=0.0,eps21=0.0,eps22=0.0;
	double chi_new=0.0,chi_total_new=0.0;//X_new[2][2]={0.0},Y_new[2][2]={0.0},eps_new[2][2]={0.0}
	
	double detY_old=0.0,detY_inv_old=0.0;
	double chi_old=0.0,chi_total_old=0.0;//X_old[2][2]={0.0},Y_old[2][2]={0.0},eps_old[2][2]={0.0},
	N = sqrt(totalN);
	//---------------------------------------------------------------------------------------------------------------//
		
	for(j=0;j<nbr_sq+1;j++) //where particle 0 consists the moved_particle id
	{
	
	//**************************************CHI-NEW************************************************************//
		//p[sit].chi_temp = 0.0;
      chi_new = 0;X11=0.0;X12=0.0;X21=0.0;X22=0.0;Y11=0.0;Y12=0.0;Y21=0.0;Y22=0.0;
      
		sit = p[moved_particle].nn_static[j];
	//-------------getting precalculated Y matrix---------//
		/*Y11 = p[sit].Ymat[0];
		Y12 = p[sit].Ymat[1];
		Y21 = p[sit].Ymat[2];
		Y22 = p[sit].Ymat[3];*/
		Y11 = Ymat[0];
		Y12 = Ymat[1];
		Y21 = Ymat[2];
		Y22 = Ymat[3];
	//-----------------------------------------------------//
		for(k=1;k<nbr_sq+1;k++)
		{	
			sit_dmy1 = p[sit].nn_static[k];
			rn0x = pbcx( p[sit_dmy1].ri[0] - p[sit].ri[0]); rn0y = pbcy(p[sit_dmy1].ri[1] - p[sit].ri[1]); 
			Rn0x = pbcx( p[sit_dmy1].R[0] -  p[sit].R[0] ); Rn0y = pbcy(p[sit_dmy1].R[1]  - p[sit].R[1] );
			
			X11 += (rn0x)*(Rn0x);
			X12 += (rn0x)*(Rn0y);
			X21 += (rn0y)*(Rn0x);
			X22 += (rn0y)*(Rn0y);
		} 
	//---------------------------------------------------//
		
		detY_new = (Y11*Y22) - (Y12*Y21) ; detY_inv_new = 1.0/detY_new;

		eps11	=	(-X12*Y21 + X11*Y22)*detY_inv_new - 1;
		eps12	=	( X12*Y11 - X11*Y12)*detY_inv_new;
		eps21	=	(-X22*Y21 + X21*Y22)*detY_inv_new;
		eps22	=	( X22*Y11 - X21*Y12)*detY_inv_new - 1;
	//-----------this for loop calc the chi_sit--------//
		for(k=1;k<nbr_sq+1;k++)
		{ 
			sit_dmy1 = p[sit].nn_static[k];
			rn0x = pbcx( p[sit_dmy1].ri[0] - p[sit].ri[0]); rn0y = pbcy(p[sit_dmy1].ri[1] - p[sit].ri[1]); 
			Rn0x = pbcx( p[sit_dmy1].R[0] -  p[sit].R[0] ); Rn0y = pbcy(p[sit_dmy1].R[1]  - p[sit].R[1] );
			
			chi_new += pow((rn0x-(1.0+eps11)*(Rn0x)-(eps12)*(Rn0y)),2)+pow((rn0y-(1 +eps22)*(Rn0y)-(eps21)*(Rn0x)),2);		
		}
		p[sit].chi_temp = chi_new;
		chi_total_new+= chi_new;
		chi_total_old += p[p[moved_particle].nn_static[j]].chi;

	//*************************************change in bias term *****************************************************************//

		//change_in_bias += h_field(PBCr(p[moved_particle].R[0], p[moved_particle].R[1]))*((p[p[moved_particle].nn_static[j]].chi_temp - p[p[moved_particle].nn_static[j]].chi));
		change_in_bias += h_field(PBCr(p[sit].R[0], p[sit].R[1]))*((p[sit].chi_temp - p[sit].chi));
		//printf("%lf\t%lf\n",p[p[moved_particle].nn_static[j]].chi_temp,p[p[moved_particle].nn_static[j]].chi );
	}
	*jstdelChi = chi_total_new - chi_total_old;
	//printf("%lf\n",chi_total_new-chi_total_old);
return change_in_bias;
}

////////////////////////////////
void update_chi__rectangle(int mvd_part)
{
	int sit=0,j=0;
	int nbr_sq =8;
	for(j =0;j<nbr_sq+1;j++)
	{
		sit = p[mvd_part].nn_static[j];
		p[sit].chi = p[sit].chi_temp;
	}

}