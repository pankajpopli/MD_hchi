//------------------ potentail  polar form -------------------//
#include "global.h"

/////////// local function for pot_array_pointer ///////////////////////////////////////////////////////////////////////////
double pot_GCM(double distance, int mvdpart);
double pot_lj(double distance, int mvdpart);
double pot_r12(double distance,int mvdpart);
double pot_sq_ntwrk(double distance, int mvdpart);
double pot_sq_ntwrk_plus_bending(double distance, int mvdpart);
double pot_warning(double distance, int mvdpart);
double pot_no_pot(double distance, int mvdpart);

double pot_triang_ntwrk(double distance, int mvdpart);
double pot_triang_ntwrk_plus_bending(double distance, int mvdpart);
double pot_honeycomb_ntwrk(double distance, int mvdpart);
double pot_honeycomb_ntwrk_plus_bending(double distance, int mvdpart);
double pot_kagome_ntwrk(double distance, int mvdpart);

double changeE_sq_network(int mvdpart);

double changeE_triang_network(int mvdpart);
double changeE_triang_network_plus_bending(int mvdpart);

double changeE_honeycomb_network(int mvdpart);
double changeE_honeycomb_network_plus_bending(int mvdpart);

double changeE_kagome_network(int mvdpart);



////////////////////////////////////////// pointer array definition ///////////////////////////////////////////////////////
double (*pot_ptr_arr[])(double distance, int mvdpartIndex)  = {pot_warning, pot_no_pot, pot_triang_ntwrk,pot_triang_ntwrk_plus_bending, pot_sq_ntwrk, pot_honeycomb_ntwrk,pot_honeycomb_ntwrk_plus_bending, pot_kagome_ntwrk, pot_r12, pot_lj, pot_GCM};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





//#######################################################################################################################
//####################### function to be called from main ###############################################################
//#######################################################################################################################
//[M1] functions to be called from main @@@@
double potential(double distance, int mvdpartIndex)
{	
	return (*pot_ptr_arr[pot_slctr])(distance,mvdpartIndex);
}
//#######################################################################################################################





//#######################################################################################################################
//################################### define potentials below ##########################################################
//#######################################################################################################################
double pot_GCM(double distance, int mvdpart)
{
	return epsilon*exp(-pow((distance/sigma),2));
}
double pot_lj(double distance, int mvdpart)
{
	return 4.0*epsilon*( pow((sigma/distance),12) - pow((sigma/distance),6) );
}
double pot_r12(double distance,int mvdpart)
{
	return epsilon*(pow((sigma/distance),12));
}
double pot_sq_ntwrk(double distance, int mvdpart)
{
	return changeE_sq_network(mvdpart);
}
double pot_no_pot(double distance, int mvdpart)
{
   return 0.0;
}
double pot_warning(double distance, int mvdpart)
{
	printf(RED"\n Error ! Wrong choice of potential in Global.c\n"RESET);
	exit(0);
	return -1;
}
double pot_triang_ntwrk(double distance, int mvdpart)
{
	return changeE_triang_network(mvdpart) ;
}
double pot_triang_ntwrk_plus_bending(double distance, int mvdpart)
{
	return changeE_triang_network_plus_bending(mvdpart) ;
}
double pot_honeycomb_ntwrk(double distance, int mvdpart)
{	
	return changeE_honeycomb_network(mvdpart) ;
}
double pot_honeycomb_ntwrk_plus_bending(double distance, int mvdpart)
{	
	return changeE_honeycomb_network_plus_bending(mvdpart) ;
}
double pot_kagome_ntwrk(double distance, int mvdpart)
{	
	return changeE_kagome_network(mvdpart) ;
}
//#######################################################################################################################




/// [2] change in energy in sq lattice_network
double inline changeE_sq_network(int mvdpart) /// only for square
{
	int bondno=0,Dbondno=0,dirctn=0;
	double E_old=0.0,E_new=0.0;
	
   // energy before displacment of index particle
	for(dirctn = 1;dirctn<5;dirctn++){
		bondno =	p[mvdpart].nn_static[dirctn];
		E_old +=  (0.5*k1* pow((PBCr( (p[bondno].ri0[0] - p[mvdpart].ri0[0]),  (p[bondno].ri0[1] - p[mvdpart].ri0[1]) )-a),2));
		E_new += (0.5*k1* pow((PBCr( (p[bondno].ri0[0] - p[mvdpart].ri[0]),  (p[bondno].ri0[1] - p[mvdpart].ri[1]) )-a),2));
	}
	for(dirctn = 5;dirctn<9;dirctn++){
		Dbondno =	p[mvdpart].nn_static[dirctn];
		E_old +=  (0.5*k2* pow((PBCr( (p[Dbondno].ri0[0] - p[mvdpart].ri0[0]),  (p[Dbondno].ri0[1] - p[mvdpart].ri0[1]) )-(sqrt(2.0)*a)),2));
		E_new +=  (0.5*k2* pow((PBCr( (p[Dbondno].ri0[0] - p[mvdpart].ri[0]),  (p[Dbondno].ri0[1] - p[mvdpart].ri[1]) )-(sqrt(2.0)*a)),2));
	}
return (E_new-E_old);
}


// [3] change in energy in triangular lattice_network
double inline changeE_triang_network(int mvdpart) 
{
	int bondno=0,dirctn=0;
	double E_old=0.0,E_new=0.0;

	// energy before displacment of index particle
	for(dirctn = 1;dirctn<7;dirctn++){
		bondno =	p[mvdpart].nn_static[dirctn];
		E_old += (0.5*k1* pow((PBCr( (p[bondno].ri0[0] - p[mvdpart].ri0[0]),  (p[bondno].ri0[1] - p[mvdpart].ri0[1]) )-a),2));
		E_new += (0.5*k1* pow((PBCr( (p[bondno].ri0[0] - p[mvdpart].ri[0]),  (p[bondno].ri0[1] - p[mvdpart].ri[1]) )-a),2));
	}

return (E_new-E_old);
}
double inline changeE_triang_network_plus_bending(int mvdpart){
		
	int bondno=0,Dbondno=0,dirctn=0;
	double E_old=0.0,E_new=0.0;

	double E_old_angl=0.0,E_new_angl=0.0;
	int bondno1=0,xid=0,yid=0;
	double rikx0=0.0,rjkx0=0.0,riky0=0.0,rjky0=0.0,rikx=0.0,rjkx=0.0,riky=0.0,rjky=0.0,co=0.0,cn=0.0,thetao=0.0,thetan=0.0;
	int mvdpart2=0;

	// energy before displacment of index particle
	for(dirctn = 1;dirctn<7;dirctn++){
		bondno =	p[mvdpart].nn_static[dirctn];
		E_old += (0.5*k1* pow((PBCr( (p[bondno].ri0[0] - p[mvdpart].ri0[0]),  (p[bondno].ri0[1] - p[mvdpart].ri0[1]) )-a),2));
		E_new += (0.5*k1* pow((PBCr( (p[bondno].ri0[0] - p[mvdpart].ri[0]),  (p[bondno].ri0[1] - p[mvdpart].ri[1]) )-a),2));
	}
	///------------------bending potential----------------------------------//
	
	for(dirctn = 1;dirctn<6;dirctn++)
	{
		bondno =	p[mvdpart].nn_static[dirctn];
		bondno1 =	p[mvdpart].nn_static[dirctn+1];
		rikx0=pbcx(p[bondno].ri0[0]-p[mvdpart].ri0[0]);
		rjkx0=pbcx(p[bondno1].ri0[0]-p[mvdpart].ri0[0]);
		riky0=pbcy(p[bondno].ri0[1]-p[mvdpart].ri0[1]);
		rjky0=pbcy(p[bondno1].ri0[1]-p[mvdpart].ri0[1]);
		
		rikx=pbcx(p[bondno].ri0[0]-p[mvdpart].ri[0]);
		rjkx=pbcx(p[bondno1].ri0[0]-p[mvdpart].ri[0]);
		riky=pbcy(p[bondno].ri0[1]-p[mvdpart].ri[1]);
		rjky=pbcy(p[bondno1].ri0[1]-p[mvdpart].ri[1]);
		
		co = (rikx0*rjkx0 + riky0*rjky0)/(sqrt(pow(rikx0,2)+pow(riky0,2))*sqrt(pow(rjkx0,2)+pow(rjky0,2)));
		cn = (rikx*rjkx + riky*rjky)/(sqrt(pow(rikx,2)+pow(riky,2))*sqrt(pow(rjkx,2)+pow(rjky,2)));
		if (co<-1) co=-1;if (co>1) co=1;
		if (cn<-1) cn=-1;if (cn>1) cn=1;
		
		thetao = acos(co);thetan = acos(cn);
		
		E_old_angl += 0.5*k2*pow((thetao-(60.0*PI/180.0)),2);
		E_new_angl += 0.5*k2*pow((thetan-(60.0*PI/180.0)),2);
	}
	
		bondno =	p[mvdpart].nn_static[6];
		bondno1 =	p[mvdpart].nn_static[1];
		rikx0=pbcx(p[bondno].ri0[0]-p[mvdpart].ri0[0]);
		rjkx0=pbcx(p[bondno1].ri0[0]-p[mvdpart].ri0[0]);
		riky0=pbcy(p[bondno].ri0[1]-p[mvdpart].ri0[1]);
		rjky0=pbcy(p[bondno1].ri0[1]-p[mvdpart].ri0[1]);
		
		rikx=pbcx(p[bondno].ri0[0]-p[mvdpart].ri[0]);
		rjkx=pbcx(p[bondno1].ri0[0]-p[mvdpart].ri[0]);
		riky=pbcy(p[bondno].ri0[1]-p[mvdpart].ri[1]);
		rjky=pbcy(p[bondno1].ri0[1]-p[mvdpart].ri[1]);

		co = (rikx0*rjkx0 + riky0*rjky0)/(sqrt(pow(rikx0,2)+pow(riky0,2))*sqrt(pow(rjkx0,2)+pow(rjky0,2)));
		cn = (rikx*rjkx + riky*rjky)/(sqrt(pow(rikx,2)+pow(riky,2))*sqrt(pow(rjkx,2)+pow(rjky,2)));
		if (co<-1) co=-1;if (co>1) co=1;
		if (cn<-1) cn=-1;if (cn>1) cn=1;
	
		thetao = acos(co);thetan = acos(cn);
		
		E_old_angl += 0.5*k2*pow((thetao-(60.0*PI/180.0)),2);
		E_new_angl += 0.5*k2*pow((thetan-(60.0*PI/180.0)),2);
		
		for(dirctn=1;dirctn<6;dirctn++)
		{
			
			bondno =	p[mvdpart].nn_static[dirctn];
			bondno1 =	p[mvdpart].nn_static[dirctn+1];
			rikx0 = 	pbcx(p[mvdpart].ri0[0]-p[bondno].ri0[0]);
			rjkx0 = 	pbcx(p[bondno1].ri0[0]-p[bondno].ri0[0]);
			riky0 = 	pbcy(p[mvdpart].ri0[1]-p[bondno].ri0[1]);
			rjky0 = 	pbcy(p[bondno1].ri0[1]-p[bondno].ri0[1]);
			
			rikx = 	pbcx(p[mvdpart].ri[0]-p[bondno].ri[0]);
			rjkx = 	pbcx(p[bondno1].ri[0]-p[bondno].ri[0]);
			riky = 	pbcy(p[mvdpart].ri[1]-p[bondno].ri[1]);
			rjky = 	pbcy(p[bondno1].ri[1]-p[bondno].ri[1]);
			
			co = ((rikx0*rjkx0 + riky0*rjky0)/(sqrt(pow(rikx0,2) + pow(riky0,2))*sqrt(pow(rjkx0,2) + pow(rjky0,2))));
			cn = ((rikx*rjkx + riky*rjky)/(sqrt(pow(rikx,2) + pow(riky,2))*sqrt(pow(rjkx,2) + pow(rjky,2))));
			if (co<-1) co=-1;if (co>1) co=1;
			if (cn<-1) cn=-1;if (cn>1) cn=1;
			
			thetao = acos(co);thetan = acos(cn);
			
			E_old_angl += 0.5*k2*pow((thetao-(60.0*PI/180.0)),2);
			E_new_angl += 0.5*k2*pow((thetan-(60.0*PI/180.0)),2);
			
			bondno =	p[mvdpart].nn_static[7-dirctn];
			bondno1 =	p[mvdpart].nn_static[7-dirctn-1];
			rikx0 = 	pbcx(p[mvdpart].ri0[0]-p[bondno].ri0[0]);
			rjkx0 = 	pbcx(p[bondno1].ri0[0]-p[bondno].ri0[0]);
			riky0 = 	pbcy(p[mvdpart].ri0[1]-p[bondno].ri0[1]);
			rjky0 = 	pbcy(p[bondno1].ri0[1]-p[bondno].ri0[1]);
			
			rikx = 	pbcx(p[mvdpart].ri[0]-p[bondno].ri[0]);
			rjkx = 	pbcx(p[bondno1].ri[0]-p[bondno].ri[0]);
			riky = 	pbcy(p[mvdpart].ri[1]-p[bondno].ri[1]);
			rjky = 	pbcy(p[bondno1].ri[1]-p[bondno].ri[1]);
			
			co = ((rikx0*rjkx0 + riky0*rjky0)/(sqrt(pow(rikx0,2) + pow(riky0,2))*sqrt(pow(rjkx0,2) + pow(rjky0,2))));
			cn = ((rikx*rjkx + riky*rjky)/(sqrt(pow(rikx,2) + pow(riky,2))*sqrt(pow(rjkx,2) + pow(rjky,2))));
			if (co<-1) co=-1;if (co>1) co=1;
			if (cn<-1) cn=-1;if (cn>1) cn=1;
			
			thetao = acos(co);thetan = acos(cn);
			
			E_old_angl += 0.5*k2*pow((thetao-(60.0*PI/180.0)),2);
			E_new_angl += 0.5*k2*pow((thetan-(60.0*PI/180.0)),2);
		}
		
			bondno =	p[mvdpart].nn_static[1];
			bondno1 =	p[mvdpart].nn_static[6];
			rikx0 = 	pbcx(p[mvdpart].ri0[0]-p[bondno].ri0[0]);
			rjkx0 = 	pbcx(p[bondno1].ri0[0]-p[bondno].ri0[0]);
			riky0 = 	pbcy(p[mvdpart].ri0[1]-p[bondno].ri0[1]);
			rjky0 = 	pbcy(p[bondno1].ri0[1]-p[bondno].ri0[1]);
			
			rikx = 	pbcx(p[mvdpart].ri[0]-p[bondno].ri[0]);
			rjkx = 	pbcx(p[bondno1].ri[0]-p[bondno].ri[0]);
			riky = 	pbcy(p[mvdpart].ri[1]-p[bondno].ri[1]);
			rjky = 	pbcy(p[bondno1].ri[1]-p[bondno].ri[1]);
			
			co = ((rikx0*rjkx0 + riky0*rjky0)/(sqrt(pow(rikx0,2) + pow(riky0,2))*sqrt(pow(rjkx0,2) + pow(rjky0,2))));
			cn = ((rikx*rjkx + riky*rjky)/(sqrt(pow(rikx,2) + pow(riky,2))*sqrt(pow(rjkx,2) + pow(rjky,2))));
			if (co<-1) co=-1;if (co>1) co=1;
			if (cn<-1) cn=-1;if (cn>1) cn=1;
			
			thetao = acos(co);thetan = acos(cn);
			
			E_old_angl += 0.5*k2*pow((thetao-(60.0*PI/180.0)),2);
			E_new_angl += 0.5*k2*pow((thetan-(60.0*PI/180.0)),2);
			
			
			
			bondno =	p[mvdpart].nn_static[6];
			bondno1 =	p[mvdpart].nn_static[1];
			rikx0 = 	pbcx(p[mvdpart].ri0[0]-p[bondno].ri0[0]);
			rjkx0 = 	pbcx(p[bondno1].ri0[0]-p[bondno].ri0[0]);
			riky0 = 	pbcy(p[mvdpart].ri0[1]-p[bondno].ri0[1]);
			rjky0 = 	pbcy(p[bondno1].ri0[1]-p[bondno].ri0[1]);
			
			rikx = 	pbcx(p[mvdpart].ri[0]-p[bondno].ri[0]);
			rjkx = 	pbcx(p[bondno1].ri[0]-p[bondno].ri[0]);
			riky = 	pbcy(p[mvdpart].ri[1]-p[bondno].ri[1]);
			rjky = 	pbcy(p[bondno1].ri[1]-p[bondno].ri[1]);
			
			co = ((rikx0*rjkx0 + riky0*rjky0)/(sqrt(pow(rikx0,2) + pow(riky0,2))*sqrt(pow(rjkx0,2) + pow(rjky0,2))));
			cn = ((rikx*rjkx + riky*rjky)/(sqrt(pow(rikx,2) + pow(riky,2))*sqrt(pow(rjkx,2) + pow(rjky,2))));
			if (co<-1) co=-1;if (co>1) co=1;
			if (cn<-1) cn=-1;if (cn>1) cn=1;
			
			thetao = acos(co);thetan = acos(cn);
			
			E_old_angl += 0.5*k2*pow((thetao-(60.0*PI/180.0)),2);
			E_new_angl += 0.5*k2*pow((thetan-(60.0*PI/180.0)),2);

return (E_new-E_old)+(E_new_angl-E_old_angl);
}

// [4] change in energy in honeycomb lattice_network
double inline changeE_honeycomb_network(int mvdpart) 
{
	int bondno=0,Dbondno=0,dirctn=0;
	double E_old=0.0,E_new=0.0;

	double E_old_angl=0.0,E_new_angl=0.0;
	int bondno1=0,xid,yid;
	double rikx0=0.0,rjkx0=0.0,riky0=0.0,rjky0=0.0,rikx=0.0,rjkx=0.0,riky=0.0,rjky=0.0,co=0.0,cn=0.0,thetao=0.0,thetan=0.0;
	int mvdpart2=0;
	// energy before displacment of index particle
	for(dirctn = 1;dirctn<4;dirctn++)
	{
		//bondno =	p[mvdpart].nbr_lst[dirctn];
		bondno =	p[mvdpart].nn_static[dirctn];
		E_old += (0.5*k1*pow((PBCr( (p[bondno].ri0[0] - p[mvdpart].ri0[0]),  (p[bondno].ri0[1] - p[mvdpart].ri0[1]) )-a),2));
		E_new += (0.5*k1* pow((PBCr( (p[bondno].ri0[0] - p[mvdpart].ri[0]),  (p[bondno].ri0[1] - p[mvdpart].ri[1]) )-a),2));
		
	}
	for(dirctn = 4;dirctn<10;dirctn++)
	{
		//Dbondno =	p[mvdpart].nbr_lst[dirctn];
		Dbondno =	p[mvdpart].nn_static[dirctn];
		E_old += (0.5*k2* pow((PBCr( (p[Dbondno].ri0[0] - p[mvdpart].ri0[0]),  (p[Dbondno].ri0[1] - p[mvdpart].ri0[1]) )-(sqrt(3.0)*a)),2));
		E_new += (0.5*k2* pow((PBCr( (p[Dbondno].ri0[0] - p[mvdpart].ri[0]),  (p[Dbondno].ri0[1] - p[mvdpart].ri[1]) )-(sqrt(3.0)*a)),2));
		
	}

return (E_new-E_old);
}
double inline changeE_honeycomb_network_plus_bending(int mvdpart){
//------------------------------------------bending energy change -----------------------------------------------------------//
	
	int bondno=0,Dbondno=0,dirctn=0;
	double E_old=0.0,E_new=0.0;

	double E_old_angl=0.0,E_new_angl=0.0;
	int bondno1=0,xid,yid;
	double rikx0=0.0,rjkx0=0.0,riky0=0.0,rjky0=0.0,rikx=0.0,rjkx=0.0,riky=0.0,rjky=0.0,co=0.0,cn=0.0,thetao=0.0,thetan=0.0;
	int mvdpart2=0;
	// energy before displacment of index particle
	for(dirctn = 1;dirctn<4;dirctn++)
	{
		//bondno =	p[mvdpart].nbr_lst[dirctn];
		bondno =	p[mvdpart].nn_static[dirctn];
		E_old += (0.5*k1*pow((PBCr( (p[bondno].ri0[0] - p[mvdpart].ri0[0]),  (p[bondno].ri0[1] - p[mvdpart].ri0[1]) )-a),2));
		E_new += (0.5*k1* pow((PBCr( (p[bondno].ri0[0] - p[mvdpart].ri[0]),  (p[bondno].ri0[1] - p[mvdpart].ri[1]) )-a),2));
		
	}
	for(dirctn = 4;dirctn<10;dirctn++)
	{
		//Dbondno =	p[mvdpart].nbr_lst[dirctn];
		Dbondno =	p[mvdpart].nn_static[dirctn];
		E_old += (0.5*k2* pow((PBCr( (p[Dbondno].ri0[0] - p[mvdpart].ri0[0]),  (p[Dbondno].ri0[1] - p[mvdpart].ri0[1]) )-(sqrt(3.0)*a)),2));
		E_new += (0.5*k2* pow((PBCr( (p[Dbondno].ri0[0] - p[mvdpart].ri[0]),  (p[Dbondno].ri0[1] - p[mvdpart].ri[1]) )-(sqrt(3.0)*a)),2));
		
	}
	
	int bond_i,bond_j,bond_k;
	double kb=1.0;
	
	xid = mvdpart%(int)(sqrt(totalN));
	yid = floor(mvdpart/(int)(sqrt(totalN)));


	if((mvdpart+yid)%2==1 )	
	{	
		int bondi_list[4]={3,1,0,0};
		int bondj_list[4]={2,2,0,0};
		for(dirctn = 0;dirctn<2;dirctn++) //first three bonds
		{
			bond_i =	p[mvdpart].nn_static[bondi_list[dirctn]];
			bond_j =	p[mvdpart].nn_static[bondj_list[dirctn]];
			bond_k = 	mvdpart;
			//printf("%d\t%d\t%d\n",bond_i,bond_k,bond_j);
			rikx0=pbcx(p[bond_i].ri0[0]-p[bond_k].ri0[0]);
			rjkx0=pbcx(p[bond_j].ri0[0]-p[bond_k].ri0[0]);
			riky0=pbcy(p[bond_i].ri0[1]-p[bond_k].ri0[1]);
			rjky0=pbcy(p[bond_j].ri0[1]-p[bond_k].ri0[1]);
			
			rikx=pbcx(p[bond_i].ri0[0]-p[bond_k].ri[0]);
			rjkx=pbcx(p[bond_j].ri0[0]-p[bond_k].ri[0]);
			riky=pbcy(p[bond_i].ri0[1]-p[bond_k].ri[1]);
			rjky=pbcy(p[bond_j].ri0[1]-p[bond_k].ri[1]);
			
			co = (rikx0*rjkx0 + riky0*rjky0)/(sqrt(pow(rikx0,2)+pow(riky0,2))*sqrt(pow(rjkx0,2)+pow(rjky0,2)));
			cn = (rikx*rjkx + riky*rjky)/(sqrt(pow(rikx,2)+pow(riky,2))*sqrt(pow(rjkx,2)+pow(rjky,2)));
			if (co<-1) co=-1;if (co>1) co=1;
			if (cn<-1) cn=-1;if (cn>1) cn=1;
			
			thetao = acos(co);thetan = acos(cn);
			
			E_old_angl += 0.5*kb*pow((thetao-(120.0*PI/180.0)),2);
			E_new_angl += 0.5*kb*pow((thetan-(120.0*PI/180.0)),2);
		}
		  bondi_list[0]=6;bondi_list[1]=7;bondi_list[2]=4;bondi_list[3]=5;
		 int bondk_list[4]={3,1,2,2};
		for(dirctn = 0;dirctn<4;dirctn++)	// last 6 bonds
		{
			bond_i =	p[mvdpart].nn_static[bondi_list[dirctn]];
			bond_j =	mvdpart;
			bond_k = 	p[mvdpart].nn_static[bondk_list[dirctn]];
			//printf("%d\t%d\t%d\n",bond_i,bond_k,bond_j);
			rikx0=pbcx(p[bond_i].ri0[0]-p[bond_k].ri0[0]);
			rjkx0=pbcx(p[bond_j].ri0[0]-p[bond_k].ri0[0]);
			riky0=pbcy(p[bond_i].ri0[1]-p[bond_k].ri0[1]);
			rjky0=pbcy(p[bond_j].ri0[1]-p[bond_k].ri0[1]);
			
			rikx=pbcx(p[bond_i].ri0[0]-p[bond_k].ri0[0]);
			rjkx=pbcx(p[bond_j].ri[0]-p[bond_k].ri0[0]);
			riky=pbcy(p[bond_i].ri0[1]-p[bond_k].ri0[1]);
			rjky=pbcy(p[bond_j].ri[1]-p[bond_k].ri0[1]);
			
			co = (rikx0*rjkx0 + riky0*rjky0)/(sqrt(pow(rikx0,2)+pow(riky0,2))*sqrt(pow(rjkx0,2)+pow(rjky0,2)));
			cn = (rikx*rjkx + riky*rjky)/(sqrt(pow(rikx,2)+pow(riky,2))*sqrt(pow(rjkx,2)+pow(rjky,2)));
			if (co<-1) co=-1;if (co>1) co=1;
			if (cn<-1) cn=-1;if (cn>1) cn=1;
			
			thetao = acos(co);thetan = acos(cn);
			
			E_old_angl += 0.5*kb*pow((thetao-(120.0*PI/180.0)),2);
			E_new_angl += 0.5*kb*pow((thetan-(120.0*PI/180.0)),2);
		}

	}else if((mvdpart+yid)%2==0 )	
	{	
		int bondi_list[4]={2,1,0,0};
		int bondj_list[4]={3,3,0,0};
		for(dirctn = 0;dirctn<2;dirctn++) //first three bonds
		{
			bond_i =	p[mvdpart].nn_static[bondi_list[dirctn]];
			bond_j =	p[mvdpart].nn_static[bondj_list[dirctn]];
			bond_k = 	mvdpart;
			//printf("%d\t%d\t%d\n",bond_i,bond_k,bond_j);
			rikx0=pbcx(p[bond_i].ri0[0]-p[bond_k].ri0[0]);
			rjkx0=pbcx(p[bond_j].ri0[0]-p[bond_k].ri0[0]);
			riky0=pbcy(p[bond_i].ri0[1]-p[bond_k].ri0[1]);
			rjky0=pbcy(p[bond_j].ri0[1]-p[bond_k].ri0[1]);
			
			rikx=pbcx(p[bond_i].ri0[0]-p[bond_k].ri[0]);
			rjkx=pbcx(p[bond_j].ri0[0]-p[bond_k].ri[0]);
			riky=pbcy(p[bond_i].ri0[1]-p[bond_k].ri[1]);
			rjky=pbcy(p[bond_j].ri0[1]-p[bond_k].ri[1]);
			
			co = (rikx0*rjkx0 + riky0*rjky0)/(sqrt(pow(rikx0,2)+pow(riky0,2))*sqrt(pow(rjkx0,2)+pow(rjky0,2)));
			cn = (rikx*rjkx + riky*rjky)/(sqrt(pow(rikx,2)+pow(riky,2))*sqrt(pow(rjkx,2)+pow(rjky,2)));
			if (co<-1) co=-1;if (co>1) co=1;
			if (cn<-1) cn=-1;if (cn>1) cn=1;
			
			thetao = acos(co);thetan = acos(cn);
			
			E_old_angl += 0.5*kb*pow((thetao-(120.0*PI/180.0)),2);
			E_new_angl += 0.5*kb*pow((thetan-(120.0*PI/180.0)),2);
		}
		 bondi_list[0]=4;bondi_list[1]=5;bondi_list[2]=6;bondi_list[3]=7;
		 int bondk_list[4]={2,1,3,3};
		for(dirctn = 0;dirctn<4;dirctn++)	// last 6 bonds
		{
			bond_i =	p[mvdpart].nn_static[bondi_list[dirctn]];
			bond_j =	mvdpart;
			bond_k = 	p[mvdpart].nn_static[bondk_list[dirctn]];
			//printf("%d\t%d\t%d\n",bond_i,bond_k,bond_j);
			rikx0=pbcx(p[bond_i].ri0[0]-p[bond_k].ri0[0]);
			rjkx0=pbcx(p[bond_j].ri0[0]-p[bond_k].ri0[0]);
			riky0=pbcy(p[bond_i].ri0[1]-p[bond_k].ri0[1]);
			rjky0=pbcy(p[bond_j].ri0[1]-p[bond_k].ri0[1]);
			
			rikx=pbcx(p[bond_i].ri0[0]-p[bond_k].ri0[0]);
			rjkx=pbcx(p[bond_j].ri[0]-p[bond_k].ri0[0]);
			riky=pbcy(p[bond_i].ri0[1]-p[bond_k].ri0[1]);
			rjky=pbcy(p[bond_j].ri[1]-p[bond_k].ri0[1]);
			
			co = (rikx0*rjkx0 + riky0*rjky0)/(sqrt(pow(rikx0,2)+pow(riky0,2))*sqrt(pow(rjkx0,2)+pow(rjky0,2)));
			cn = (rikx*rjkx + riky*rjky)/(sqrt(pow(rikx,2)+pow(riky,2))*sqrt(pow(rjkx,2)+pow(rjky,2)));
			if (co<-1) co=-1;if (co>1) co=1;
			if (cn<-1) cn=-1;if (cn>1) cn=1;
			
			thetao = acos(co);thetan = acos(cn);
			
			E_old_angl += 0.5*kb*pow((thetao-(120.0*PI/180.0)),2);
			E_new_angl += 0.5*kb*pow((thetan-(120.0*PI/180.0)),2);
		}

	}
return ((E_new-E_old)+(E_new_angl-E_old_angl));
}

// [5] change in energy in kagome lattice_network
double inline changeE_kagome_network(int mvdpart) /// only for square
{
	int bondno=0,Dbondno=0,dirctn=0;
	double E_old=0.0,E_new=0.0;
	
   // energy before displacment of index particle
	for(dirctn = 1;dirctn<5;dirctn++){
		//bondno =	p[mvdpart].nbr_lst[dirctn];
		bondno =	p[mvdpart].nn_static[dirctn];
		E_old += (0.5*k1*pow((PBCr( (p[bondno].ri0[0] - p[mvdpart].ri0[0]),  (p[bondno].ri0[1] - p[mvdpart].ri0[1]) )-(a*0.5)),2));
		E_new += (0.5*k1* pow((PBCr( (p[bondno].ri0[0] - p[mvdpart].ri[0]),  (p[bondno].ri0[1] - p[mvdpart].ri[1]) )-(a*0.5)),2));
	}
	for(dirctn = 5;dirctn<9;dirctn++){
		//Dbondno =	p[mvdpart].nbr_lst[dirctn];
		Dbondno =	p[mvdpart].nn_static[dirctn];
		E_old += (0.5*k2* pow((PBCr( (p[Dbondno].ri0[0] - p[mvdpart].ri0[0]),  (p[Dbondno].ri0[1] - p[mvdpart].ri0[1]) )-(sqrt(3.0)*a*0.5)),2));
		E_new += (0.5*k2* pow((PBCr( (p[Dbondno].ri0[0] - p[mvdpart].ri[0]),  (p[Dbondno].ri0[1] - p[mvdpart].ri[1]) )-(sqrt(3.0)*a*0.5)),2));
	}
   
return (E_new-E_old);
}
//------------------------------------------------------------------------------------------------------------------------------------------------//






























/*
///[1] different potentials,
double potential(double distance, char *str2)
{
	double potE=0.0;
	if(strcmp(str2,"harmonic")==0)
	{
		potE = 0.5*k1*( pow((distance - a),2) );
	}
	else if(strcmp(str2,"harmonic-square")==0)	//---> analytical expression for this is
	{															// 0.5* (k1 + 4*(k2-k1)(cos(theta)sin(theta))^2)*(|rij| - |Rij|)^2
		if((distance <= ((sqrt(2)*a)+(a/10.0))) && (distance >= (sqrt(2)*a)-(a/10.0)))
		{	
			potE = 0.5*k2*( pow( (distance - (sqrt(2.0)*a)) ,2) );
		}
		else if ((distance <= (a+(a/10.0))) && (distance >= (a-(a/10.0))))
		{
			potE = 0.5*k1*( pow((distance - a),2) );
		}else
		{
			printf("\n Error !, flutuation is too large(>a/10) to handle in this form of potential.\n");
			printf(" for large flutuation or high Temp use pot_type 'square-network'\n" );
			exit(0);
			//potE = 0.0;
		}
	}
	else if(strcmp(str2,"GCM")==0)
	{
		potE = epsilon*exp(-pow((distance/sigma),2));
	}
	else if(strcmp(str2,"lj")==0)
	{
		potE = 4.0*epsilon*( pow((sigma/distance),12) - pow((sigma/distance),6) );
	}else
	{
		printf("Wrong choice of potential, check global.c\n"); 
		exit(0);
	}
	
return potE;
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------//
// becoz function potential has no argument of mvd particle index we require below functions for netwroks and networks/harmonic springs requires 
// NN appproximation, which cant be made in potential function also even if we choose potential function, a cutoff will make the buisness messy as harmonic 
// bond could go out cutoff region and still should be counted for energy change.
//------------------------------------------------------------------------------------------------------------------------------------------------------------//

/// [2] change in energy in sq lattice_network
double changeE_sq_network(int index) /// only for square
{
	int bondno=0,Dbondno=0,dirctn=0;
	double E_old=0.0,E_new=0.0;
	// energy before displacment of index particle
	for(dirctn = 1;dirctn<5;dirctn++)
	{
		bondno =	p[index].nn_static[dirctn];
		E_old = E_old + (0.5*k1* pow((PBCr( (p[bondno].ri0[0] - p[index].ri0[0]),  (p[bondno].ri0[1] - p[index].ri0[1]) )-a),2));
		
	}
	for(dirctn = 5;dirctn<9;dirctn++)
	{
		Dbondno =	p[index].nn_static[dirctn];
		E_old = E_old + (0.5*k2* pow((PBCr( (p[Dbondno].ri0[0] - p[index].ri0[0]),  (p[Dbondno].ri0[1] - p[index].ri0[1]) )-(sqrt(2.0)*a)),2));
		
	}
	// energy after displacment of index particle
	for(dirctn = 1;dirctn<5;dirctn++)
	{
		bondno =	p[index].nn_static[dirctn];
		E_new = E_new+(0.5*k1* pow((PBCr( (p[bondno].ri0[0] - p[index].ri[0]),  (p[bondno].ri0[1] - p[index].ri[1]) )-a),2));
		
	}
	for(dirctn = 5;dirctn<9;dirctn++)
	{
		Dbondno =	p[index].nn_static[dirctn];
		E_new = E_new + (0.5*k2* pow((PBCr( (p[Dbondno].ri0[0] - p[index].ri[0]),  (p[Dbondno].ri0[1] - p[index].ri[1]) )-(sqrt(2.0)*a)),2));
	}
return (E_new-E_old);
}


/// [3] change in energy in triangular lattice_network
double changeE_triang_network(int index) 
{
	int bondno=0,Dbondno=0,dirctn=0;
	double E_old=0.0,E_new=0.0;
	// energy before displacment of index particle
	for(dirctn = 1;dirctn<7;dirctn++)
	{
		bondno =	p[index].nn_static[dirctn];
		E_old = E_old + (0.5*k1* pow((PBCr( (p[bondno].ri0[0] - p[index].ri0[0]),  (p[bondno].ri0[1] - p[index].ri0[1]) )-a),2));
		
	}
	// energy after displacment of index particle
	for(dirctn = 1;dirctn<7;dirctn++)
	{
		bondno =	p[index].nn_static[dirctn];
		E_new = E_new+(0.5*k1* pow((PBCr( (p[bondno].ri0[0] - p[index].ri[0]),  (p[bondno].ri0[1] - p[index].ri[1]) )-a),2));
		
	}
return (E_new-E_old);
}
//------------------------------------------------------------------------------------------------------------------------------------------------//



///[4] Non-affine field h{Ri}
double h_field(double distance)
{
	double ref_dist = distance;
	if(distance <= 5*sqrt(2)*a)
	{
		return Vchi;
		//return (hchi + Vchi*exp( -((fabsf(ref_dist))/(lambda))    ) );
	}else
	{
		//return (hchi + Vchi*exp( -((fabsf(ref_dist))/(lambda))    ) );
		return Vchi;
	}

}

*/


//--------------------------------------------------------END of Program--------------------------------------------------------------//
