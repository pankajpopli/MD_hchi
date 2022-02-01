#include "global.h"



///////////////////////// local function for c_lisst_delE pointer array /////////////////////////////////
double c_list_delE_ntwrk(int idx, int *c_updt_yn);
double c_list_delE_pot(int idx, int *c_updt_yn);
double c_list_delE_warng(int idx, int *c_updt_yn);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void   c_list_force_ntwrk();
void   c_list_force_pot();
void   c_list_force_war();
//////////////////////////////////// c_list_delE pointer array definition /////////////////////////////////
double(*clistDelE_ptr_arr[])(int idx , int *c_updt_check) = {c_list_delE_warng, c_list_delE_ntwrk,c_list_delE_pot};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void(*clist_force_ptr_arr[])() = {c_list_force_war, c_list_force_ntwrk,c_list_force_pot};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////Local functions and array///////////////////////////////////////////////////////////////////////
int cnn[9];
void sort_acc_cid(int pid2cid[][2]);
void cnbz(int id_cell);
double c_list_delE_partial(int idx, int new_cid);
_Bool cid_change_check(double cell_x, double cell_y, int nc_x, int nc_y, int idx, int* new_c);
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







//#####################################################################################################################################//
// [ 1 ] ########################################## Create cell list ##################################################################//
//#####################################################################################################################################//
void create_c_list()
{
	int pid2cid[totalN][2];
	int i,cx,cy,cid;
	
	//caculatin parameters for cell_list
	ncell_x = floor(lbox_x/cutoff);
	ncell_y = floor(lbox_y/cutoff);
	ncell	  = ncell_x*ncell_y;
	//printf("%d\n",ncell);exit(0);
	lcell_x = lbox_x/((double)(ncell_x)); 
	lcell_y = lbox_y/((double)(ncell_y));

	//calculate pid and cid
	for(i=0;i<totalN;i++){
      
		cx = floor((( (p[i].R[0]+(lbox_x*0.5))/(lcell_x) )));
		cy = floor((( (p[i].R[1]+(lbox_y*0.5))/(lcell_y) )));
		//cid = cx + ncell_x*cy;
		cid = (cx+ncell_x)%ncell_x + (((cy+ncell_y)%ncell_y)*(ncell_x));
      if(cid<0){printf(RED" cell list failed,update_cell_list_full\n"RESET);exit(0);}
		map[i].cid = cid;    //struct map contains pid and corresponding cid
		map[i].pid = i;
		pid2cid[i][0] = i;   //array pid2cid contains pid and corresponding cid
		pid2cid[i][1] = cid;
	}

   // sorting according to cid, cid is in sequence but pid is not.
   sort_acc_cid(pid2cid);

   // initializing array containg length and starting point of a cell
   // clen tells the length of each cell. cbgn tells the starting adress of each cell.
   clen = malloc((ncell)*sizeof(int));
   if (clen == NULL){
  		printf( "malloc failed\n");
  		exit(0);
		}
   cbgn = malloc((ncell)*sizeof(int));
   if (cbgn == NULL){
  		printf( "malloc failed\n");
  		exit(0);
		}
    
   //creating cell list
   for(i=0;i<ncell;i++){
   	clen[i] = 0;
   	cbgn[i] = 0;
   }
   for(i=0;i<totalN;i++){
      c_list[i] = pid2cid[i][0];    //contain pid such that cid is in increasing order/sequence.
    	cid = pid2cid[i][1];
    	clen[cid]++;
   }
   for(i=1;i<ncell;i++){
      cbgn[i] = cbgn[i-1] + clen[i-1];
   }

}

// This is function updates the full cell list if pot is chosen.
void update_c_list_full()
{
	int dummy;
   
   if(clist_ntwrk_pot_slctr==1){
      dummy=1;
   }
   else if(clist_ntwrk_pot_slctr==2){
   
      int pid2cid[totalN][2];
      //int cell_sorted[totalN];
      int i,cx,cy,cid;
 
      //calculate pid and cid
      for(i=0;i<totalN;i++){
         cx = floor((( (p[i].ri[0]+(lbox_x*0.5))/(lcell_x) )));
         cy = floor((( (p[i].ri[1]+(lbox_y*0.5))/(lcell_y) )));
         //cid = cx + ncell_x*cy;
         cid = (cx+ncell_x)%ncell_x + (((cy+ncell_y)%ncell_y)*(ncell_x));
         if(cid<0){printf(RED" cell list failed,update_cell_list_full\n"RESET);exit(0);}
         map[i].cid = cid;    //struct map contains pid and corresponding cid
         map[i].pid = i;
         pid2cid[i][0] = i;   //array pid2cid contains pid and corresponding cid
         pid2cid[i][1] = cid;
      }
      sort_acc_cid(pid2cid);

      // initializing array containg length and starting point of a cell
      // clen tells the length of each cell. cbgn tells the starting adress of each cell.

      //creating cell list
      for(i=0;i<ncell;i++){
         clen[i] = 0;
         cbgn[i] = 0;
      }
      for(i=0;i<totalN;i++){
         c_list[i] = pid2cid[i][0];    //contain pid such that cid is in increasing order/sequence.
         cid = pid2cid[i][1];
         clen[cid]++;
      }
      for(i=1;i<ncell;i++){
         cbgn[i] = cbgn[i-1] + clen[i-1];
         
      }
   }else{
      printf(RED" Invalid potential or ntwrk type chosen\n"RESET);
   }
}
//#####################################################################################################################################//






//#####################################################################################################################################//
// [ 2 ] ################################## delta E using cell_list (5sep2017) #######################################################//
//######################################## using pointer function for ntwrk and normal pot ##########################################//
//###################################################################################################################################//

//[M2] function to be called from main---------
double c_list_delE(int idx, int *c_updt_yn)
{
	double return_val=0;
	int c_updt_dummy;
	return_val = (*clistDelE_ptr_arr[clist_ntwrk_pot_slctr])(idx , &c_updt_dummy);
	*c_updt_yn = c_updt_dummy;
	return return_val;
}


//[M2R0] routine to generate warning if neither pot or ntwrk chosen @@@@
double c_list_delE_warng(int idx, int *c_updt_yn)
{
	printf(RED"neither pot or ntwrk, fuckoff\n"RESET);
   return -1;exit(0);
}

//[M2R1] routine to calc change in ntwrk energy given displaced particle @@@@
double c_list_delE_ntwrk(int idx, int *c_updt_chck)
{
	double distance =0;
	*c_updt_chck = 0;
	return potential(distance,idx);
}

//[M2R2] routine to calc change in system pot energy given displaced particle index and pot @@@@
double c_list_delE_pot(int idx, int *c_updt_chck)
{
	int cid_of_mvd_pid=0,pid_in_cell=0,pid_mvd=0,new_cc=0;
	double energyO=0.0,energyN=0.0,dis_old=0.0,dis_new=0.0,result=0.0;
	int cc,pc;
	
	{
		if( cid_change_check(lcell_x, lcell_y, ncell_x, ncell_y, idx, &new_cc)==true){
			*c_updt_chck = 1;
			result =  c_list_delE_partial( idx,  new_cc);
		}
		else{ //if( cid_change_check(lcell_x, lcell_y, ncell_x, ncell_y, index, &new_cc)==false)
			*c_updt_chck = 0;
			cid_of_mvd_pid = map[idx].cid;
			cnbz(cid_of_mvd_pid); //calculate cell nbr's index of mvd particle's cell
         
			for( cc=0;cc<9;cc++){	
				for( pc=0;pc<clen[cnn[cc]];pc++){	
               if(cnn[cc]>=ncell){printf("cell list allocation exceded in pot\n\n%d\t%d\n",ncell,cnn[cc]);exit(0);}
					pid_in_cell = c_list[cbgn[cnn[cc]]+pc]; pid_mvd = idx;
					if(pid_mvd != pid_in_cell){
                  dis_old = PBCr( (p[pid_in_cell].ri0[0]-p[pid_mvd].ri0[0]), (p[pid_in_cell].ri0[1]-p[pid_mvd].ri0[1]) );
						dis_new = PBCr( (p[pid_in_cell].ri[0]-p[pid_mvd].ri[0]), (p[pid_in_cell].ri[1]-p[pid_mvd].ri[1]) );
						//dis_new = PBCr( (p[pid_mvd].ri[0]-p[pid1].ri[0]), 	(p[pid_mv].ri[1]-p[pid1].ri[1]) 	 );
						if(dis_old < cutoff){
							energyO = energyO + potential( dis_old,idx);
						}
						if(dis_new < cutoff){
							energyN = energyN + potential( dis_new,idx);
						}
					}
				}
			}
			result =  energyN-energyO;
		}
	}
return result;
}

//[Aux M2R2.1] to be used by [R2] to calc deltaE if cell id is changed @@@@
double c_list_delE_partial(int idx, int new_cid)
{
	
	int cid_of_mvd_pid=0,pid_in_cell=0,pid_mvd=0;//new_c=0;
	double energyO=0.0,energyN=0.0,dis_old=0.0,dis_new=0.0;
	int cc,pc;
	
		cid_of_mvd_pid = map[idx].cid;
		cnbz(cid_of_mvd_pid);
		// energy if particle in original cell
		for( cc=0;cc<9;cc++){	
			for( pc=0;pc<clen[cnn[cc]];pc++){	
				if(cnn[cc]>=ncell){printf("cell list allocation exceded in partial\n\n%d\t%d\n",ncell,cnn[cc]);exit(0);}
				pid_in_cell = c_list[cbgn[cnn[cc]]+pc]; pid_mvd = idx;
				if(pid_mvd != pid_in_cell){	
					dis_old = PBCr( (p[pid_in_cell].ri0[0]-p[pid_mvd].ri0[0]), (p[pid_in_cell].ri0[1]-p[pid_mvd].ri0[1]) );
					if(dis_old < cutoff){
						energyO = energyO + potential( dis_old,idx);
					}
					
				}
			}
		}
		//energy if particle in new cell
		cnbz(new_cid);
		for( cc=0;cc<9;cc++){	
			for( pc=0;pc<clen[cnn[cc]];pc++){	
				if(cnn[cc]>=ncell){	
					//printf("%d\n",new_cid);
					printf("cell list allocation exceded in partial 2\n\n%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",ncell,cnn[cc],cc,cnn[0],cnn[1],cnn[2],cnn[3],cnn[4],cnn[5],cnn[6],cnn[7],cnn[8],ncell_x,ncell_y);
					exit(0);
				}
				pid_in_cell = c_list[cbgn[cnn[cc]]+pc]; pid_mvd = idx;
				if(pid_mvd != pid_in_cell){	
					dis_new = PBCr( (p[pid_in_cell].ri[0]-p[pid_mvd].ri[0]), (p[pid_in_cell].ri[1]-p[pid_mvd].ri[1]) );
					if(dis_new < cutoff){
						energyN = energyN + potential( dis_new,idx);
					}
					
				}
			}
		}
		
return energyN-energyO;
}

// [Aux M2R2.2] to check whther cell list needs to be updated or not @@@
_Bool cid_change_check(double cell_x, double cell_y, int nc_x, int nc_y, int idx, int* new_c)
{
	int check,new_cid,new_cx,new_cy;
	//double mvd_cord_x = p[index].ri[0];
	//double mvd_cord_y = p[index].ri[1];
	new_cx = floor((( (p[idx].ri[0]+(lbox_x*0.5))/(cell_x) )));
	new_cy = floor((( (p[idx].ri[1]+(lbox_y*0.5))/(cell_y) )));
	//new_cid = (new_cx+nc_x)%nc_x + (((new_cy+nc_x)%nc_x)*(nc_x));
	new_cid = (new_cx+nc_x)%nc_x + (((new_cy+nc_y)%nc_y)*(nc_x));
	*new_c = new_cid;
	if( new_cid != map[idx].cid){
		check = true;
	}else{
		check = false;
	}
return check;
}
//###################################################################################################################################//






//#####################################################################################################################################//
// [ 3 ] ################################## force using cell_list (28june2018) #######################################################//
//######################################## using pointer function for ntwrk and normal pot ##########################################//
//###################################################################################################################################//
//[M3] function to be called from main
void c_list_forceC()
{
   (*clist_force_ptr_arr[clist_ntwrk_pot_slctr])();
}

//[M3R0] routine to calculate force in a ntwrk for all particles @@@@
void  c_list_force_war()
{
   printf("neither pot or ntwrk, fuckoff\n");exit(0);
}

//[M3R1] routine to calculate force in a ntwrk for all particles @@@@
void  c_list_force_ntwrk()
{  
   double fx=0.0,fy=0.0;
   force_intr_atmc(0.0, 0,&fx ,&fy);
}
//[M3R2] routine to calculate force in a potnl for all particles @@@@
void  c_list_force_pot()
{
   int i=0,cid_of_pid=0,cc=0,pc=0,pid_in_cell=0;
   double distance=0.0,fx=0.0,fy=0.0,sumx=0.0,sumy=0.0;
   
   for(i=0;i<totalN;i++){
      sumx=0.0,sumy=0.0;
      cid_of_pid = map[i].cid;
      cnbz(cid_of_pid); //calculate cell nbr's index of mvd particle's cell
      for( cc=0;cc<9;cc++){	
         for( pc=0;pc<clen[cnn[cc]];pc++){	
            if(cnn[cc]>=ncell){
               printf("cell list allocation exceded in pot\n\n%d\t%d\n",ncell,cnn[cc]);exit(0);
            }
            pid_in_cell = c_list[cbgn[cnn[cc]]+pc];
            if(i != pid_in_cell){	
               //dis_old = PBCr( (p[pid_in_cell].ri0[0]-p[index].ri0[0]), (p[pid_in_cell].ri0[1]-p[index].ri0[1]) );
               distance = PBCr( (p[pid_in_cell].ri[0]-p[i].ri[0]), (p[pid_in_cell].ri[1]-p[i].ri[1]) );
               //dis_new = PBCr( (p[pid_mvd].ri[0]-p[pid1].ri[0]), 	(p[pid_mv].ri[1]-p[pid1].ri[1]) 	 );
               if(distance < cutoff){
                  force_intr_atmc(distance, i, &fx, &fy);
                  sumx += (pbcx(p[i].ri[0]-p[pid_in_cell].ri[0])/distance)*fx;
                  sumy += (pbcy(p[i].ri[1]-p[pid_in_cell].ri[1])/distance)*fy;
               }
               //if(dis_new < cutoff)
               //{
               //   energyN = energyN + potential( dis_new,index);
               //}
            }
         }
      }
      p[i].forceC[0] = sumx; p[i].forceC[1] = sumy;
	}
}


//#####################################################################################################################################//
// [ 3 ] ############################################### totalE (5sep2017) ###########################################################//
//###################################################################################################################################//
double c_list_totalE(char *string)
{
	double totalE=0,dis=0,totalEharm1=0.0,totalEharm2=0.0;
	int i=0,pid_in_cell=0,pid_mvd=0,bondno=0,Dbondno=0,dirctn=0;
	int dummy_index=0,pc,cc;
	
	if(strcmp(string,"oldCord")==0){	//checek for old config
		if(strcmp(pot_type,"sq_ntwrk")==0){ //check for pot_type square-network, need to put traing-network if want to genralize
			for(i=0;i<totalN;i++){	
				for(dirctn = 1;dirctn<5;dirctn++){	
					bondno =	p[i].nn_static[dirctn];
					totalEharm1 = totalEharm1 + 0.5*k1*pow((PBCr( (p[bondno].ri0[0]-p[i].ri0[0]) , (p[bondno].ri0[1]-p[i].ri0[1]))-a),2) ;        
				}
				for(dirctn = 5;dirctn<9;dirctn++){
					Dbondno =	p[i].nn_static[dirctn];
					totalEharm2 = totalEharm2 + (0.5*k2*pow((PBCr( (p[Dbondno].ri0[0]-p[i].ri0[0]) , (p[Dbondno].ri0[1]-p[i].ri0[1]))-(sqrt(2.0)*a)),2));
				}
			}
			totalE = totalEharm1+totalEharm2;
		}
      else if(strcmp(pot_type,"no_pot")==0){
         totalE = 0.0;
      }else if(strcmp(pot_type,"triang_ntwrk")==0){
         for(i=0;i<totalN;i++){	
            for(dirctn = 1;dirctn<7;dirctn++){	
               bondno =	p[i].nn_static[dirctn];
               totalEharm1 = totalEharm1 + 0.5*k1*pow((PBCr( (p[bondno].ri0[0]-p[i].ri0[0]) , (p[bondno].ri0[1]-p[i].ri0[1]))-a),2) ;        
            }
            totalE = totalEharm1;
         }
      }else	if(strcmp(pot_type,"honeycomb_ntwrk")==0){ //check for pot_type square-network, need to put traing-network if want to genralize
         for(i=0;i<totalN;i++){	
            for(dirctn = 1;dirctn<4;dirctn++){	
               //bondno =	p[i].nbr_lst[dirctn];
               bondno =	p[i].nn_static[dirctn];
               totalEharm1 = totalEharm1 + 0.5*k1*pow((PBCr( (p[bondno].ri0[0]-p[i].ri0[0]) , (p[bondno].ri0[1]-p[i].ri0[1]))-a),2) ;        
            }
            for(dirctn = 4;dirctn<10;dirctn++){
               //Dbondno =	p[i].nbr_lst[dirctn];
               Dbondno =	p[i].nn_static[dirctn];
               totalEharm2 = totalEharm2 + (0.5*k2*pow((PBCr( (p[Dbondno].ri0[0]-p[i].ri0[0]) , (p[Dbondno].ri0[1]-p[i].ri0[1]))-(sqrt(3.0)*a)),2));
            }
         }
         totalE = totalEharm1+totalEharm2;
      }else	if(strcmp(pot_type,"kagome_ntwrk")==0){ //check for pot_type square-network, need to put traing-network if want to genralize
         for(i=0;i<totalN;i++){	
            for(dirctn = 1;dirctn<5;dirctn++){	
               //bondno =	p[i].nbr_lst[dirctn];
               bondno =	p[i].nn_static[dirctn];
               totalEharm1 = totalEharm1 + 0.5*k1*pow((PBCr( (p[bondno].ri0[0]-p[i].ri0[0]) , (p[bondno].ri0[1]-p[i].ri0[1]))-(a*0.5)),2) ;        
            }
            for(dirctn = 5;dirctn<9;dirctn++){
               //Dbondno =	p[i].nbr_lst[dirctn];
               Dbondno =	p[i].nn_static[dirctn];
               totalEharm2 = totalEharm2 + (0.5*k2*pow((PBCr( (p[Dbondno].ri0[0]-p[i].ri0[0]) , (p[Dbondno].ri0[1]-p[i].ri0[1]))-(sqrt(3.0)*a*0.5)),2));
            }
         }
         totalE = totalEharm1+totalEharm2;
      }
      else{ // do it for other defined pot_type
			for(i=0;i<totalN;i++){
				cnbz(map[i].cid);
				for(cc=0;cc<9;cc++){	
					for(pc=0;pc<clen[cnn[cc]];pc++){	
						pid_in_cell = c_list[cbgn[cnn[cc]]+pc];	pid_mvd = i;
						if(pid_mvd != pid_in_cell){ //to avoid V(r_ii)
							//dis = PBCr( (p[pid1].ri0[0]-p[pid_mv].ri0[0]), (p[pid1].ri0[1]-p[pid_mv].ri0[1]) );
							dis = PBCr( (p[pid_mvd].ri0[0]-p[pid_in_cell].ri0[0]), 	(p[pid_mvd].ri0[1]-p[pid_in_cell].ri0[1]) 	 );
							if(dis < cutoff){
								totalE = totalE + potential( dis,dummy_index);
							}
						}
					}
				}
			}
		}
      printf(WHT"\n Initial configurational energy of the system is % .8lf\n", totalE*0.5);
	}
	else if(strcmp(string,"newCord")==0){	//check for new config
		if(strcmp(pot_type,"sq_ntwrk")==0){
			for(i=0;i<totalN;i++){	
				for(dirctn = 1;dirctn<5;dirctn++){
					bondno =	p[i].nn_static[dirctn];
					totalEharm1 = totalEharm1 + 0.5*k1*pow((PBCr( (p[bondno].ri[0]-p[i].ri[0]) , (p[bondno].ri[1]-p[i].ri[1]))-a),2) ;        
					//printf("%d\t%d\t%lf\n",i,bondno,pow((PBCr( (p[bondno].ri[0]-p[i].ri[0]) , (p[bondno].ri[1]-p[i].ri[1]))-a),2)  );
				}
				for(dirctn = 5;dirctn<9;dirctn++){
					Dbondno =	p[i].nn_static[dirctn];
					totalEharm2 = totalEharm2 + (0.5*k2*pow((PBCr( (p[Dbondno].ri[0]-p[i].ri[0]) , (p[Dbondno].ri[1]-p[i].ri[1]))-(sqrt(2.0)*a)),2));
				}
			}
         totalE = totalEharm1+totalEharm2;
		}else if(strcmp(pot_type,"no_pot")==0){
         totalE = 0.0;
      }else if(strcmp(pot_type,"triang_ntwrk")==0){
         double rij_x=0.0,rij_y=0.0;
         for(i=0;i<totalN;i++){	
            for(dirctn = 1;dirctn<7;dirctn++){	
               bondno =	p[i].nn_static[dirctn];
               rij_x = (p[i].ri[0] - p[bondno].ri[0]);
               rij_y = (p[i].ri[1] - p[bondno].ri[1]);
               totalEharm1 += 0.5*k1*pow((PBCr(rij_x,rij_y)-a),2); // apprx working
            }
            totalE = totalEharm1;
         }
      }else if(strcmp(pot_type,"honeycomb_ntwrk")==0){
         for(i=0;i<totalN;i++){
            for(dirctn = 1;dirctn<4;dirctn++){
               //bondno =	p[i].nbr_lst[dirctn];
               bondno =	p[i].nn_static[dirctn];
               totalEharm1 = totalEharm1 + 0.5*k1*pow((PBCr( (p[bondno].ri[0]-p[i].ri[0]) , (p[bondno].ri[1]-p[i].ri[1]))-a),2) ;        
               //printf("%d\t%d\t%lf\n",i,bondno,pow((PBCr( (p[bondno].ri[0]-p[i].ri[0]) , (p[bondno].ri[1]-p[i].ri[1]))-a),2)  );
            }
            for(dirctn = 4;dirctn<10;dirctn++){
               //Dbondno =	p[i].nbr_lst[dirctn];
               Dbondno =	p[i].nn_static[dirctn];
               totalEharm2 = totalEharm2 + (0.5*k2*pow((PBCr( (p[Dbondno].ri[0]-p[i].ri[0]) , (p[Dbondno].ri[1]-p[i].ri[1]))-(sqrt(3.0)*a)),2));
            }
         }
         totalE = totalEharm1+totalEharm2;
      }else if(strcmp(pot_type,"kagome_ntwrk")==0){
         for(i=0;i<totalN;i++){	
            for(dirctn = 1;dirctn<5;dirctn++){	
               //bondno =	p[i].nbr_lst[dirctn];
               bondno =	p[i].nn_static[dirctn];
               totalEharm1 = totalEharm1 + 0.5*k1*pow((PBCr( (p[bondno].ri[0]-p[i].ri[0]) , (p[bondno].ri[1]-p[i].ri[1]))-(a*0.5)),2) ;        
               //printf("%d\t%d\t%lf\n",i,bondno,pow((PBCr( (p[bondno].ri[0]-p[i].ri[0]) , (p[bondno].ri[1]-p[i].ri[1]))-a),2)  );
            }
            for(dirctn = 5;dirctn<9;dirctn++){
               //Dbondno =	p[i].nbr_lst[dirctn];
               Dbondno =	p[i].nn_static[dirctn];
               totalEharm2 = totalEharm2 + (0.5*k2*pow((PBCr( (p[Dbondno].ri[0]-p[i].ri[0]) , (p[Dbondno].ri[1]-p[i].ri[1]))-(sqrt(3.0)*a*0.5)),2));
            }
         }
         totalE = totalEharm1+totalEharm2;
      }
      else{	// do it for other defined pot_type
			for(i=0;i<totalN;i++){	
				cnbz(map[i].cid);
				for(cc=0;cc<9;cc++){	
					for(pc=0;pc<clen[cnn[cc]];pc++){
						pid_in_cell = c_list[cbgn[cnn[cc]]+pc];	pid_mvd = i;
						if(pid_mvd != pid_in_cell){	//to avoid V(r_ii)
							//dis = PBCr( (p[pid1].ri0[0]-p[pid_mv].ri0[0]), (p[pid1].ri0[1]-p[pid_mv].ri0[1]) );
							dis = PBCr( (p[pid_mvd].ri[0]-p[pid_in_cell].ri[0]), 	(p[pid_mvd].ri[1]-p[pid_in_cell].ri[1]) 	 );
							if(dis < cutoff){
								totalE = totalE + potential( dis,dummy_index);
							}
						}
					}
				}
			}
		}
   }
//printf(WHT"\n Initial configurational energy of the system is % .8lf\n", totalE*0.5);
//printf(MAG" Warning ! energy calculated is using cell_list method, only sq_ntwrk, triang_ntwrk and polar potential with cutoff is supported, rest is garbag.\n"RESET);
return (totalE*0.5);
}
//#####################################################################################################################################//








//#####################################################################################################################################//
// [ 4 ] ########################################## cell_list update partial #########################################################//
//###################################################################################################################################//
void c_list_update_partial(int idx)
{
	double new_cord_x=0.0, new_cord_y=0.0;
	int new_cx=0, new_cy=0, new_cid=0, old_cid=0, pos_in_c_list=0, temp=0,i=0,cc=0,dummy;
	
   new_cord_x = p[idx].ri[0];
   new_cord_y = p[idx].ri[1];
   new_cx = floor((( (p[idx].ri[0]+(lbox_x*0.5))/(lcell_x) )));
   new_cy = floor((( (p[idx].ri[1]+(lbox_y*0.5))/(lcell_y) )));
   //new_cid = (new_cx+ncell_x)%ncell_x + (((new_cy+ncell_x)%ncell_x)*(ncell_x));
   new_cid = (new_cx+ncell_x)%ncell_x + (((new_cy+ncell_y)%ncell_y)*(ncell_x));	//works for rectnglr system also
   old_cid =	map[idx].cid;
   //printf("new cid %d\t%d\n",new_cid,old_cid);
   //get the postn of particle in the cell.
   for(i=cbgn[old_cid];i<cbgn[old_cid]+clen[old_cid];i++){
      if(map[idx].pid == c_list[i]){
         pos_in_c_list = i;
         break;
      }
   }
   if(new_cid < old_cid){
      temp = map[idx].pid;
      for(i=pos_in_c_list;i> cbgn[new_cid]+clen[new_cid]; i--){	
         c_list[i] = c_list[i-1];
      }
      c_list[(cbgn[new_cid]+clen[new_cid])] = temp;
      clen[new_cid]++;
      clen[old_cid]--;
      for(cc=old_cid;cc>new_cid;cc--){
         cbgn[cc]++;
      }
      map[idx].cid = new_cid;
   }else if(new_cid > old_cid){
      temp = map[idx].pid;
      for( i=pos_in_c_list;i< cbgn[new_cid]-1; i++){	
         c_list[i] = c_list[i+1];
      }
      c_list[cbgn[new_cid]-1] = temp;
      clen[new_cid]++;
      clen[old_cid]--;
      for(cc=new_cid;cc>old_cid;cc--){
         cbgn[cc]--;
      }
   
      map[idx].cid = new_cid;
   }else if(new_cid == old_cid){
      dummy = 3;
      //printf("Error in cid");
   }		
}

//#####################################################################################################################################//




//*****************************************************************************************************************//
//********************************************* [ Aux funcs in general ]*******************************************//
//***************************************************************************************************************//
// [1] 
void sort_acc_cid(int pid2cid[][2])
{
	int g=0,h=0;
   
	for( g=1;g<totalN;g++){
      for( h=0; h<(totalN-1); h++){
         if(pid2cid[h][1]>pid2cid[h+1][1]){
            int temp=pid2cid[h][1];
            pid2cid[h][1]=pid2cid[h+1][1];
            pid2cid[h+1][1]=temp;
            int temp1 = pid2cid[h][0];
            pid2cid[h][0] = pid2cid[h+1][0];
            pid2cid[h+1][0]=temp1;
         }
      }
   }
}

// [2] cell list nbrs
void cnbz(int id_cell)
{ 
	//int N = sqrt(ncell);
	//int x = id_cell%N;
	//int y = floor(id_cell/N);
	
	int x = (id_cell%ncell_x); int y = floor(id_cell/ncell_x);
	cnn[0]	=	x + ncell_x*y;
	cnn[1]	=	(x+ncell_x)%ncell_x 		+ ((y+ncell_y+1)%ncell_y)*ncell_x;	//| E
	cnn[2]	=	(x+ncell_x)%ncell_x 		+ ((y+ncell_y-1)%ncell_y)*ncell_x;	//| W
	cnn[3]	=	(x+ncell_x+1)%ncell_x 	+ ((y+ncell_y)%ncell_y)*ncell_x;	   //| N
	cnn[4]	=	(x+ncell_x-1)%ncell_x 	+ ((y+ncell_y)%ncell_y)*ncell_x;	   //| S
	cnn[5]	=	(x+ncell_x+1)%ncell_x 	+ ((y+ncell_y+1)%ncell_y)*ncell_x;	//| NE
	cnn[6]	=	(x+ncell_x-1)%ncell_x 	+ ((y+ncell_y+1)%ncell_y)*ncell_x;	//| SE
	cnn[7]	=	(x+ncell_x+1)%ncell_x 	+ ((y+ncell_y-1)%ncell_y)*ncell_x;	//| NW
	cnn[8]	=	(x+ncell_x-1)%ncell_x 	+ ((y+ncell_y-1)%ncell_y)*ncell_x;	//| SW
	
}
/*****************************************************************************************************************/
//***************************************************************************************************************//
//***************************************************************************************************************//


//////////////////////////////////////////////// End of Cell list ////////////////////////////////////////////////////////
