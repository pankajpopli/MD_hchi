///////////// calculates p(chi) for square lattice /////////////////
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<stddef.h>
#define sigma 1.0
#define N 32
#define red_rho 1.0
#define a (sigma/(sqrt(red_rho)))
#define skipconfig 300
#define timesteps 1999
#define dimen 2
#define totalN	(N*N)
#define x_size (N*a)
#define y_size (N*a)
struct Atom{
				double R[dimen];
				double r[dimen];
};
int nn[8];
void nearN(int x, int y);
double pbcx(double dx);
double pbcy(double dy);


int main()
{
	struct Atom atom[totalN]={0.0};
	FILE *f1;
	char file[128];
	int NN[totalN*6][2];
	int i,j,k,xi,yi,neighbours=8;
	double detY,detY_inv;
	double chi,X[2][2]={0.0},Y[2][2]={0.0},eps[2][2]={0.0};
	
	f1 = fopen("cordinates.dat","r");
	for(i=0;i<totalN;i++)
		{
			fscanf(f1,"%*d%lf\t%lf\n",&atom[i].R[0],&atom[i].R[1]);
		}
	fclose(f1);
	


//-------------------------------
	for(i=skipconfig;i<timesteps;i++)
	{
		sprintf (file,"time/new_cordinates_%d.dat", i);
		f1 = fopen (file,"r");
		for(j=0;j<totalN;j++)
		{
			fscanf(f1,"%*d\t%lf\t%lf\t%*f\t%*f\t%*f\t%*f\t%*f\n",&atom[j].r[0],&atom[j].r[1]);
			//printf("%lf\t%lf\n",atom[j].r[0],atom[j].r[1]);
		}fclose(f1);
		
		for(j=0;j<totalN;j++)
		{
			X[0][0] = X[1][0] = X[0][1] = X[1][1] = 0;
         Y[0][0] = Y[1][0] = Y[0][1] = Y[1][1] = 0;
         chi = 0;
         
			xi = j%N;
			yi = floor(j/N);
			nearN(xi,yi); //find nearN and storres in nn[]
			for(k=0;k<neighbours;k++)
			{
				X[0][0] += (pbcx(atom[nn[k]].r[0] - atom[j].r[0]))*(pbcx(atom[nn[k]].R[0] - atom[j].R[0]));
				X[0][1] += (pbcx(atom[nn[k]].r[0] - atom[j].r[0]))*(pbcy(atom[nn[k]].R[1] - atom[j].R[1]));
				X[1][0] += (pbcy(atom[nn[k]].r[1] - atom[j].r[1]))*(pbcx(atom[nn[k]].R[0] - atom[j].R[0]));
				X[1][1] += (pbcy(atom[nn[k]].r[1] - atom[j].r[1]))*(pbcy(atom[nn[k]].R[1] - atom[j].R[1]));
				
				Y[0][0] += (pbcx(atom[nn[k]].R[0] - atom[j].R[0]))*(pbcx(atom[nn[k]].R[0] - atom[j].R[0]));
				Y[0][1] += (pbcx(atom[nn[k]].R[0] - atom[j].R[0]))*(pbcy(atom[nn[k]].R[1] - atom[j].R[1]));
				Y[1][0] += (pbcy(atom[nn[k]].R[1] - atom[j].R[1]))*(pbcx(atom[nn[k]].R[0] - atom[j].R[0]));
				Y[1][1] += (pbcy(atom[nn[k]].R[1] - atom[j].R[1]))*(pbcy(atom[nn[k]].R[1] - atom[j].R[1]));
			}
			
			detY = (Y[0][0]*Y[1][1]) - (Y[0][1]*Y[1][0]) ; detY_inv = 1.0/detY;
			
			eps[0][0]	=	-X[0][1]*Y[1][0]*detY_inv + X[0][0]*Y[1][1]*detY_inv - 1;
			eps[0][1]	=	 X[0][1]*Y[0][0]*detY_inv - X[0][0]*Y[0][1]*detY_inv;
			eps[1][0]	=	-X[1][1]*Y[1][0]*detY_inv + X[1][0]*Y[1][1]*detY_inv;
			eps[1][1]	=	 X[1][1]*Y[0][0]*detY_inv - X[1][0]*Y[0][1]*detY_inv - 1;
			
			for(k=0;k<neighbours;k++)
			{
				chi += pow((pbcx(atom[nn[k]].r[0]-atom[j].r[0]) - (1.0 +eps[0][0])*(pbcx(atom[nn[k]].R[0] - atom[j].R[0])) -(eps[0][1])*(pbcy(atom[nn[k]].R[1] - atom[j].R[1]))),2) +  pow((pbcy(atom[nn[k]].r[1]-atom[j].r[1]) - (1 +eps[1][1])*(pbcy(atom[nn[k]].R[1] - atom[j].R[1])) -(eps[1][0])*(pbcx(atom[nn[k]].R[0] - atom[j].R[0]))),2);
			}
			printf("%.10lf\n",chi);
		}
	 }
}

//########################################################### AUXILAARY FUNCTIONS ################################################################
// 1) find nearest neighbour traingular lattice
void nearN(int x, int y)
{ 
	//nn[0]	=	i + N*j;
	nn[0]	=	(x+N)%N + ((y+N+1)%N)*N;	//E
	nn[1]	=	(x+N)%N + ((y+N-1)%N)*N;	//W
	nn[2]	=	(x+N+1)%N + ((y+N)%N)*N;	//N
	nn[3]	=	(x+N-1)%N + ((y+N)%N)*N;	//S
	nn[4]	=	(x+N+1)%N + ((y+N+1)%N)*N;	//NE
	nn[5]	=	(x+N-1)%N + ((y+N+1)%N)*N;	//SE
	nn[6]	=	(x+N+1)%N + ((y+N-1)%N)*N;	//NW
	nn[7]	=	(x+N-1)%N + ((y+N-1)%N)*N;	//SW
	
}

// 2) periodic pundary condition for coordinates 
double pbcx(double dx)
{
	if(fabsf(dx)>(x_size*0.5)) 
	{
		if(dx > x_size*0.5)
			{ dx = dx - x_size; }
		if(dx <= -x_size*0.5)
			{ dx = dx + x_size; }
	}
return dx;
}
double pbcy(double dy)
{ 
	if(fabsf(dy)>(y_size*0.5)) 
	{
		if(dy > y_size*0.5)
			{ dy = dy - y_size; }
		if(dy <= -y_size*0.5)
			{ dy = dy + y_size; }
	}
return dy;
}
//--------------------------------------------- END OF PROGRAM ---------------------------------
