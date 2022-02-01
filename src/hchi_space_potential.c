#include"global.h"

///[1] Non-affine field h{Ri}
double h_field(double distance)
{
	double patch_size=1.5*a, Vchi=0, lambda=1;
	//return (hchi + Vchi*exp( -((fabsf(distance))/(lambda))    ) );
	double ref_dist = distance;
	if(distance <= patch_size*sqrt(2)*a){
		//return Vchi;
		return (hchi_0 + Vchi*exp( -((fabsf(ref_dist))/(lambda))    ) );
	}else{
		//return (hchi + Vchi*exp( -((fabsf(ref_dist))/(lambda))    ) );
		//return Vchi;
		return 0.0;
	}
	printf(RED"we are not using space varying h_field | potential_polar_form.c,"RESET);exit(0);
}

void spatial_random_h(const gsl_rng *ran3){
	printf(GRN"\nInitialising random hchi values\n"RESET);
	int i=0;
	double std_a=1;
		for(i=0;i<totalN;i++){
			std_a = (gsl_ran_gaussian_ziggurat(ran3, hchi_width) + hchi_mean);
			p[i].h_xy	= std_a;
		}
}