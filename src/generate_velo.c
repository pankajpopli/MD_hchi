#include"global.h"

void generate_velo(const gsl_rng *ran)
{
   int i=0;
   double rnd_angle=0.0;
   FILE *f;
   
   for(i=0;i<totalN;i++){
      //rnd_angle = 2.0*PI*gsl_rng_uniform(ran);
      rnd_angle = 2.0*PI*((gsl_rng_uniform_int(ran,(totalN*10+1))/(totalN*10.0)));
      p[i].vel[0] = Vck0*cos(rnd_angle);
      p[i].vel[1] = Vck0*sin(rnd_angle);
   }
}
