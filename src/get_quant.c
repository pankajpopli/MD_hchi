#include"global.h"

double calc_Temp()
{
   int i=0;
   double sum=0.0;
      
   for(i=0;i<totalN;i++){
      sum += mass*(((p[i].vel[0]*p[i].vel[0]) + (p[i].vel[1]*p[i].vel[1]))*mass)/(2.0*K_B);
   }

return (sum/totalN);
}
