#include "global.h"

double pbcx(double dx)
{
	
   if(fabsf(dx)>(lbox_x*0.5)) {
		if(dx > lbox_x*0.5)
         { dx = dx - lbox_x; }
		if(dx <= -lbox_x*0.5)
			{ dx = dx + lbox_x; }
	}
  
return dx;
}
double pbcy(double dy)
{ 
	
   if(fabsf(dy)>(lbox_y*0.5)) {
		if(dy > lbox_y*0.5)
			{ dy = dy - lbox_y; }
		if(dy <= -lbox_y*0.5)
			{ dy = dy + lbox_y; }
	}

return dy;
}

double PBCr(double rx, double ry)
{
	
	//periodicity in x
	if(fabsf(rx)>(lbox_x*0.5)){
		if(rx > lbox_x*0.5)
			{ rx = rx - lbox_x; }
		if(rx <= -lbox_x*0.5)
			{ rx = rx + lbox_x; }
	}
	//periodicity in y 
	if(fabsf(ry)>(lbox_y*0.5)){
		if(ry > lbox_y*0.5)
			{ ry = ry - lbox_y; }
		if(ry <= -lbox_y*0.5)
			{ ry = ry + lbox_y; }
	}
	
return sqrt( pow(rx,2) + pow(ry,2) );
}

