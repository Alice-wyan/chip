#pragma once

#include "lib.h"


class CLine
{
public:
	CLine(){

	}
void clear(){
	//cout<<"clear line"<<endl;
       s11Xy  =0;
       s12Xy  =0;
       s22Xy  =0;
       g1Xy   =0;
       g2Xy   =0;     	
       ddXy=0;
	   a1Xy=0;
	   a2Xy =0;	
	   nHits=0;
}	
double k(){
return a2Xy;
}	
double b(){
return a1Xy;
}	

int vertical(){
if(ddXy==0) return 1; 
return 0;
}	

void add(CEdgeHit *a){
	add(a->x,a->y);
}
void add( double x,double y){

  s11Xy = s11Xy +  1 ;
  s12Xy = s12Xy +  x ;
  s22Xy = s22Xy +  x * x ;
  g1Xy  = g1Xy  +  y ;
  g2Xy  = g2Xy  +  x * y ;
  nHits++;

  if ( nHits >= 2  )
  {
     ddXy  = s11Xy * s22Xy -  s12Xy*s12Xy  ;
     if ( ddXy != 0 ) {
        a1Xy  = ( g1Xy * s22Xy - g2Xy * s12Xy ) / ddXy ;//b
        a2Xy  = ( g2Xy * s11Xy - g1Xy * s12Xy ) / ddXy ;//k

		//cout<<"k="<<a2Xy<<endl;
		//cout<<"b="<<a1Xy<<endl;
     }
  }	
}	



double check( CEdgeHit *a){
     return check(a->x,a->y);
}

double check( double x,double y){

     if ( ddXy != 0 ) {
   temp = (a2Xy * x - y + a1Xy) ;
   dxy  = temp * temp / ( a2Xy * a2Xy + 1.F ) ;
     }
     else {
   dxy=fabs(s12Xy/s11Xy-x);
   dxy*=dxy;
     }
	return dxy;
}

	
   double    s11Xy  ;
   double    s12Xy  ;
   double    s22Xy  ;
   double    g1Xy   ;
   double    g2Xy   ;     	
   double    ddXy, a1Xy, a2Xy ;
   double    ddXx;
   int nHits;
   double   temp;
   double   dxy;

   CEdgeHit *first;
   CEdgeHit *last;

};

