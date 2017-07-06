#pragma once

#include "lib.h"

class CPointProjection
{
public:
	CPointProjection(){

	}

void setup(CPointProjection * pcpp){
	setup(pcpp->dx,pcpp->dy,pcpp->sin,pcpp->cos,pcpp->order,pcpp->pv);
}	

void setup(double _dx,double _dy,double _sin,double _cos,int _order,CEdgeHitChain *_pv){
	dx=_dx;
	dy=_dy;
	sin=_sin;
	cos=_cos;
	order=_order;
	pv=_pv;
	
setupFirst();	
	
}	

void setupFirst(){

if(order==0){
    sf.x=dx+pv->first->x;
    sf.y=dy+pv->first->y;
}else{
    sf.x=dx+pv->last->x;
    sf.y=dy+pv->last->y;	
}	
	
}	


void expect(double xp,double yp){
if(order==0){
	x=(xp-pv->first->x)*cos-(yp-pv->first->y)*sin+sf.x;
	y=(xp-pv->first->x)*sin+(yp-pv->first->y)*cos+sf.y;
}else{
	x=(xp-pv->last->x)*cos-(yp-pv->last->y)*sin+sf.x;
	y=(xp-pv->last->x)*sin+(yp-pv->last->y)*cos+sf.y;	
}
}
	
public:
	
double dx;
double dy;
double sin;
double cos;
int order;	
CEdgeHit sf;
CEdgeHitChain *pv;

double x;
double y;

};

