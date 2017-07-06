#pragma once

#include "lib.h"

class CRegionMatch
{
public:
	CRegionMatch(){
    regionLength=3;
	ndx=10;
	ndy=10;
	nds=10;
	dx=1;
	dy=1;
	ds=PI/20;
	}

void setup(CFrame *_psfr,CFrame *_pdfr,CEdgeHitChain *_pv,CPointProjection * _pcpp){
	psfr=_psfr;
	pdfr=_pdfr;
	pv=_pv;
	pcpp=_pcpp;
}

double minDiff(){

	double min=100000;

	cpp.setup(pcpp);
	double s=atan2(cpp.sin,cpp.cos);
	double cpps;
	for(int k=-nds/2;k<nds/2;k++){
		cpps=s+k*ds;
		pcpp->sin=sin(cpps);
		pcpp->cos=cos(cpps);    
	for(int i=-ndx/2;i<ndx/2;i++){
		pcpp->dx=cpp.dx+i*dx;
	for(int j=-ndy/2;j<ndy/2;j++){
		pcpp->dy=cpp.dy+j*dy;
		pcpp->setupFirst();
		diff();
		//if(i==0&&j==0&&k==0){
		//cout<<"sum="<<sum<<endl;
		//}
		if(sum<min){
		min=sum;
		minx=pcpp->dx;
		miny=pcpp->dy;
		mins=pcpp->sin;
		minc=pcpp->cos;

		}
	}}}
	pcpp->setup(&cpp);
	return min;
}

double diff(){
 sumb=0;
 sumb2=0;
 sumg=0;
 sumg2=0;
 sumr=0;
 sumr2=0;
 sum=0;
 n=0;
 double b;
 double g;
 double r;
 
	for(unsigned int i=0;i<x.size();i++){
		pcpp->expect(x[i],y[i]);
	if(pcpp->x<0) continue;
	if(pcpp->x>=psfr->width) continue;
	if(pcpp->y<0) continue;
	if(pcpp->y>=psfr->height) continue;	
    b=psfr->b(x[i],y[i])-pdfr->bf(pcpp->x,pcpp->y);
    g=psfr->g(x[i],y[i])-pdfr->gf(pcpp->x,pcpp->y);
    r=psfr->r(x[i],y[i])-pdfr->rf(pcpp->x,pcpp->y);
	sumb+=b;
	sumb2+=b*b;
	sumg+=g;
	sumg2+=g*g;
	sumr+=r;
	sumr2+=r*r;
	n++;
	}
	sumb/=n;
	sumb2/=n;
	sumg/=n;
	sumg2/=n;
	sumr/=n;
	sumr2/=n;
	
	sum=sumb2-sumb*sumb+sumg2-sumg*sumg+sumr2-sumr*sumr;
	sum=sqrt(sum);
	return sum;
    
}


void generate(){
	x.clear();
	y.clear();
	
for(CEdgeHit*a=pv->first;a!=NULL;a=a->next){
for(int i=1;i<regionLength;i++){
	if(!find(a->x,a->y)) fill(a->x,a->y);
	if(!find(a->x+i,a->y)) fill(a->x+i,a->y);
	if(!find(a->x-i,a->y)) fill(a->x-i,a->y);
	if(!find(a->x,a->y+i)) fill(a->x,a->y+i);
	if(!find(a->x,a->y-i)) fill(a->x,a->y-i);
}	
}

}
bool find(int _x,int _y){
	for(unsigned int i=0;i<x.size();i++){
		if(x[i]==_x&&y[i]==_y) return true;
	}
	return false;
}	

bool fill(int _x,int _y){
	if(_x<0) return false;
	if(_x>=psfr->width) return false;
	if(_y<0) return false;
	if(_y>=psfr->height) return false;
	x.push_back(_x);
	y.push_back(_y);
	return true;
}
	
public:
	
CFrame *psfr;
CFrame *pdfr;
CEdgeHitChain *pv;

vector<int> x;
vector<int> y;
int regionLength;
CPointProjection * pcpp;
CPointProjection cpp;

double sumb;
double sumb2;
double sumg;
double sumg2;
double sumr;
double sumr2;
double sum;
double n;

int ndx;
int ndy;
int nds;
double dx;
double dy;
double ds;
double minx;
double miny;
double mins;
double minc;

};

