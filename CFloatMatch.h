#pragma once

#include "lib.h"


class CFloatMatch
{
public:
CFloatMatch(){

};
~CFloatMatch(){

};

double cross(CPoint &a,CPoint &b){
	return a.x*b.y-b.x*a.y;
}

CPoint center(CProfile *cp){
	double sumx=0;
	double sumw=0;
	double sumy=0;

	
	for(int i=0;i<cp->nPoints;i++){
		sumx+= cp->PointDistance(i,(i+1)%cp->nPoints)*(cp->Point[i].x+cp->Point[(i+1)%cp->nPoints].x)/2;
		sumy+= cp->PointDistance(i,(i+1)%cp->nPoints)*(cp->Point[i].y+cp->Point[(i+1)%cp->nPoints].y)/2;
		sumw+=cp->PointDistance(i,(i+1)%cp->nPoints);
	}
	CPoint c;
	c.x=sumx/sumw;
	c.y=sumy/sumw;
	return c;
}

void setup(CProfile *in1O,CProfile *in2O){
    in0=in1O;
    in1=in2O;
}


double transform(int nsita){
	min=nsita*10;
	return match_contour(nsita,min);
}

double match_contour(int nsita,double &min){

	CPoint c0=center(in0);
	CPoint c1=center(in1);

	xc0=c0.x;
	yc0=c0.y;
	xc1=c1.x;
	yc1=c1.y;


	TProfile *tp_0;
	TProfile *tp_1;


	tp_0=new TProfile("tpname0","tptitle0",nsita,-PI,PI);
	for(int i=0;i<in0->nPoints;i++){
		a.x=in0->Point[i].x-xc0;
		a.y=in0->Point[i].y-yc0;
		b.x=in0->Point[(i+1)%in0->nPoints].x-xc0;
		b.y=in0->Point[(i+1)%in0->nPoints].y-yc0;

		tp_0->Fill(atan2(a.y+b.y,a.x+b.x),cross(a,b));
	}
	tp_1=new TProfile("tpname1","tptitle1",nsita,-PI,PI);
	for(int i=0;i<in1->nPoints;i++){
		a.x=in1->Point[i].x-xc0;
		a.y=in1->Point[i].y-yc0;
		b.x=in1->Point[(i+1)%in1->nPoints].x-xc0;
		b.y=in1->Point[(i+1)%in1->nPoints].y-yc0;

		tp_1->Fill(atan2(a.y+b.y,a.x+b.x),cross(a,b));
	}
	min_pos=match_profile(tp_0,tp_1,min);
	tp_0->Delete();
	tp_1->Delete();

	return double(min_pos*2*PI)/double(nsita);

}


int match_profile(TProfile *t0,TProfile *t1,double &min){
	;
	int n0=t0->GetNbinsX();
	int n1=t1->GetNbinsX();

	if(n0!=n1) return -1;
	int min_pos;
	double f0;
	double f1;
	double sum;


	min_pos=-1;
	for(int i=0;i<n0;i++){
    sum=0;
	for(int j=0;j<n0;j++){
	f0=t0->GetBinContent((i+j+n0)%n0+1);
	f1=t1->GetBinContent(j+1);
	sum+=(f0-f1)*(f0-f1);
	}	
	sum/=n0;
	sum=sqrt(sum);
	if(sum<min){
	min=sum;
	min_pos=i;
	}
	}
	return min_pos;
}


void print(){

	cout<<"CFloatMatch::min_pos="<<min_pos<<endl;
	cout<<"CFloatMatch::min="<<min<<endl;

}

public:

CProfile *in0;
CProfile *in1;

double	min;
int min_pos;
double	xc0;
double	yc0;
double	xc1;
double	yc1;
CPoint a;
CPoint b;

};