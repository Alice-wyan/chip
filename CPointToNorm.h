#pragma once

#include "lib.h"

class CPointToNorm
{
public:
CPointToNorm(){
	ini();

};


~CPointToNorm(){
	clear(); 

};

void ini(){

	outFlag=0;
	frIn.f=NULL;
}

void setup(IplImage *pFrame,int max=20){
	setupFrame(pFrame);
	setupOut(max);
}



void setupFrame(IplImage *pFrame){
	frIn.set(pFrame);
}


void setupOut(int max){
 normLength=max;

 if(outFlag==1) delete []out.Point;
  out.Point=new CPoint[normLength];
  outFlag=1;
 
}

void clear(){

	 if(outFlag==1) delete []out.Point;  
	outFlag=0;

}





void transform(CPoint &p,int chl=0){

	if(frIn.f==NULL){
	cout<<"CPointToNorm::frIn.f==NULL"<<endl;
	exit(0);
	}

	frIn.grad(p.x,p.y,chl,vx,vy);
	double v=sqrt(vx*vx+vy*vy);
	if(v==0){
	vx=sqrt(0.5);
	vy=sqrt(0.5);
	}else{
	vx=vx/v;
	vy=vy/v;
	}
	for(int i=0;i<normLength;i++){
		out.Point[i].x=p.x+vx*(i-normLength/2);  
		out.Point[i].y=p.y+vy*(i-normLength/2);
	}
	out.nPoints=normLength;
}



void transform(CPoint &p,double n,double dx,double dy){

	if(frIn.f==NULL){
	cout<<"CPointToNorm::frIn.f==NULL"<<endl;
	exit(0);
	}


	double v=sqrt(dx*dx+dy*dy);
	if(v==0){
	vx=sqrt(0.5);
	vy=sqrt(0.5);
	}else{
	vx=dx/v;
	vy=dy/v;
	}
	for(int i=0;i<normLength;i++){
		out.Point[i].x=p.x-vy*(i-normLength/2)+n*vx;   //clockwise 90 degree
		out.Point[i].y=p.y+vx*(i-normLength/2)+n*vy;
	}
	out.nPoints=normLength;
}


void print(CFrame &fr){

	for(int i=0;i<out.nPoints;i++){
		if(out.Point[i].x>=fr.width) continue;
		if(out.Point[i].x<0) continue;
		if(out.Point[i].y>=fr.height) continue;
		if(out.Point[i].y<0) continue;

		*(fr.pb(int(out.Point[i].x),int(out.Point[i].y)))=255;
		*(fr.pg(int(out.Point[i].x),int(out.Point[i].y)))=0;
		*(fr.pr(int(out.Point[i].x),int(out.Point[i].y)))=0;

	
	}

}

public:
CProfile out;
int outFlag;
int normLength;
CFrame frIn;
double vx;
double vy;

};