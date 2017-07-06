#pragma once

#include "lib.h"

class CPointToCircle
{
public:
CPointToCircle(){
	ini();

};


~CPointToCircle(){
	clear(); 

};

void ini(){

	outFlag=0;
	frIn.f=NULL;
}

void setup(IplImage *pFrame,int max=20,double r=3){
	setupFrame(pFrame);
	setupOut(max);
	radius=r;
}



void setupFrame(IplImage *pFrame){
	frIn.set(pFrame);
}


void setupOut(int max){
 circleLength=max;

 if(outFlag==1) delete []out.Point;
  out.Point=new CPoint[circleLength];
  outFlag=1;
 
}

void clear(){

	 if(outFlag==1) delete []out.Point;  
	outFlag=0;

}





void transform(CPoint &p){

	if(frIn.f==NULL){
	cout<<"CPointToCircle::frIn.f==NULL"<<endl;
	exit(0);
	}
	double ds=2*PI/circleLength;

	for(int i=0;i<circleLength;i++){
		out.Point[i].x=p.x+radius*cos(i*ds);  
		out.Point[i].y=p.y+radius*sin(i*ds); 
		if(out.Point[i].x<0) out.Point[i].x=0;
		if(out.Point[i].x>frIn.width-1) out.Point[i].x=frIn.width-1;
		if(out.Point[i].y<0) out.Point[i].y=0;
		if(out.Point[i].y>frIn.height-1) out.Point[i].y=frIn.height-1;

	}
	out.nPoints=circleLength;
}



void print(CFrame &fr){

	for(int i=0;i<out.nPoints;i++){
		if(out.Point[i].x>=fr.width) continue;
		if(out.Point[i].x<0) continue;
		if(out.Point[i].y>=fr.height) continue;
		if(out.Point[i].y<0) continue;

		*(fr.pb(int(out.Point[i].x),int(out.Point[i].y)))=0;
		*(fr.pg(int(out.Point[i].x),int(out.Point[i].y)))=255;
		*(fr.pr(int(out.Point[i].x),int(out.Point[i].y)))=0;

	
	}

}


void print(CFrame &fr,CPoint &p,double peakPos[2]){
	double ds=2*PI/circleLength;
	int x;
	int y;
	for(int i=0;i<2;i++){
		x=int(p.x+radius*cos(peakPos[i]*ds));
		y=int(p.y+radius*sin(peakPos[i]*ds));

		if(x>=fr.width) continue;
		if(x<0) continue;
		if(y>=fr.height) continue;
		if(y<0) continue;
		if(fr.channel==3){
		*fr.pb(x,y)=0;
		*fr.pg(x,y)=0;
		*fr.pr(x,y)=255;}
		else if(fr.channel==1){
		*fr.pb(x,y)=255;
		} 


	
	}

}



public:
CProfile out;
int outFlag;
int circleLength;
double radius;
CFrame frIn;


};