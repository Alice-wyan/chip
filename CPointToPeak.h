#pragma once

#include "lib.h"

#define MAXORM 100

class CPointToPeak
{



public:

	CPointToPeak(){
	peakToOriginMax=3;
	peakDisMax=2;
	peakMin=30;
	diffLengthL=1;
	diffLengthR=1;

	peakLengthL=1;
	peakLengthR=1;

	}

//double diffX(int x,int y){
//	int scale=0;
//	int sum=0;
//	if(x<diffLength){
//	scale=diffLength+x;
//	p1 = cvGet2D(frIn.f,y, 0);
//	p2 = cvGet2D(frIn.f,y,x+diffLength);
//	sum=abs(p1.val[0]-p2.val[0])+abs(p1.val[1]-p2.val[1])+abs(p1.val[2]-p2.val[2]);
//	return double(sum)/double(scale);
//	}
//	if(x+diffLength>=frIn.width){
//    scale=diffLength+frIn.width-1-x;
//
//	p1 = cvGet2D(frIn.f,y, x-diffLength);
//	p2 = cvGet2D(frIn.f,y, frIn.width-1);
//	sum=abs(p1.val[0]-p2.val[0])+abs(p1.val[1]-p2.val[1])+abs(p1.val[2]-p2.val[2]);
//	return double(sum)/double(scale);
//	}
//	scale=diffLength*2;
//	p1 = cvGet2D(frIn.f,y, x-diffLength);
//	p2 = cvGet2D(frIn.f,y, x+diffLength);
//	sum=abs(p1.val[0]-p2.val[0])+abs(p1.val[1]-p2.val[1])+abs(p1.val[2]-p2.val[2]);
//	//sum=(p1.val[0]-p2.val[0])+(p1.val[1]-p2.val[1])+(p1.val[2]-p2.val[2]);
//	return double(sum)/double(scale);
//}
//
//double diffY(int x,int y){
//	int scale=0;
//	int sum=0;
//	if(y<diffLength){
//	scale=diffLength+y;
//	p1 = cvGet2D(frIn.f,0, x);
//	p2 = cvGet2D(frIn.f,y+diffLength,x);
//	sum=abs(p1.val[0]-p2.val[0])+abs(p1.val[1]-p2.val[1])+abs(p1.val[2]-p2.val[2]);;
//	return double(sum)/double(scale);
//	}
//	if(y+diffLength>=frIn.height){
//    scale=diffLength+frIn.height-1-y;
//
//	p1 = cvGet2D(frIn.f,y-diffLength, x);
//	p2 = cvGet2D(frIn.f,frIn.height-1, x);
//	sum=abs(p1.val[0]-p2.val[0])+abs(p1.val[1]-p2.val[1])+abs(p1.val[2]-p2.val[2]);;
//	return double(sum)/double(scale);
//
//	}
//	scale=diffLength*2;
//	p1 = cvGet2D(frIn.f,y-diffLength, x);
//	p2 = cvGet2D(frIn.f,y+diffLength, x);
//	sum=abs(p1.val[0]-p2.val[0])+abs(p1.val[1]-p2.val[1])+abs(p1.val[2]-p2.val[2]);;
//	return double(sum)/double(scale);
//
//}




double diffX(int x,int y){
	int scale=0;
	int sum=0;
	if(x<diffLengthL){
	scale=diffLengthR+x;
	p1 = cvGet2D(frIn.f,y, 0);
	p2 = cvGet2D(frIn.f,y,x+diffLengthR);
	sum=abs(p1.val[0]-p2.val[0])+abs(p1.val[1]-p2.val[1])+abs(p1.val[2]-p2.val[2]);
	return double(sum)/double(scale);
	}
	if(x+diffLengthR>=frIn.width){
    scale=diffLengthL+frIn.width-1-x;

	p1 = cvGet2D(frIn.f,y, x-diffLengthL);
	p2 = cvGet2D(frIn.f,y, frIn.width-1);
	sum=abs(p1.val[0]-p2.val[0])+abs(p1.val[1]-p2.val[1])+abs(p1.val[2]-p2.val[2]);
	return double(sum)/double(scale);
	}
	scale=diffLengthL+diffLengthR;
	p1 = cvGet2D(frIn.f,y, x-diffLengthL);
	p2 = cvGet2D(frIn.f,y, x+diffLengthR);
	sum=abs(p1.val[0]-p2.val[0])+abs(p1.val[1]-p2.val[1])+abs(p1.val[2]-p2.val[2]);
	//sum=(p1.val[0]-p2.val[0])+(p1.val[1]-p2.val[1])+(p1.val[2]-p2.val[2]);
	return double(sum)/double(scale);
}

double diffY(int x,int y){
	int scale=0;
	int sum=0;
	if(y<diffLengthL){
	scale=diffLengthR+y;
	p1 = cvGet2D(frIn.f,0, x);
	p2 = cvGet2D(frIn.f,y+diffLengthR,x);
	sum=abs(p1.val[0]-p2.val[0])+abs(p1.val[1]-p2.val[1])+abs(p1.val[2]-p2.val[2]);;
	return double(sum)/double(scale);
	}
	if(y+diffLengthR>=frIn.height){
    scale=diffLengthL+frIn.height-1-y;

	p1 = cvGet2D(frIn.f,y-diffLengthL, x);
	p2 = cvGet2D(frIn.f,frIn.height-1, x);
	sum=abs(p1.val[0]-p2.val[0])+abs(p1.val[1]-p2.val[1])+abs(p1.val[2]-p2.val[2]);;
	return double(sum)/double(scale);

	}
	scale=diffLengthL+diffLengthR;
	p1 = cvGet2D(frIn.f,y-diffLengthL, x);
	p2 = cvGet2D(frIn.f,y+diffLengthR, x);
	sum=abs(p1.val[0]-p2.val[0])+abs(p1.val[1]-p2.val[1])+abs(p1.val[2]-p2.val[2]);;
	return double(sum)/double(scale);

}

double peakX(int x,int y,double & pos){

	int start=x-peakToOriginMax;
	int stop =x+peakToOriginMax;

	if(x<peakToOriginMax) start=0;
	else if(x+peakToOriginMax>=frIn.width) stop=frIn.width-1;

	double max=diffX(start,y);
	double p,pl,pr;
	int pm=start;
	for(int i=start+1;i<=stop;i++){
		p=diffX(i,y);
		if(p>max){
		max=p;
		pm=i;
		}
	}

	if(pm-peakLengthL<=start||pm+peakLengthL>=stop) return -1;
	pl=max*pm;
	pr=max;
	for(int i=1;i<=peakLengthL;i++){
		pos=diffX(pm-i,y);
	pl+=pos*(pm-i);
	pr+=pos;
	}

	for(int i=1;i<=peakLengthR;i++){
		pos=diffX(pm+i,y);
	pl+=pos*(pm+i);
	pr+=pos;
	}
	pos=pl/pr;
	return max;
}


double peakY(int x,int y,double & pos){

	int start=y-peakToOriginMax;
	int stop =y+peakToOriginMax;

	if(y<peakToOriginMax) start=0;
	else if(y+peakToOriginMax>=frIn.height) stop=frIn.height-1;

	double max=diffX(x,start);
	double p,pl,pr;
	int pm=start;
	for(int i=start+1;i<=stop;i++){
		p=diffY(x,i);
		if(p>max){
		max=p;
		pm=i;
		}
	}

	if(pm-peakLengthL<=start||pm+peakLengthL>=stop) return -1;
	pl=max*pm;
	pr=max;
	for(int i=1;i<=peakLengthL;i++){
		pos=diffY(x,pm-i);
	pl+=pos*(pm-i);
	pr+=pos;
	}

	for(int i=1;i<=peakLengthR;i++){
		pos=diffY(x,pm+i);
	pl+=pos*(pm+i);
	pr+=pos;
	}
	pos=pl/pr;
	return max;
}


int connect(int x1,int y1,int x2,int y2){
	int dx=x1-x2;
	int dy=y1-y2;

	double pk,ps;
	double nearp;
	int d;
	int sig;
	if(abs(dx)>abs(dy)){
		d=abs(dx);
		if(dx>0) sig=-1;else sig=1;

		peak.clear();
		posx.clear();
		posy.clear();
		nearp=y1;
		for(int i=0;i<=d;i++){
			pk=peakY(x1+i*sig,int(nearp),nearp);

			if(pk<0) return -1;
			peak.push_back(pk);
			posx.push_back(x1+i*sig);
			posy.push_back(nearp);
		}

	}else{
		d=abs(dy);
		if(dy>0) sig=-1;else sig=1;	
	    if(dy==0) return -1;

		peak.clear();
		posx.clear();
		posy.clear();
		nearp=x1;
		for(int i=0;i<=d;i++){
			pk=peakX(int(nearp),y1+i*sig,nearp);
			//cout<<nearp<<endl;
			if(pk<0) return -1;
			peak.push_back(pk);
			posx.push_back(nearp);
			posy.push_back(y1+i*sig);
		}

	}


	return d;

}


int tryConnect(int x1,int y1,int x2,int y2){
	int dx=x1-x2;
	int dy=y1-y2;

	double pk,ps;
	double nearp;
	int d;
	int sig;
	if(abs(dx)>abs(dy)){
		d=abs(dx);
		if(dx>0) sig=-1;else sig=1;

		//peak.clear();
		//posx.clear();
		//posy.clear();
		nearp=y1;
		for(int i=0;i<=d;i++){
			ps=nearp;
			pk=peakY(x1+i*sig,int(nearp),nearp);
			if(fabs(ps-nearp)>peakDisMax) return -1;
			if(pk<peakMin) return -1;
			//peak.push_back(pk);
			//posx.push_back(x1+i*sig);
			//posy.push_back(nearp);
		}
		if(fabs(nearp-y2)>peakDisMax)return -1;
	}else{
		d=abs(dy);
		if(dy>0) sig=-1;else sig=1;	
	    if(dy==0) return -1;

		//peak.clear();
		//posx.clear();
		//posy.clear();
		nearp=x1;
		for(int i=0;i<=d;i++){
			ps=nearp;
			pk=peakX(int(nearp),y1+i*sig,nearp);
			//cout<<nearp<<endl;
			if(fabs(ps-nearp)>peakDisMax) return -1;
			if(pk<peakMin) return -1;
			//peak.push_back(pk);
			//posx.push_back(nearp);
			//posy.push_back(y1+i*sig);
		}
		if(fabs(nearp-x2)>peakDisMax)return -1;

	}

	return d;


}


void print(){

	for(unsigned int i=0;i<peak.size();i++){
	*(frIn.pb(int(posx[i]),int(posy[i])))=0;
	*(frIn.pg(int(posx[i]),int(posy[i])))=255;
	*(frIn.pr(int(posx[i]),int(posy[i])))=0;
	}

}

private:
CvScalar p1;
CvScalar p2;
public:

CFrame frIn;
int peakToOriginMax;
double peakDisMax;
double peakMin;

int diffLengthR;
int diffLengthL;

int peakLengthR;
int peakLengthL;


vector<double> peak;
vector<double> posx;
vector<double> posy;

};