#pragma once

#include "lib.h"

#define MAXORM 100

class CCircleToPeak
{
public:
CCircleToPeak(){
	ini();

};


~CCircleToPeak(){
	clear(); 

};

void ini(){


vc=new CfChannel[MAXORM];
diffData=new double*[3];
diffData[0]=new double[MAXORM];
diffData[1]=new double[MAXORM];
diffData[2]=new double[MAXORM];

}

void setup(IplImage *pFrame,CProfile *inO,int l){
	setupFrame(pFrame);
	setupProfile(inO);
	setupLength(2);
}



void setupFrame(IplImage *pFrame){
	frIn.set(pFrame);
}


void setupProfile(CProfile *inO){
in=inO;
}

void setupLength(int l){
diffLength=l;
}

void setupChannel(){
	for(int i=0;i<in->nPoints;i++){
		vc[i].b=frIn.bf(in->Point[i].x,in->Point[i].y);
		vc[i].g=frIn.gf(in->Point[i].x,in->Point[i].y);
		vc[i].r=frIn.rf(in->Point[i].x,in->Point[i].y);
	}
}

void diff(){
	for(int i=0;i<in->nPoints;i++){
		diffData[0][i]=vc[(i+diffLength)%in->nPoints].b-vc[(i-diffLength+in->nPoints)%in->nPoints].b;
		diffData[1][i]=vc[(i+diffLength)%in->nPoints].g-vc[(i-diffLength+in->nPoints)%in->nPoints].g;
		diffData[2][i]=vc[(i+diffLength)%in->nPoints].r-vc[(i-diffLength+in->nPoints)%in->nPoints].r;
	}
}



void findPeak(){
	height[0]=0;
	height[1]=0;
	for(int i=0;i<in->nPoints;i++){
		if(diffData[0][i]>height[1]){
		height[1]=diffData[0][i];
		maxpos=i;
		}
		if(diffData[0][i]<height[0]){
		height[0]=diffData[0][i];
		minpos=i;
		}		
	}
}


void clear(){

	delete [] vc;
	delete [] diffData[0];
	delete [] diffData[1];
	delete [] diffData[2];
	delete [] diffData;

}



void setStartStop(int min,int max){

int istart=(max+min)/2;
int istop=(istart+in->nPoints/2)%in->nPoints;
		start.b=vc[istart].b;
		start.g=vc[istart].g;
		start.r=vc[istart].r;
		stop.b=vc[istop].b;
		stop.g=vc[istop].g;
		stop.r=vc[istop].r;
}


void uniformDiff(){

		if(start.b>stop.b){
			for(int i=0;i<in->nPoints;i++){
			diffData[0][i]*=-1;
			}
		}

		if(start.g>stop.g){
			for(int i=0;i<in->nPoints;i++){
			diffData[1][i]*=-1;
			}
		}

		if(start.r>stop.r){
			for(int i=0;i<in->nPoints;i++){
			diffData[2][i]*=-1;
			}
		}

		for(int i=0;i<in->nPoints;i++){
		diffData[0][i]+=diffData[1][i];
		diffData[0][i]+=diffData[2][i];
		}


}


void transform(){

		setupChannel();
		diff();
		setStartStop(0,0);
		uniformDiff();
		findPeak();
		setStartStop(minpos,maxpos);
		

       double wei;
       int nj;
	   for(int i=0;i<2;i++){
		peakPos[i]=0;
		wei=0;
		int spos;
		if(i==0)spos=minpos;
		else spos=maxpos;

		for(int j=spos-diffLength;j<=spos+diffLength;j++){
			nj=(j+in->nPoints)%in->nPoints;

		    peakPos[i]+=diffData[0][nj]*j;
		    wei+=diffData[0][nj];
		}
		peakPos[i]=peakPos[i]/wei;

		if(peakPos[i]<0)peakPos[i]+=(double)in->nPoints;
		if(peakPos[i]>=in->nPoints)peakPos[i]-=(double)in->nPoints;

	   }

}


void print(){


cout<<" CCircleToPeak::height "<<height[0]<<" "<<height[1]<<endl;
cout<<"CCircleToPeak::peakPos "<<peakPos[0]<<" "<<peakPos[1]<<endl;


}



void print(TH1D *h,TH1D *raw,int chl=0){
	if(h==NULL){		
		return;}

	for(int i=0;i<in->nPoints;i++){
		//cout<<diffData[chl][i]<<endl;
		h->SetBinContent(i+1,diffData[chl][i]);

	}
	for(int i=0;i<in->nPoints;i++){
		//cout<<vc[i].b<<endl;
		raw->SetBinContent(i+1,vc[i].b);
	}

}

double colorDiff(){
	return fabs(start.b-stop.b)+fabs(start.g-stop.g)+fabs(start.r-stop.r);
}

public:
CProfile *in;

CFrame frIn;
CfChannel * vc;
double ** diffData;

int diffLength;

////////////

CfChannel start;
CfChannel stop;




double height[2];
double peakPos[2];
int minpos;
int maxpos;

};