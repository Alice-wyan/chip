#pragma once

#include "lib.h"

#define MAXORM 100

class CNormToPeak
{
public:
CNormToPeak(){
	ini();

};


~CNormToPeak(){
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

		if(i-diffLength<0||i+diffLength>=in->nPoints){
		diffData[0][i]=0;
		diffData[1][i]=0;
		diffData[2][i]=0;		
		}else{
		diffData[0][i]=vc[i+diffLength].b-vc[i-diffLength].b;
		diffData[1][i]=vc[i+diffLength].g-vc[i-diffLength].g;
		diffData[2][i]=vc[i+diffLength].r-vc[i-diffLength].r;



		}
	}

}


int firstPeak(double threshold,double &height, int chl=0, int checkLength=2){
	
	int found=1;

	for(int i=checkLength;i<in->nPoints-checkLength;i++){
		if(diffData[chl][i]<threshold) continue;
		found=1;
		for(int j=0;j<checkLength;j++){
			if(diffData[chl][i]<diffData[chl][i+j]){found=0;break; }
			if(diffData[chl][i]<diffData[chl][i-j]){found=0;break; }
		}
		if(found==1) {height=diffData[chl][i];return i;}
	}
	
	height=diffData[chl][in->nPoints/2];

	return in->nPoints/2;
}


int middlePeak(double threshold,double &height,int chl=0, int checkLength=2){
	
	int found=1;

	int last=-1;
	int first=-1;

	for(int i=in->nPoints/2;i<in->nPoints-checkLength;i++){
		if(diffData[chl][i]<threshold) continue;
		found=1;
		for(int j=0;j<checkLength;j++){
			if(diffData[chl][i]<diffData[chl][i+j]){found=0;break; }
			if(diffData[chl][i]<diffData[chl][i-j]){found=0;break; }
		}
		if(found==1) {last=i;break;}
	}

	for(int i=in->nPoints/2;i>=checkLength;i--){
		if(diffData[chl][i]<threshold) continue;
		found=1;
		for(int j=0;j<checkLength;j++){
			if(diffData[chl][i]<diffData[chl][i+j]){found=0;break; }
			if(diffData[chl][i]<diffData[chl][i-j]){found=0;break; }
		}
		if(found==1) {first=i;break;}
	}

	if(first==-1&&last==-1){
	height=diffData[chl][in->nPoints/2];
    return in->nPoints/2;
	}else if(abs(last-in->nPoints/2)>abs(first-in->nPoints/2)){
	height=diffData[chl][first];
    return first;	
	}else{
	height=diffData[chl][last];
    return last;		
	}
	
}

int findPeak(int chl,double &height, int checkLength=2){
	
	double max=0;
	int found=in->nPoints/2;

	for(int i=checkLength;i<in->nPoints-checkLength;i++){
		if(diffData[chl][i]>max){
		max=diffData[chl][i];
		found=i;
		}
	}
	height=max;
	return found;
}


void clear(){

	delete [] vc;
	delete [] diffData[0];
	delete [] diffData[1];
	delete [] diffData[2];
	delete [] diffData;

}



void setStartStop(){

		start.b=(vc[0].b+vc[1].b)/2;
		start.g=(vc[0].g+vc[1].g)/2;
		start.r=(vc[0].r+vc[1].r)/2;

		stop.b=(vc[in->nPoints-1].b+vc[in->nPoints-2].b)/2;
		stop.g=(vc[in->nPoints-1].g+vc[in->nPoints-2].g)/2;
		stop.r=(vc[in->nPoints-1].r+vc[in->nPoints-2].r)/2;

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


void transform(double threshold=-1){
	double height;
		setupChannel();
		diff();
		setStartStop();


		uniformDiff();

		if(threshold<0)
	   pos=findPeak(0,height);
		else
	   pos=firstPeak(threshold,height);

	   //cout<<"pos="<<pos<<endl;

		peakPos=0;
		wei=0;
		for(int j=pos-diffLength;j<=pos+diffLength;j++){
			if(j<0) continue;
			if(j>=in->nPoints) continue;

		    peakPos+=diffData[0][j]*j;
		    wei+=diffData[0][j];
		}
		peakPos=peakPos/wei;
		 //cout<<"peakPos="<<peakPos<<endl;
		pos=int(peakPos);
		wei=peakPos-pos;

		if(pos+1>=in->nPoints){
			peakInProfile=in->Point[in->nPoints-1];
		
		}else{

		peakInProfile.x=in->Point[pos].x*(1-wei)+in->Point[pos+1].x*(wei);
		peakInProfile.y=in->Point[pos].y*(1-wei)+in->Point[pos+1].y*(wei);
		}

}


void print(CFrame &fr){


	if(int(peakInProfile.x)<0)peakInProfile.x=0;
	if(int(peakInProfile.x)>=fr.width) peakInProfile.x=fr.width-1;
	if(int(peakInProfile.y)<0)peakInProfile.y=0;
	if(int(peakInProfile.y)>=fr.height) peakInProfile.y=fr.height-1;

		*(fr.pb(int(peakInProfile.x),int(peakInProfile.y)))=0;
		*(fr.pg(int(peakInProfile.x),int(peakInProfile.y)))=0;
		*(fr.pr(int(peakInProfile.x),int(peakInProfile.y)))=255;



}

void print(TH1D *h,TH1D *raw,int chl=0){
	if(h==NULL){		
		return;}

	cout<<"peakPos="<<peakPos<<endl;
	cout<<"peakInProfile.x="<<peakInProfile.x<<"  peakInProfile.y="<<peakInProfile.y<<endl;
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


double peakPos;
double wei;
int pos;
CPoint peakInProfile;


};