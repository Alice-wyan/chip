#pragma once

#include "lib.h"

#define MAXFLOAT 1000

class CEdgeToFloat
{
public:
CEdgeToFloat(){
	ini();

};


~CEdgeToFloat(){
	clear(); 

};

void ini(){
	cmFlag=0;
	outFlag=0;
	total=0;
}

void setup(IplImage *in,IplImage *edgef,CProfile *edg){
	frIn.set(in);
	frEdge.set(edgef);

	pn.setup(in,10);
	ap.setup(in,&pn.out,2);

	edge=edg;

	if(cmFlag==1) {cvReleaseImage(&curveMap);cmFlag=0;}

	curveMap = cvCreateImage( cvGetSize(in), 8, 1 );
	cmFlag=1;

	cvZero(curveMap);

	setupOut();

}



void setupOut(){
	if(outFlag==0){
	out=new CProfile[MAXFLOAT];
	outFlag=1;
	}
}

void clear(){
	if(cmFlag==1) {cvReleaseImage(&curveMap);cmFlag=0;}

	if(outFlag==1){
	for(int i=0;i<total;i++){
		out[i].clear();
	}
	delete [] out;
	outFlag=0;
	total=0;
	}
}






void transform(){

	cvZero(frEdge.f);

	double vx;
	double vy;
	CPoint base;

	for(int i=0;i<edge->nPoints;i++){

		if(frEdge.b(int(edge->Point[i].x),int(edge->Point[i].y))==255) continue;
	
		pn.transform(edge->Point[i]);

		ap.transform();
		if(ap.diffData[0][int(ap.peakPos)]<200) continue;

		base=edge->Point[i];
		vx=-pn.vy;
		vy=pn.vx;  

		for(int j=0;j<100;j++){
			pn.transform(base,1,vx,vy);
		    ap.transform();
		    ap.print(pn.frIn);
			vx=ap.peakInProfile.x-base.x;
			vy=ap.peakInProfile.y-base.y;
			base=ap.peakInProfile;

		}

		break;


		


		//pn.print(pn.frIn);

	
	}
}


void print(CFrame &fr){
}

void print(){

}

public:
CFrame frIn;
CFrame frEdge;
CProfile* edge;

IplImage * curveMap;
int cmFlag;

CPointToNorm pn;
CNormToPeak ap;


CProfile* out;
int outFlag;

int total;


};