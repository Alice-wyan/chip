#pragma once

#include "lib.h"

#define MAXFLOAT 1000

class CEdgeToCircleChain
{
public:
CEdgeToCircleChain(){
	ini();

};


~CEdgeToCircleChain(){
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
	pc.setup(in);
	cp.setup(in,&pc.out,2);

	edge=edg;

	if(cmFlag==1) {cvReleaseImage(&curveMap);cmFlag=0;}

	curveMap = cvCreateImage( cvGetSize(in), 8, 1 );
	cmFlag=1;
    frCm.set(curveMap);
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

//	cvZero(frEdge.f);
	cvZero(curveMap);

	CPoint base;
	double pdir=0;
	double dir0=0;
	double dir1=0;
	double dir=0;
	double ds=2*PI/pc.circleLength;

	for(int i=0;i<edge->nPoints;i++){
		base=edge->Point[i];
		pdir=0;



		for(int j=0;j<5;j++){

		pc.transform(base);
		cp.transform();
		if(cp.colorDiff()<100) break;
			pdir=pdir+0.5*pc.circleLength;
			if(pdir>=pc.circleLength)dir0-=pc.circleLength;

			dir0=fabs(pdir-cp.peakPos[0]);
			dir1=fabs(pdir-cp.peakPos[1]);

			if(dir0>dir1){
			pdir=cp.peakPos[0];
			}else{
			pdir=cp.peakPos[1];			
			}
			base.x+=pc.radius*cos(pdir*ds);
			base.y+=pc.radius*sin(pdir*ds);

		    pc.print(frCm,base,cp.peakPos);
		}
		

		//*frIn.pr(int(base.x),int(base.y))=255;

		break;
	
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
CFrame frCm;

CPointToNorm pn;
CPointToCircle pc;
CCircleToPeak cp;
CNormToPeak ap;


CProfile* out;
int outFlag;

int total;


};