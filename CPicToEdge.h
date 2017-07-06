#pragma once

#include "lib.h"

class CPicToEdge
{
public:
CPicToEdge(){
	ini();

};


~CPicToEdge(){
	clear(); 

};

void ini(){
	memFlag=0;
	outFlag=0;
	in=NULL;
	dst=NULL;
}

void setup(IplImage *pFrame,int max=20000,int gap=4){
	setupFrame(pFrame);
	setupOut(max,gap);
}

void setupFile(char *fn){
	if(memFlag==1) cvReleaseImage(&in);   
	in = cvLoadImage(fn, 1);
	memFlag=1;
	if(dst!=NULL)  cvReleaseImage(&dst);   
	dst = cvCreateImage( cvGetSize(in), 8, 1 );
	frDst.set(dst);
}

void setupFrame(IplImage *pFrame){
	in=pFrame;
	memFlag=0;
	if(dst!=NULL)  cvReleaseImage(&dst);   
	dst = cvCreateImage( cvGetSize(in), 8, 1 );
	frDst.set(dst);
}


void setupOut(int max,int gap){
 maxOut=max;
 bin=gap;
 if(outFlag==1) delete []out.Point;
  out.Point=new CPoint[maxOut];
  outFlag=1;
 
}

void clear(){
	//cout<<"memFlag="<<memFlag<<endl;
	//cout<<"outFlag="<<outFlag<<endl;
	//cout<<"dst="<<int(dst)<<endl;
	if(memFlag==1) cvReleaseImage(&in);   
	memFlag=0;
	 if(outFlag==1) delete []out.Point;  
	outFlag=0;

	if(dst!=NULL)  cvReleaseImage(&dst);   
	dst=NULL;
}



void transform(){
	out.nPoints=0;
	cvCanny( in, dst, 50, 200, 3 );
	for(int i=0;i<frDst.width;i=i+bin){
		for(int j=0;j<frDst.height;j=j+bin){
			if(frDst.b(i,j)>0){
				out.Point[out.nPoints].x=double(i);
				out.Point[out.nPoints].y=double(j);
				out.nPoints++;}
	}}


}

void transformHit(IplImage * pf){//pe.transformHit(newEdge);
	posToIndex.set(pf);//for  read
	outHit.clear();
	CEdgeHit eh;//序列

//	cvCanny( in, dst, 50, 200, 3 );
	cvCanny( in, dst, 50, 200, 3 );//cvCanny边缘检测    
	for(int i=0;i<frDst.width;i=i+bin){//bin=1  frDst==dst即轮廓
		for(int j=0;j<frDst.height;j=j+bin){
			if(frDst.b(i,j)>0){//b通道>0说明什么
				eh.x=i;
				eh.y=j;
				((int*)posToIndex.pb(i,j))[0]=outHit.size();
				outHit.push_back(eh);//outHit里存储着CEdgeHit eh

			}
	}}


}

void transformHitMask(IplImage * pf,vector <CvPoint> & mask){
	posToIndex.set(pf);
	outHit.clear();
	CEdgeHit eh;

//	cvCanny( in, dst, 50, 200, 3 );
	cvCanny( in, dst, 50, 200, 3 );
	for(int i=0;i<frDst.width;i=i+bin){//bin=1
		for(int j=0;j<frDst.height;j=j+bin){


			if(frDst.b(i,j)>0){

			int find=0;
			for(int k=0;k<mask.size();k++){
				if(i==mask[k].x&&j==mask[k].y) {find=1;break;}
			}
			if(find==0) continue;


				eh.x=i;
				eh.y=j;
				((int*)posToIndex.pb(i,j))[0]=outHit.size();
				outHit.push_back(eh);

			}
	}}


}


void print(CFrame &fr){
	for(int i=0;i<out.nPoints;i++){
	*(fr.pb(int(out.Point[i].x),int(out.Point[i].y)))=0;
	*(fr.pg(int(out.Point[i].x),int(out.Point[i].y)))=255;
	*(fr.pr(int(out.Point[i].x),int(out.Point[i].y)))=0;
	}
}

void print(){
cout<<"CPicToEdge::out.nPoints="<<out.nPoints<<endl;
}

public:
IplImage *in;
IplImage * dst;
CProfile out;
int maxOut;
int bin;

CFrame frDst;
vector<CEdgeHit> outHit;
CFrame posToIndex;

private:
int memFlag;
int outFlag;



};