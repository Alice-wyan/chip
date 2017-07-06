#pragma once

#include "lib.h"

class CPicToBin
{
public:
CPicToBin(){
	ini();

};


CPicToBin(IplImage *pFrame){
	ini();
	setup(pFrame);
};

~CPicToBin(){
	clear(); 
	if(outFlag==1)cvReleaseImage(&out);  
};

void ini(){
	memFlag=0;
	outFlag=0;
	in=NULL;
}

void setup(char *fn,CColorCut *ccO){
	if(memFlag==1) cvReleaseImage(&in);   
	in = cvLoadImage(fn, 1);
	memFlag=1;
	cc=ccO;
}

void clear(){
	if(memFlag==1) cvReleaseImage(&in);   
	memFlag=0;
}

void setup(IplImage *pFrame){
	in=pFrame;
	memFlag=0;
}


void transform(){
	if(outFlag==1)cvReleaseImage(&out); 
	if(in==NULL){
	cout<<"in==NULL"<<endl;
	exit(0);
	}
	out=cvCreateImage(cvGetSize(in),8,1);
	outFlag=1;

	frIn.set(in);
	frOut.set(out);
	BlackNum=0;
	WhiteNum=0;	
		for(int i=0;i<frIn.width;i++)
		{
			for(int j=0;j<frIn.height;j++)
			{
				if(cc->cut(frIn.pb(i,j))==1)
				{
					*frOut.pb(i,j)=255;
					//BinaryImageTemp.x = i; 
					//BinaryImageTemp.y = j; 
					//BinaryImage.push_back(BinaryImageTemp);
					WhiteNum++;
				}
				else
				{
					*frOut.pb(i,j)=0;
					BlackNum++;
				}
			}
		}


}


void print(){
cout<<"CPicToBin::WhiteNum="<<WhiteNum<<endl;
cout<<"CPicToBin::BlackNum="<<BlackNum<<endl;

}

public:
IplImage *in;
IplImage *out;
int memFlag;
int outFlag;
CFrame frIn;
CFrame frOut;
CColorCut *cc;

int WhiteNum;
int BlackNum;

};