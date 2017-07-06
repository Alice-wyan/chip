#pragma once

#include "lib.h"

class CBinToProfile
{
public:
CBinToProfile(){
	ini();

};
CBinToProfile(char *fn){
	ini();

	setup(fn);
};

CBinToProfile(IplImage *pFrame){
	ini();
	setup(pFrame);
};

~CBinToProfile(){
	clear(); 
	clearOutMem();
	//cvReleaseMemStorage(&storage);
};

void ini(){
	storageFlag=0;
	maxPoint=5000;
    //storage = cvCreateMemStorage(0); 	
	memFlag=0;
	in=NULL;
	transFlag=0;
}

void setup(char *fn){
	if(memFlag==1) cvReleaseImage(&in);   
	in = cvLoadImage(fn, 1);
	memFlag=1;
}

void setupMaxPoint(int mp){
	maxPoint=mp;
}

void clear(){
	if(memFlag==1) cvReleaseImage(&in);   
	memFlag=0;
}

void clearOutMem(){
if(transFlag==1){
	for(int i=0;i<totals;i++){
		out[i].clearMem();
	}
	delete [] out;
}
}


void setup(IplImage *pFrame){
	in=pFrame;
	memFlag=0;
}


void transform(){
	if(in==NULL){
	cout<<"in==NULL"<<endl;
	exit(0);
	}

clearOutMem();
storage = cvCreateMemStorage(0);



totals = cvFindContours(in,storage,&contours,sizeof(CvContour),CV_RETR_CCOMP,CV_CHAIN_APPROX_NONE,cvPoint(0,0));
storageFlag=1;



out=new CProfile[totals];

    int count=0;
    int i=0;  
    for(;contours != 0; contours = contours -> h_next)  /// 这样可以访问每一个轮廓  ====横向轮廓  
    {  
	    
        out[count].setup(contours->total);
		out[count].nPoints=contours->total;
		
        for(i = 0; i < contours -> total; i++)    // 提取一个轮廓的所有坐标点  
        {  
            CvPoint *pt = (CvPoint*) cvGetSeqElem(contours, i);   // 得到一个轮廓中一个点的函数cvGetSeqElem  
			out[count].Point[i].x=(double)pt->x;
			out[count].Point[i].y=(double)pt->y;

        } 
		count++;

		

        CvSeq *InterCon = contours->v_next;     // 访问每个轮廓的纵向轮廓  
        for(; InterCon != 0; InterCon = InterCon ->h_next)  
        {  
        out[count].setup(InterCon->total);
		out[count].nPoints=InterCon->total;

            for(i = 0; i < InterCon->total; i++ )  
            {  
                CvPoint *pt = (CvPoint*)cvGetSeqElem(InterCon, i);  
			out[count].Point[i].x=(double)pt->x;
			out[count].Point[i].y=(double)pt->y;

            }
			count++;

        }  
    }  

cvReleaseMemStorage(&storage);

}


void print(IplImage *can){
	CFrame fr;
	fr.set(can);
	for(int i=0;i<totals;i++){
		cout<<"CBinToProfile::out[i].nPoints="<<out[i].nPoints<<endl;
		for(int j=0;j<out[i].nPoints;j++){
			
			(*fr.pb(int(out[i].Point[j].x),int(out[i].Point[j].y)))=255;
		}
	
	}

	cout<<"CBinToProfile::totals="<<totals<<endl;

}

public:

IplImage *in;
CProfile* out;

int memFlag;
int maxPoint;
int transFlag;
int storageFlag;

CFrame frIn;
CvMemStorage *storage;
CvSeq *contours;
int totals;

};