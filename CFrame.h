#pragma once

#include "lib.h"

class CFrame
{
public:

	char * p;
	int width;
	int height;
	int channel;
	int widthStep;
	IplImage* f;

	CFrame(){};
	CFrame(IplImage* pFrame){
	set(pFrame);
	}
void set(IplImage* pFrame){
		p=pFrame->imageData;  //指向排列的图像数据的首地址
		width=pFrame->width;
		height=pFrame->height;
		channel=pFrame->nChannels;
		widthStep=pFrame->widthStep;
		f=pFrame;
		
	}//for read
	unsigned char r(int x,int y){return (unsigned char)*(p+y*widthStep+x*channel+2);};//原点是顶左
	unsigned char g(int x,int y){return (unsigned char)*(p+y*widthStep+x*channel+1);};//y、x分别代表行和列,返回xy那个点G的值
	unsigned char b(int x,int y){return (unsigned char)*(p+y*widthStep+x*channel+0);};//字节排列顺序为BGR

	double bf(double x,double y){return field(x,y,0);};//原点是顶左
	double gf(double x,double y){return field(x,y,1);};//y、x分别代表行和列,返回xy那个点G的值
	double rf(double x,double y){return field(x,y,2);};//字节排列顺序为BGR


	double field(double x,double y,int chl){
	int ix=int(x);
	if(ix>=width-1)ix=width-2;
	if(ix<0)ix=0;
	int iy=int(y);
	if(iy>=height-1)iy=height-2;
	if(iy<0)iy=0;

	double dx=x-int(x);
	double dy=y-int(y);	
	double xy  =double(*(pb(ix  ,iy  )+chl));	
	double x1y =double(*(pb(ix+1,iy  )+chl));	
	double xy1 =double(*(pb(ix  ,iy+1)+chl));	
	double x1y1=double(*(pb(ix+1,iy+1)+chl));	
	
	return xy*(1-dx)*(1-dy)+x1y*dx*(1-dy)+xy1*(1-dx)*dy+x1y1*dx*dy;

	}
	void  grad(double x,double y,int chl,double &gx,double &gy){
	int ix=int(x);
	if(ix>=width-1)ix=width-2;
	if(ix<0)ix=0;
	int iy=int(y);
	if(iy>=height-1)iy=height-2;
	if(iy<0)iy=0;

	double dx=x-int(x);
	double dy=y-int(y);	
	double xy  =double(*(pb(ix  ,iy  )+chl));	
	double x1y =double(*(pb(ix+1,iy  )+chl));	
	double xy1 =double(*(pb(ix  ,iy+1)+chl));	
	double x1y1=double(*(pb(ix+1,iy+1)+chl));	
	gx=xy*(-1)*(1-dy)+x1y*1*(1-dy)+xy1*(-1)*dy+x1y1*1*dy;
	gy=xy*(1-dx)*(-1)+x1y*dx*(-1)+xy1*(1-dx)*1+x1y1*dx*1;
	
	}

	unsigned char* pr(int x,int y){return (unsigned char*)(p+y*widthStep+x*channel+2);};//for write
	unsigned char* pg(int x,int y){return (unsigned char*)(p+y*widthStep+x*channel+1);};//返回xy那个点，G的值的指针，然后就可以修改
	unsigned char* pb(int x,int y){return (unsigned char*)(p+y*widthStep+x*channel+0);};

	unsigned char* pbyte(int x,int y){return pb(x,y);};
	float * pf(int x,int y){return (float*)(pbyte(x,y));};
};