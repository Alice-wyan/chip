#pragma once

#include "lib.h"

#define MAXNORM 50
#define USEDNORM 25

class CIntProfileToFloat
{
public:
CIntProfileToFloat(){
	setFlag=0;
	outMemFlag=0;
	interval=4;

};
~CIntProfileToFloat(){
	if(setFlag==1) {
		out->clearMem();
		delete out;}
	setFlag=0;
};




void CalcClosePoint(CPoint & NormalPoint){
		
		ClosePoint.x = (int)NormalPoint.x;
		ClosePoint.y = (int)NormalPoint.y;
		
		decimal.x = NormalPoint.x - ClosePoint.x;
		decimal.y = NormalPoint.y - ClosePoint.y;

}

double BilinearInterpolation(double u,double v,double *pixel)
{	
	return (1-u)*(1-v)*pixel[0] + u*(1-v)*pixel[1] + (1-u)*v*pixel[2] + u*v*pixel[3];
}

CfChannel CalcInterpolationChannel()
{
	CfChannel Channel;


		ClosePointB[0] = (double)fr.b(int(ClosePoint.x),int(ClosePoint.y));
		ClosePointG[0] = (double)fr.g(int(ClosePoint.x),int(ClosePoint.y));
		ClosePointR[0] = (double)fr.r(int(ClosePoint.x),int(ClosePoint.y));

		ClosePointB[1] = (double)fr.b(int(ClosePoint.x)+1,int(ClosePoint.y));
		ClosePointG[1] = (double)fr.g(int(ClosePoint.x)+1,int(ClosePoint.y));
		ClosePointR[1] = (double)fr.r(int(ClosePoint.x)+1,int(ClosePoint.y));

		ClosePointB[2] = (double)fr.b(int(ClosePoint.x)+1,int(ClosePoint.y)+1);
		ClosePointG[2] = (double)fr.g(int(ClosePoint.x)+1,int(ClosePoint.y)+1);
		ClosePointR[2] = (double)fr.r(int(ClosePoint.x)+1,int(ClosePoint.y)+1);

		ClosePointB[3] = (double)fr.b(int(ClosePoint.x),int(ClosePoint.y)+1);
		ClosePointG[3] = (double)fr.g(int(ClosePoint.x),int(ClosePoint.y)+1);
		ClosePointR[3] = (double)fr.r(int(ClosePoint.x),int(ClosePoint.y)+1);
	
	Channel.b = BilinearInterpolation(decimal.x,decimal.y,ClosePointB);
	Channel.g = BilinearInterpolation(decimal.x,decimal.y,ClosePointG);
	Channel.r = BilinearInterpolation(decimal.x,decimal.y,ClosePointR);
	
	return Channel;
}

double CalcFloatPointPixel(CPoint *Channel)
{
	double TotalLength=0;
	double TotalPixel=0;
	
	for(int i=0;i<normal->nPoints-2;i++)
	{
		FirstDerivative[i].x = (double)i+1;
		FirstDerivative[i].y = fabs(Channel[i+2].y - Channel[i].y);	
	}

	for(int i=0;i<normal->nPoints-2;i++)
	{
		TotalLength += FirstDerivative[i].x * FirstDerivative[i].y;
		TotalPixel += FirstDerivative[i].y;		
	}
	double ap = (double)TotalLength/TotalPixel;
	
	return ap;
}


void CalcFloatEdgePoint()
{
	CfChannel Channel;

	double b_point = 0;
	double g_point = 0;
	double r_point = 0;
	CPoint decimal={0,0};

	for(int i=0;i<normal->nPoints;i++)
	{
		CalcClosePoint(normal->Point[i]);
		Channel = CalcInterpolationChannel();
		ChannelB[i].x = (double)i;
		ChannelB[i].y = Channel.b;
		ChannelG[i].x = (double)i;
		ChannelG[i].y = Channel.g;
		ChannelR[i].x = (double)i;
		ChannelR[i].y = Channel.r;
	}

	b_point = CalcFloatPointPixel(ChannelB);
	g_point = CalcFloatPointPixel(ChannelG);
	r_point = CalcFloatPointPixel(ChannelR);
	avePoint = (b_point + g_point + r_point)/3 + 0;

}


void CalcFloatCoordinate(CPoint &Coordinate)
{


	int LeftPoint = (int) avePoint;

	if(avePoint - LeftPoint == 0)
	{
		Coordinate=normal->Point[LeftPoint];
	}

	int RightPoint = LeftPoint + 1;
	double delta = avePoint - LeftPoint;

	double hypotenuse = sqrt(pow(normal->Point[RightPoint].x - normal->Point[LeftPoint].x,2)+pow(normal->Point[RightPoint].y - normal->Point[LeftPoint].y,2));
	double cosine = (normal->Point[RightPoint].x - normal->Point[LeftPoint].x)/hypotenuse;
	double sine = (normal->Point[RightPoint].y - normal->Point[LeftPoint].y)/hypotenuse;

	Coordinate.x = normal->Point[LeftPoint].x + delta*cosine;
	Coordinate.y = normal->Point[LeftPoint].y + delta*sine;
	
}

void setup(CProfile *inO,int index,IplImage *srcO,CCubicBSpline *cbO){
	in=&(inO[index]);
	src=srcO;
	cb=cbO;
	fr.set(src);
	normal=&(cb->norm);


	if(outMemFlag==1){
	out->clearMem();
	out->setup(in->nPoints);
	}else{
	out=new CProfile(in->nPoints);
	outMemFlag=1;
	}

	setFlag=1;
}



void transform(){
	if(setFlag==0) {
	cout<<"setFlag==0"<<endl;
	exit(0);
	};

	for(int i=0;i<in->nPoints;i++){
	
	pos[0]=i-interval/2-interval;  //min interval is 2
	pos[1]=i-interval/2;
	pos[2]=i+interval/2;  //min interval is 2
	pos[3]=i+interval/2+interval;

	pos[0]=(pos[0]+in->nPoints)%(in->nPoints);
	pos[1]=(pos[1]+in->nPoints)%(in->nPoints);
	pos[2]=(pos[2]+in->nPoints)%(in->nPoints);
	pos[3]=(pos[3]+in->nPoints)%(in->nPoints);

	cb->fitPoint[0]=&(in->Point[pos[0]]);
	cb->fitPoint[1]=&(in->Point[pos[0]]);
	cb->fitPoint[2]=&(in->Point[pos[0]]);
	cb->fitPoint[3]=&(in->Point[pos[0]]);


	cb->fit();

	cb->setNormal(USEDNORM);

	CalcFloatEdgePoint();
	
	CalcFloatCoordinate(out->Point[i]);

	}
	out->nPoints=in->nPoints;


}


void print(){


}

public:

CProfile *in;
IplImage *src;
CCubicBSpline *cb;
CProfile *normal;

int setFlag;
int outMemFlag;

int interval;

CFrame fr;
CProfile *out;



/////////////
CPoint decimal;
CPoint ClosePoint;
double ClosePointB[4];
double ClosePointG[4];
double ClosePointR[4];
int  pos[4];

double avePoint;

CPoint	ChannelB[MAXNORM];
CPoint	ChannelG[MAXNORM];
CPoint	ChannelR[MAXNORM];
CPoint	FirstDerivative[MAXNORM];


};