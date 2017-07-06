#pragma once

#include "lib.h"

class CPoint
{
	public:
		double x;
		double y;

		void print(){
		cout<<"CPoint:: x="<<x<<" y="<<y<<endl;
		}
};

class CfChannel
{
	public:
		double b;
		double g;
		double r;
//	double alpha;
};

class CProfile
{
	public:

	CProfile(){
		memFlag=0;
		nPoints=0;

	};	
	CProfile(int maxMem){
		memFlag=0;
		nPoints=0;
		setup(maxMem);
	};	

	~CProfile(){
		clearMem();
	};

	void setup(int max){
	clearMem();
	Point=new CPoint[max];
	memFlag=1;	
	}

    void clear()
	{
		nPoints=0;		
	}

	void print(){
	cout<<"CProfile::nPoints="<<nPoints<<endl;
	for(int i=0;i<nPoints;i++){
		Point[i].print();
	}
	}

    void clearMem()
	{
	if(memFlag==1) delete []Point;
	nPoints=0;
	memFlag=0;
	}

double PointDistance(CPoint &a,CPoint &b){

return sqrt(pow(a.x-b.x,2)+pow(a.y-b.y,2));

}

double PointDistance(int i,int j){

return PointDistance(Point[(i+nPoints)%nPoints],Point[(j+nPoints)%nPoints]);

}

double BeforeAngle(int i,int dis){
	return atan2(Point[(i+nPoints)%nPoints].y-Point[(i-dis+nPoints)%nPoints].y,Point[(i+nPoints)%nPoints].x-Point[(i-dis+nPoints)%nPoints].x);
}
double AfterAngle(int i,int dis){
	return atan2(Point[(i+dis+nPoints)%nPoints].y-Point[(i+nPoints)%nPoints].y,Point[(i+dis+nPoints)%nPoints].x-Point[(i+nPoints)%nPoints].x);
}

	CPoint *Point;
	int nPoints;
	int memFlag;

};

