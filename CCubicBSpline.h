#pragma once

#include "lib.h"

class CCubicBSpline
{
public:

	CCubicBSpline(){
		norm.setup(50);
	}
	~CCubicBSpline(){
	norm.clearMem();

	}


void setNormal(int maxPoint)
{

	CPoint normal = {0,0};
	for(int i=0;i<=maxPoint-1;i++)
	{			
		if(b[1]==0&&b[2]==0&&b[3]==0)//垂直于X轴
		{
			normal.x = a[0] + a[1]*(0.5) + a[2]*0.25 + a[3]*0.125;														 
			normal.y = b[0] + (i-maxPoint/2.0);	
		}				
		else if(a[1]==0&&a[2]==0&&a[3]==0)//垂直于Y轴
		{
			normal.x = a[0] + (i-maxPoint/2.0);
			normal.y = b[0] + b[1]*(0.5) + b[2]*0.25 + b[3]*0.125;
		}
		else
		{										
			normal.x = a[0] + a[1]*(0.5) + a[2]*0.25 + a[3]*0.125 + (i-maxPoint/2.0)*(b[1]+2*b[2]*0.5+3*b[3]*0.25)/\
								 sqrt((b[1]+2*b[2]*0.5+3*b[3]*0.25)*(b[1]+2*b[2]*0.5+3*b[3]*0.25)+\
								 (a[1]+2*a[2]*0.5+3*a[3]*0.25)*(a[1]+2*a[2]*0.5+3*a[3]*0.25));
											
			normal.y = b[0] + b[1]*(0.5) + b[2]*0.25 + b[3]*0.125 - (i-maxPoint/2.0)*(a[1]+2*a[2]*0.5+3*a[3]*0.25)/\
								 sqrt((b[1]+2*b[2]*0.5+3*b[3]*0.25)*(b[1]+2*b[2]*0.5+3*b[3]*0.25)+\
								 (a[1]+2*a[2]*0.5+3*a[3]*0.25)*(a[1]+2*a[2]*0.5+3*a[3]*0.25));
		}
		norm.Point[i] = normal;
	}
	norm.nPoints = maxPoint;
}

void fit()
{
	a[0] =  1.0/6.0*( fitPoint[0]->x + 4*fitPoint[1]->x +   fitPoint[2]->x );
	a[1] = -1.0/2.0*( fitPoint[0]->x -  								  	fitPoint[2]->x );
	a[2] =  1.0/2.0*( fitPoint[0]->x - 2*fitPoint[1]->x +   fitPoint[2]->x );
	a[3] = -1.0/6.0*( fitPoint[0]->x - 3*fitPoint[1]->x + 3*fitPoint[2]->x - fitPoint[3]->x );
                                                                                              
	b[0] =  1.0/6.0*( fitPoint[0]->y + 4*fitPoint[1]->y +   fitPoint[2]->y );                  
	b[1] = -1.0/2.0*( fitPoint[0]->y -  						  	  	fitPoint[2]->y );                  
	b[2] =  1.0/2.0*( fitPoint[0]->y - 2*fitPoint[1]->y +   fitPoint[2]->y );                  
	b[3] = -1.0/6.0*( fitPoint[0]->y - 3*fitPoint[1]->y + 3*fitPoint[2]->y - fitPoint[3]->y );
}


double CalcSegmentCurvature()//计算拟合线段的曲率
	{
		double SCurvature=0;
		SCurvature = (2*a[1]*b[2]-2*a[2]*b[1]+(6*a[1]*b[3]-6*a[3]*b[1])*0.5+(6*a[2]*b[3]-6*a[3]*b[2])*0.25)/\
											pow(sqrt(pow(a[1]+2*a[2]*0.5+3*a[3]*0.25,2)+pow(b[1]+2*b[2]*0.5+3*b[3]*0.25,2)),3);
		return SCurvature;
	}

double CalcSegmentArcLength()//计算拟合线段的弧长
	{
		double ArcLen=0;
		for(int i=0;i<=2;i++)
		{
			ArcLen += sqrt(pow((a[1]+2*a[2]*(1.0/2*i)+3*a[3]*pow((1.0/2*i),2)),2)+\
			pow((b[1]+2*b[2]*(1.0/2*i)+3*b[3]*pow((1.0/2*i),2)),2))*(1.0/3);
		}
		return ArcLen;
	}



void print()
{
	cout<<"a="<<a[0]<<" "<<a[1]<<" "<<a[2]<<" "<<a[3]<<endl;
	cout<<"b="<<b[0]<<" "<<b[1]<<" "<<b[2]<<" "<<b[3]<<endl;	
}

public:
CPoint *fitPoint[4];
double a[4];//拟合曲线多项式X轴系数
double b[4];//拟合曲线多项式Y轴系数
CProfile norm;
};