#pragma once

#include "lib.h"


class CPolyTag
{
public:
	CPolyTag(){
	
	}
	


void transform(vector <CEdgeHit*> pc)
{
int n=pc.size();
double *ta=new double[n-1];
int *ti=new int[n-1];
range.clear();
ang.clear();
dis.clear();
for(int i=0;i<n;i++){// for any point in n
//if(i==1) break;
	double max=0;
	double tm=0;
for(int j=0;j<n;j++){// vector between other point and base point
	if(i==j) continue;
	if(j<i){
	tm=mod(pc[i],pc[j]);		
	}else{
	tm=mod(pc[i],pc[j]);				
	}

	if(tm>max)max=tm;

}	
dis.push_back(max);

for(int j=0;j<n;j++){// vector between other point and base point
	if(i==j) continue;
	if(j<i){
	ta[j]=angle(pc[i],pc[j]);		
	ti[j]=j;	
	}else{
	ta[j-1]=angle(pc[i],pc[j]);		
	ti[j-1]=j;			
	}
}	

BubbleSort(ta,ti,n-1);
int m=maxIncludedAngleIndex(ta,n-1);


double mrange=ta[(m+1)%(n-1)]-ta[m];
if(m==n-2) mrange+=2*PI;
//cout<<"m="<<m<<"   mrange="<<mrange<<endl;
range.push_back(mrange);

double mang=(ta[(m+1)%(n-1)]+ta[m])/2;
if(m==n-2){
if(mang>=0) mang-=PI; else mang+=PI;
} 
//cout<<"m="<<m<<"   mang="<<mang<<endl;
ang.push_back(mang);
}// for(int i   for any point in n

delete []ta;
delete []ti;

};

int maxIncludedAngleIndex(double * R, int n){
	double max=0;
	int maxi=-1;
	
	if(n==1) return 0;
	
	for(int i=0;i<n-1;i++){
		if(R[i+1]-R[i]>max){
			max=R[i+1]-R[i];
			maxi=i;
		}
	}
	
	if(R[0]-R[n-1]+2*PI>max){
		max=R[0]-R[n-1]+2*PI;
		maxi=n-1;
	}	
	return maxi;
}

void BubbleSort(double * R,int * ind, int n,int print=0)
{
	double temp;
	int    temi;
    for (int i = 0; i < n - 1; i++)
    {
        bool noswap = true;
        for (int j = 0; j < n - 1-i; j++)
        {
            if (R[j] > R[j + 1])
            {
                temp= R[j];
                R[j] = R[j + 1];
                R[j + 1] = temp;
				
                temi = ind[j];
                ind[j] = ind[j + 1];
                ind[j + 1] = temi;				
				
                noswap = false;
            }
        }
        if (noswap)
        {
            break;
        }
    }
	
	if(print==1){
		cout<<"=====print in BubbleSort  start========="<<endl;
	for(int i=0;i<n;i++){
		cout<<R[i]<<"  "<<ind[i]<<endl;
	}
		cout<<"=====print in BubbleSort  end========="<<endl;
	
	}
}
void rangeFilter(double mid, double error){
	rangeFil.clear();
	for(int i=0;i<range.size();i++){
		if(fabs(range[i]-mid)>error) continue;
		rangeFil.push_back(i);
	}
}

int angleFilter(double mid, double error){

	for(int i=0;i<rangeFil.size();i++){
		
		if(fabs(ang[rangeFil[i]]-mid)>error) continue;
		
		return rangeFil[i];

	}
	return -1;
}

void testSort(){
	double * r=new double [5];
	int * ri=new int [5];

	r[0]=3.5;
	r[1]=4.5;
	r[2]=2.5;
	r[3]=5.5;
	r[4]=1.5;

	ri[0]=1;
	ri[1]=2;
	ri[2]=3;
	ri[3]=4;
	ri[4]=5;	
	
	BubbleSort(r,ri,5);
	
	for(int i=0;i<5;i++){
		cout<<r[i]<<"  "<<ri[i]<<endl;
	}
	
	
}

 double angle(CEdgeHit * base, CEdgeHit * target){ //angle is from base to target
	return  atan2((double)target->x-base->x,(double)target->y-base->y);
}

double mod(CEdgeHit * base, CEdgeHit * target){ 
double dx=target->x-base->x;
double dy=target->y-base->y;
	return  sqrt(dx*dx+dy*dy);
}

void print(CFrame fr){

};




private:

public:

vector <double> ang;
vector <double> range;
vector <double> dis;
vector <int> rangeFil;
};


