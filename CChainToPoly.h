#pragma once

#include "lib.h"


class CChainToPoly
{
public:
	CChainToPoly(){
		d=5;
	}
	


void transform(CEdgeHitChain ehc,double threshold=0.7,int graph=0)//自定义的transform函数
{
	int n=ehc.num;
		CEdgeHit**eh=new CEdgeHit*[n];
		double *cos=new double [n];
		double *x=new double [n];
	int ind=0;
	for(CEdgeHit*a=ehc.first;a!=NULL;a=a->next){
        eh[ind]=a;
		ind++;
	}
	double vx1;
	double vy1;
	double vx2;
	double vy2;

	for(int i=0;i<n;i++){
	vx1=eh[(i+d)%n]->x-eh[i]->x;
	vy1=eh[(i+d)%n]->y-eh[i]->y;
	vx2=eh[i]->x-eh[(i-d+n)%n]->x;
	vy2=eh[i]->y-eh[(i-d+n)%n]->y;

	    x[i]=i;
		cos[i]=fabs(dot(vx1,vy1,vx2,vy2));
	
	}

	if(graph==1){
	TGraph *g=new TGraph(n,x,cos);
	g->Write();
	delete g;
	}

vector<int> v=	circleMin(cos,n,threshold);
pc.clear();
for(int i=0;i<v.size();i++){

	pc.push_back(eh[v[i]]);

}


	delete [] eh;
	delete [] cos;
	delete [] x;
};

void print(CFrame fr){
	for(int i=0;i<pc.size();i++){
	cvCircle(fr.f,cvPoint(pc[i]->x,pc[i]->y),5,CV_RGB(0,0,255),2,8,0);
	}
};

void print(int h, CFrame fr){
	for(int i=h;i<pc.size();i++){
	cvCircle(fr.f,cvPoint(pc[i]->x,pc[i]->y),5,CV_RGB(0,0,255),2,8,0);
	break;
	}
};


vector<int> circleMin(double *a,int n,double threshold,int r=3){

vector<int> v;
int lm=1;
for(int i=0;i<n;i++){
	lm=1;
	if(a[i]>threshold)continue;

	for(int j=-r;j<=r;j++){
		if(r==0) continue;
		if(a[(i+j+n)%n]<a[i]) {
		lm=0;
		break;
		}

		if(a[(i+j+n)%n]==a[i]&&j<0) {
		lm=0;
		break;
		}



	}
	    if(lm==1){
		    v.push_back(i);
	    }


}

return v;

}

	double dot(double vx1,double vy1,double vx2,double vy2){
	double w=(vx1*vx1+vy1*vy1)*(vx2*vx2+vy2*vy2);
	return (vx1*vx2+vy1*vy2)/sqrt(w);
	}

private:

public:
	int d;
vector <CEdgeHit*> pc;
};


