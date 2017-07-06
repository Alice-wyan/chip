#pragma once

#include "lib.h"



class CCurveAngle
{
public:
CCurveAngle(){
len=10;
debug=1;
}

void getCa(CEdgeHitChain & v){
	pv=&v;

			vs.clear();
			va.clear();
			vq.clear();
			ts=0;
	if(v.resetTwoPoint(len)){
		for(int j=0;j<v.num;j++){
			va.push_back(v.twoPointAngle());
			if(v.checkTwoPoint()==0)break;	
			vs.push_back(ts);
			ts+=v.checkStep;
		}
		loop=v.checkLoop();

		if(loop==1){
		vs.push_back(ts);
		ts+=sqrt(v.loopDis());
		}else{
			va.pop_back();
		}


		int q=0;
		for(unsigned int k=0;k<ts;k++){
			if(vs[q]<k) q++;
			if(q>=va.size())q=va.size()-1;
			vq.push_back(q);
		}

}


}


double minDiff(CCurveAngle & b){
	pc=&b;
	// b is the short one
	mp=-1;
	order=-1;
	double min=10;
	double diff=0;
	double sum=0;
	double a=0;
	int sl;

	if(loop)sl=vq.size();
	else sl=vq.size()-b.vq.size();

	for(int i=0;i<sl;i++){
		diff=trim(v(i)-b.v(0));
		sum=0;
	for(unsigned int j=0;j<b.vq.size();j++){
		a=trim(trim(v(i+j)-b.v(j))-diff);
	    sum+=a*a;
	}
	sum/=b.vq.size();
	sum=sqrt(sum);

	if(sum<min){
		min=sum;
		mp=i;
		order=0;
	}
	}


	for(int i=0;i<sl;i++){
		diff=trim(v(i)-b.v(b.vq.size()-1));
		sum=0;
	for(unsigned int j=0;j<b.vq.size();j++){
		a=trim(trim(v(i+j)-b.v(b.vq.size()-1-j))-diff);
	    sum+=a*a;
	}
	sum/=b.vq.size();
	sum=sqrt(sum);

	if(sum<min){
		min=sum;
		mp=i;
		order=1;
	}
	}
	return min;
}

double v(int i){
	return va[vq[i%(vq.size())]];
}



double trim(double a){
     if(a>PI) return a-2*PI;
     if(a<=-PI) return a+2*PI;
	 return a;
}


void posDiff(CCurveAngle & b,CSixSword & v){
			 sum=0;
			 max=0;
			 double tsum;
	if(order==0){
			sfirst=pv->hitIndex(vq[mp]);
			smid=pv->hitIndex((vq[mp]+b.find)%pv->num);
			v.ccs(b.pv->first,b.fmid,sfirst,smid);

			dx=sfirst->x-b.pv->first->x;
			dy=sfirst->y-b.pv->first->y;
			s=v.s;
			c=v.c;


			CEdgeHit*bh=sfirst;
	for(CEdgeHit*a=b.pv->first;a!=NULL;a=a->next){
        v.expect(b.pv->first,a,sfirst);
		tsum=(bh->x-v.x)*(bh->x-v.x)+(bh->y-v.y)*(bh->y-v.y);
		sum+=tsum;
		if(tsum>max){
		max=tsum;
		}

		if(bh==pv->last){
			bh=pv->first;
		}else {
		    bh=bh->next;
		}
		//		if(v.x<0) continue;
		//		if(v.x>=v.fr.width) continue;
		//		if(v.y<0) continue;
		//		if(v.y>=v.fr.height) continue;
  //      if(debug>=4){
		//*v.fr.pb(v.x,v.y)=0;
		//*v.fr.pg(v.x,v.y)=0;
		//*v.fr.pr(v.x,v.y)=0;
		//}

	}
	if(debug>=2)cout<<"===="<<sqrt(sum/double(b.pv->num))<<endl;
	}


	if(order==1){

			sfirst=pv->hitIndex(vq[mp]);
			smid=pv->hitIndex((vq[mp]+(b.pv->num-b.lind-1))%pv->num);
			v.ccs(b.pv->last,b.lmid,sfirst,smid);

			dx=sfirst->x-b.pv->last->x;
			dy=sfirst->y-b.pv->last->y;
			s=v.s;
			c=v.c;
	
			CEdgeHit*bh=sfirst;
			for(CEdgeHit*a=b.pv->last;a!=NULL;a=a->prev){
        v.expect(b.pv->last,a,sfirst);

		tsum=(bh->x-v.x)*(bh->x-v.x)+(bh->y-v.y)*(bh->y-v.y);
		sum+=tsum;
		if(tsum>max){
		max=tsum;
		}

		if(bh==pv->last){
			bh=pv->first;
		}else {
		    bh=bh->next;
		}
		//		if(v.x<0) continue;
		//		if(v.x>=v.fr.width) continue;
		//		if(v.y<0) continue;
		//		if(v.y>=v.fr.height) continue;
		//		        if(debug>=4){
		//*v.fr.pb(v.x,v.y)=0;
		//*v.fr.pg(v.x,v.y)=0;
		//*v.fr.pr(v.x,v.y)=0;
		//				}
	}
	if(debug>=2)cout<<"****"<<sqrt(sum/double(b.pv->num))<<endl;	
	
	}


}

void setTwoHit(){


pv->maxDisHit(pv->first,d,find);
fmid=pv->hitIndex(find);
pv->maxDisHit(pv->last,d,lind);
lmid=pv->hitIndex(lind);

}



public:



double len;
double ts;
int loop;
int mp;
int order;
vector <double> vs;
vector <double> va;
vector <int> vq;
CEdgeHitChain *pv;

CCurveAngle * pc;

CSixSword *ps;

double d;
int find;
int lind;
CEdgeHit* fmid;
CEdgeHit* lmid;

CEdgeHit* sfirst;
CEdgeHit* smid;

double sum;
double max;
double dx;
double dy;
double s;
double c;

int debug;

};

