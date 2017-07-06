#pragma once

#include "lib.h"

class CLinePara
{
public:
	CLinePara(){

	}

void setup(CEdgeHitChain *a){
	if(a->line.size()==0) return;
	int i;
	for(i=0;i<2;i++){
		l[i].clear();
		sita[i].clear();
		dis[i].clear();
	}
	double dx;
	double dy;



	for(i=0;i<int(a->line.size());i++){
		dx=a->line[i].last->x-a->line[i].first->x;
		dy=a->line[i].last->y-a->line[i].first->y;
		l[0].push_back(sqrt(dx*dx+dy*dy));
		sita[0].push_back(atan2(dy,dx));
		dx=a->line[i].last->x-a->line[(i+1)%a->line.size()].first->x;
		dy=a->line[i].last->y-a->line[(i+1)%a->line.size()].first->y;
		dis[0].push_back(sqrt(dx*dx+dy*dy));
	}

	for(i=int(a->line.size())-1;i>=0;i--){

		l[1].push_back(l[0][i]);
		sita[1].push_back(sita[0][(i-1+a->line.size())%a->line.size()]-sita[0][i]);
		dis[1].push_back(dis[0][(i-1+a->line.size())%a->line.size()]);
	}


	double s=sita[0][0];

	for(i=0;i<int(a->line.size());i++){

		if(i<int(a->line.size())-1){
		sita[0][i]=sita[0][i]-sita[0][(i+1)%a->line.size()];
		}
		else{
		sita[0][i]=sita[0][i]-s;
		}

	}
	for(i=0;i<int(a->line.size());i++){
	if(sita[0][i]>PI)sita[0][i]-=2*PI;
	if(sita[0][i]<=-PI)sita[0][i]+=2*PI;
	if(sita[1][i]>PI)sita[1][i]-=2*PI;
	if(sita[1][i]<=-PI)sita[1][i]+=2*PI;
	}
	
}	
	
vector<int> match(CLinePara *a){
	double dl=5;
	double dsita=0.3;
	double ddis=5;
    int extend=5; 

	vector<int> ret;	
	int s=int(l[0].size());
	int sa=int(a->l[0].size());
	if(s==0) return ret;
	if(sa==0) return ret;

	int bk;
	int i0,j0,order,bkMax;

	bkMax=0;

	for(int i=0;i<s;i++){
	for(int j=0;j<sa;j++){
	if(fabs(l[0][i]-a->l[0][j])>dl) continue;
	if(fabs(sita[0][i]-a->sita[0][j])>dsita) continue;
	if(fabs(dis[0][i]-a->dis[0][j])>ddis) continue;
	bk=0;
	for(int k=1;k<1+extend;k++){
	if(fabs(l[0][(i+k)%s]-a->l[0][(j+k)%sa])>dl) break;
	if(fabs(sita[0][(i+k)%s]-a->sita[0][(j+k)%sa])>dsita) break;
	if(fabs(dis[0][(i+k)%s]-a->dis[0][(j+k)%sa])>ddis) break;	
	bk=k;
	}

	if(bk>bkMax){
	bkMax=bk;
	i0=i;
	j0=j;
	order=0;
	}

	}
	}




	for(int i=0;i<s;i++){
	for(int j=0;j<sa;j++){
	if(fabs(l[1][i]-a->l[0][j])>dl) continue;
	if(fabs(sita[1][i]-a->sita[0][j])>dsita) continue;
	if(fabs(dis[1][i]-a->dis[0][j])>ddis) continue;
	bk=0;
	for(int k=1;k<1+extend;k++){
	if(fabs(l[1][(i+k)%s]-a->l[0][(j+k)%sa])>dl) break;
	if(fabs(sita[1][(i+k)%s]-a->sita[0][(j+k)%sa])>dsita) break;
	if(fabs(dis[1][(i+k)%s]-a->dis[0][(j+k)%sa])>ddis) break;	
	bk=k;
	}

	if(bk>bkMax){
	bkMax=bk;
	i0=i;
	j0=j;
	order=1;
	}

	}
	}

	if(bkMax>0){
		ret.push_back(i0);
		ret.push_back(j0);
		ret.push_back(order);
		ret.push_back(bkMax);
	}

	return ret;
}

void print(){
	unsigned int i;
	for(i=0;i<l[0].size();i++){
    cout<<i<<" >>  "<<l[0][i]<<" "<<sita[0][i]<<" "<<dis[0][i]<<endl;	
	}
	for(i=0;i<l[0].size();i++){
    cout<<i<<" <<  "<<l[1][i]<<" "<<sita[1][i]<<" "<<dis[1][i]<<endl;	
	}
}

vector <double> l[2];
vector <double> sita[2];
vector <double> dis[2];

};

