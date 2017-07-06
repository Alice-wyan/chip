#pragma once

#include "lib.h"


class CBmpDiff
{
public:
	CBmpDiff(){

	}

	~CBmpDiff(){
		out.clear();
	}	
	
int setup(char *src,char*mod){
IplImage* psrc = cvLoadImage( src, 1 );	
if(psrc==NULL) return 0;
IplImage* pmod = cvLoadImage( mod, 1 );	
if(pmod==NULL) return 0;
CFrame fs;
CFrame fm;
fs.set(psrc);
fm.set(pmod);
CvPoint p;

out.clear();

for(int i=0;i<fs.width;i++){
	for(int j=0;j<fs.height;j++){
		if((fs.b(i,j)!=fm.b(i,j))||(fs.g(i,j)!=fm.g(i,j))||(fs.r(i,j)!=fm.r(i,j))){
			p.x=i;
			p.y=j;
			out.push_back(p);
		}
}}

cvReleaseImage(&psrc);	
cvReleaseImage(&pmod);	
return 1;
}	


void check(vector <CEdgeHitChain> &v){
	ind.clear();
	int bk=0;
	for(unsigned int i=0;i<v.size();i++){
		bk=0;
		for(CEdgeHit *a=v[i].first;a!=NULL;a=a->next){
			for(unsigned int j=0;j<out.size();j++){
				if(out[j].x==a->x&&out[j].y==a->y){
					ind.push_back(i);
					bk=1;
					break;
				}
			
			}
			if(bk==1) break;
		}
	}

}


void checkSword(vector <CEdgeHitChain> &v){
	vps.clear();
vector <CEdgeHit *> ps;
	for(unsigned int i=0;i<ind.size();i++){
		ps.clear();
		for(CEdgeHit *a=v[ind[i]].first;a!=NULL;a=a->next){
			for(unsigned int j=0;j<out.size();j++){
				if(out[j].x==a->x&&out[j].y==a->y){
					ps.push_back(a);
				}
			}
		}

		vps.push_back(ps);
	}

}


void print(CFrame fr){

	for(unsigned int i=0;i<out.size();i++){
		*(fr.pb(out[i].x,out[i].y))=0;
		*(fr.pg(out[i].x,out[i].y))=0;
		*(fr.pr(out[i].x,out[i].y))=0;	
	}

}

vector <CvPoint> out;	
vector <int> ind;	//index of CEdgeHitChain v which has diff point

vector < vector <CEdgeHit *> > vps;

};

