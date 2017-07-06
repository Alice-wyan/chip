#pragma once

#include "lib.h"



class CSeparateChain
{
public:
CSeparateChain(){

}

void setup(CEdgeHitChain &pref,CEdgeHitChain &pss,CCurveAngle &ssca)
{

	mp=ssca.mp;
	order=ssca.order;
	copyChain(pss);
	fl=pref.num;
}



void copyChain(CEdgeHitChain &pss){

	fc.clear();

	int num=0;
	for(CEdgeHit*a=pss.first;a!=NULL;a=a->next){
		fh[num]=(*a);
		fc.fill(&(fh[num]));
    num++;
	}

}



void sepChain()
{
	vOut.clear();

		CEdgeHit *start = fc.hitIndex(mp);
		if(start==NULL) return ;
		CEdgeHit *ts;
			vt.clear();

	if(fc.checkLoop()==1)
	{
		if(order==0){

			for(int i=0;i<fc.num-fl;i++){
				ts=start;
				if(start==fc.first){
					start=fc.last;
				}else{
					start=start->prev;
				}
				vt.fill(ts);
			
			}
			vOut.push_back(vt);
			
		}else{

			for(int i=0;i<fc.num-fl;i++){
				ts=start;
				if(start==fc.last){
					start=fc.first;
				}else{
					start=start->next;
				}
				vt.fill(ts);
			
			}
			vOut.push_back(vt);		
		}


	}else{
		if(order==0){

			for(int i=0;i<fc.num-fl;i++){
				ts=start;
				if(start==fc.first){
					vt.fill(start);
					vOut.push_back(vt);
					vt.clear();
					start=fc.last;
					continue;

				}else{
					start=start->prev;
				}
				vt.fill(ts);
			
			}
			vOut.push_back(vt);
			
		}else{
			for(int i=0;i<fc.num-fl;i++){
				ts=start;
				if(start==fc.last){
					vt.fill(start);
					vOut.push_back(vt);
					vt.clear();
					start=fc.first;
					continue;

				}else{
					start=start->next;
				}
				vt.fill(ts);
			
			}
			vOut.push_back(vt);		
		
		
		}	
	
	
	
	
	}
}

public:


int mp;
int order;
CEdgeHit fh[1000];


vector <CEdgeHitChain> vOut;
CEdgeHitChain fc;
CEdgeHitChain vt;

int fl;

};

