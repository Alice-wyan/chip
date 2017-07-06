#pragma once

#include "lib.h"



class CEdgeHitChain
{
public:
CEdgeHitChain(){
clear();
}
void clear(){
	first=NULL;
	last=NULL;
	num=0;
	line.clear();

}
CEdgeHit *first;
CEdgeHit *last;

CEdgeHit *one;
CEdgeHit *two;
double L;
double checkStep;


int num;

vector <CLine> line;



double loopDis(){
	return double((first->x-last->x)*(first->x-last->x)+(first->y-last->y)*(first->y-last->y));
}


int checkLoop(){
	if(loopDis()>16) return 0;
	else return 1;
}

//int checkLoop(CPointToPeak * pp){
//	if(pp==NULL) return 0;
//	if(pp->tryConnect(first->x,first->y,last->x,last->y)<0) return 0;
//	else return 1;
//}

int mergeLine(CEdgeHit**eh,double dw=2,int pre=3,int calLength=5,int goodLength=15){
	if(line.size()<2) return 0;
	if(checkLoop()==0) return 0;
	if(line[0].first!=first) return 0;
	if(line[line.size()-1].last!=last) return 0;

	int sum=0;

	for(int i=num-1;i>=0;i--){
			if(i-pre>=0){
				//cout<<i<<"  le.check(eh[i+pre])="<<le.check(eh[i+pre])<<endl;
				if(line[0].check(eh[i-pre])>dw){
					line[0].first=eh[(i+1)%num];
					break;
					
				}
				line[0].add(eh[i]);
				sum++;
			}	
	
	}
	//cout<<"in Merge line "<<sum<<endl;
	if(line[line.size()-1].nHits-sum<goodLength){line.pop_back();
	//cout<<"line pop"<<endl;
	}
	return 1;
	//cout<<line[0].k()<<" "<<line[0].b()<<" vs  "<<line[line.size()-1].k()<<" "<<line[line.size()-1].b()<<endl;
}

int checkLine(double dw=2,int pre=3,int calLength=5,int goodLength=15){
    CLine le;

	
	
	

	CEdgeHit**eh=new CEdgeHit*[num];
	int ind=0;
	for(CEdgeHit*a=first;a!=NULL;a=a->next){
        eh[ind]=a;
		ind++;
	}
	le.clear();
	//cout<<ind<<"   ind vs num  "<<num<<endl;
	for(int i=0;i<num;i++){
		le.add(eh[i]);
		if(le.nHits==1)le.first=eh[i];

		if(le.nHits>calLength){
			if(i+pre<num){
				//cout<<i<<"  le.check(eh[i+pre])="<<le.check(eh[i+pre])<<endl;
				if(le.check(eh[i+pre])>dw){
					if(le.nHits>goodLength){
					le.last=eh[i];
					line.push_back(le);
					}

					le.clear();
				}
			}else{
				//cout<<i<<"  le.check(eh[num-1])="<<le.check(eh[num-1])<<endl;
				if(le.check(eh[num-1])>dw){
					if(le.nHits>goodLength){
					le.last=eh[i];
					line.push_back(le);
					}
					le.clear();
				}else if(i==num-1){
					if(le.nHits>goodLength){
					le.last=eh[i];
					line.push_back(le);
					}

					le.clear();			
			
			}			
			
			}



		}
	}


mergeLine(eh,dw,pre,calLength,goodLength);

	delete []eh;
	return 0;
}


void fill(CEdgeHit*a){
	if(first==NULL){
	first=a;
	a->prev=NULL;
	}else{
	last->next=a;
	a->prev=last;
	}//if(first==NULL)
	last=a;
    a->next=NULL;
	num++;
	a->used=1;
	
}


void fillBack(CEdgeHit*a){
	if(last==NULL){
	last=a;
	a->next=NULL;	
	}else{
	first->prev=a;		
	a->next=first;
	}
	first=a;
	a->prev=NULL;
	num++;	
	a->used=1;
}

double dis(CEdgeHit*a,CEdgeHit*b){
	double dx=a->x-b->x;
	double dy=a->y-b->y;
return sqrt(dx*dx+dy*dy);
}
int resetTwoPoint(double len){
	L=len;
	one=first;
	two=NULL;
	for(CEdgeHit*a=first;a!=NULL;a=a->next){
		if(dis(one,a)>L){
		two=a;
		break;
		}
	}
	if(two==NULL) return 0;
	else return 1;
}

double twoPointAngle(){
	assert(one!=NULL&&two!=NULL);
	return atan2(double(two->y-one->y),double(two->x-one->x));
}



int checkTwoPoint(){
	CEdgeHit*pone=one;
	one=one->next;
	if(one==NULL) return 0;

	int disi=abs(pone->y-one->y)+abs(pone->x-one->x);
	if(disi==1) checkStep= 1;
	else if(disi==2) checkStep= 1.414213562373;
	else{
	disi=(pone->y-one->y)*(pone->y-one->y);
	disi+=(pone->x-one->x)*(pone->x-one->x);
	checkStep= sqrt(double(disi));
	}


	CEdgeHit*a=two;
	two=NULL;
	for(;a!=NULL;a=a->next){
		if(dis(one,a)>L){
		two=a;
		break;
		}
	}
	if(two!=NULL) return 1;

	if(this->checkLoop()==0) return 0;

	for(a=first;a!=NULL;a=a->next){

		if(dis(one,a)>L){
		two=a;
		break;			
	
		}
	}
	if(two==NULL) return 0;
	else return 1;

}

CEdgeHit* hitIndex(int ind){
	CEdgeHit* a=first;
	for(int i=0;i<num;i++){
	if(i==ind) return a;
	a=a->next;
	}
	return NULL;
}


CEdgeHit* maxDisHit(CEdgeHit* s,double & d,int &ind){
	CEdgeHit* a=first;
	CEdgeHit* max=first;
	double l;
	d=0;
	for(int i=0;i<num;i++){
		l=dis(s,a);
		if(l>d){
		d=l;
		max=a;
		ind=i;
		}
	a=a->next;
	}

	return max;
}


void print(CFrame fr,int drawCircle=0,int drawLine=0,int r=255,int g=0,int b=0){

	//int r=rand()%255;
	//int g=rand()%255;
	//int b=rand()%255;
	//r=255;
	//g=0;
	//b=0;

	for(CEdgeHit*a=first;a!=NULL;a=a->next){
		
	*(fr.pb(a->x,a->y))=b;
	*(fr.pg(a->x,a->y))=g;
	*(fr.pr(a->x,a->y))=r;
	}

	
if(drawCircle==1){
if(checkLoop()==0){
CEdgeHit*a=first;
cvCircle(fr.f,cvPoint(a->x,a->y),3,CV_RGB(0,0,255),2,8,0);

a=last;
cvCircle(fr.f,cvPoint(a->x,a->y),3,CV_RGB(0,255,0),2,8,0);
}
}

if(drawLine==1){

for(unsigned int i=0;i<line.size();i++){
	CvPoint pt1;
	CvPoint pt2;
	pt1.x=line[i].first->x;
	pt1.y=line[i].first->y;
	pt2.x=line[i].last->x;
	pt2.y=line[i].last->y;


cvCircle(fr.f,pt1,3,CV_RGB(0,0,255),2,8,0);
cvCircle(fr.f,pt2,3,CV_RGB(0,255,0),2,8,0);
cvLine(fr.f,pt1,pt2,CV_RGB(255,255,0));
}

}


}//void print(CFrame fr)


};

