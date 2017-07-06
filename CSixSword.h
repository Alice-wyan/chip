#pragma once

#include "lib.h"


class CSixSword
{
public:
CSixSword(){
	f=NULL;
}
~CSixSword(){
	if(f!=NULL) cvReleaseImage(&f);
}

void setup(char *s){
	  f=cvLoadImage( s, 1 );	
	  fr.set(f);
      pe.setup(f,200000,1);	
	  newEdge = cvCreateImage(cvGetSize(f),f->depth,4);
	  frNewEdge.set(newEdge);
	  pp.frIn.set(f);
}

void setup(IplImage* pF){
	  fr.set(pF);
      pe.setup(pF,200000,1);//frDst.set(dst);把dst  set 给	frDst
	  newEdge = cvCreateImage(cvGetSize(pF),pF->depth,4);
	  frNewEdge.set(newEdge);////for read
	  pp.frIn.set(pF);
	  //cvReleaseImage(&newEdge);
}

void transform(TH1D *b=NULL){


		//int nCheck=2;
		double mergeCut=4;
		int minHits=20;

        double dw=2;
        int pre=3;
        int calLength=4;
        int goodLength=10;
//////////////////////////////////////////


		//cvCvtColor(pF,gray,CV_BGR2GRAY);
		vehc.clear();
		unsigned int ind;
		int filled;
		int x;
		int y;
		//int ii;
		//int jj;
		int dx0;
		int dy0;
		int dx1;
		int dy1;
		t.debug=0;
		t.start();
		//cvZero(pe.dst);
		pe.transformHit(newEdge);//outHit里存储着CEdgeHit eh

		t.end();
		t.print("transformHit");
		//cout<<pe.outHit.size()<<endl;
		//for(int i=0;i<pe.outHit.size();i++){
		//	x=pe.outHit[i].x;
		//	y=pe.outHit[i].y;
		//	if(i!=((int*)pe.posToIndex.pb(x,y))[0])cout<<i<<"=="<<((int*)pe.posToIndex.pb(x,y))[0]<<endl;
		//}
		t.start();
		for(unsigned int i=0;i<pe.outHit.size();i++){
			if(pe.outHit[i].used==1) continue;

			ehc.clear();
			ehc.fill(&(pe.outHit[i]));//ehc=outHit[i]  链表

			x=pe.outHit[i].x;
			y=pe.outHit[i].y;

			do{
				filled=0;

				//for(int k=1;k<nCheck;k++){
				//for(int m=0;m<4;m++){
				//for(int n=-k;n<k;n++){
				//	if(m==0) {ii=-k;jj=n;}
				//	if(m==1) {ii=k;jj=-n;}
				//	if(m==2) {ii=-n;jj=-k;}
				//	if(m==3) {ii=n;jj=k;}

				for(int ii=-1;ii<=1;ii++){//-1 0 1
				for(int jj=-1;jj<=1;jj++){

					if(ii==0&&jj==0) continue;
					if(x+ii<0)continue;
					if(x+ii>=pe.frDst.width)continue;//排除一圈边框的轮廓
					if(y+jj<0)continue;
					if(y+jj>=pe.frDst.height)continue;

					if(pe.frDst.b(x+ii,y+jj)==0) continue;

					ind=((int*)pe.posToIndex.pb(x+ii,y+jj))[0];//b通道
					//if(ind<0||ind>=pe.outHit.size())cout<<pe.outHit.size()<<" "<<x<<"  "<<y<<" "<<ii<<" "<<jj<<" ind="<<ind<<endl;
					if(pe.outHit[ind].used==1) continue;

					ehc.fill(&(pe.outHit[ind]));//ehc=outHit[ind]  链表

					filled=1;
					x+=ii;
					y+=jj;
					if(filled==1) break;
				}
				if(filled==1) break;
				}
				//if(filled==1) break;
				//}
				if(filled==0) break;

			}while(filled==1);

	
			x=pe.outHit[i].x;
			y=pe.outHit[i].y;
			do{
				filled=0;


				//for(int k=1;k<nCheck;k++){
				//for(int m=0;m<4;m++){
				//for(int n=-k;n<k;n++){
				//	if(m==0) {ii=-k;jj=n;}
				//	if(m==1) {ii=k;jj=-n;}
				//	if(m==2) {ii=-n;jj=-k;}
				//	if(m==3) {ii=n;jj=k;}

				for(int ii=-1;ii<=1;ii++){
				for(int jj=-1;jj<=1;jj++){

					if(ii==0&&jj==0) continue;
					if(x+ii<0)continue;
					if(x+ii>=pe.frDst.width)continue;
					if(y+jj<0)continue;
					if(y+jj>=pe.frDst.height)continue;

					if(pe.frDst.b(x+ii,y+jj)==0) continue;
					ind=((int*)pe.posToIndex.pb(x+ii,y+jj))[0];
					//if(ind<0||ind>=pe.outHit.size())cout<<pe.outHit.size()<<" back "<<x<<"  "<<y<<" "<<ii<<" "<<jj<<" ind="<<ind<<endl;
					if(pe.outHit[ind].used==1) continue;

					ehc.fillBack(&(pe.outHit[ind]));//倒着吗？
					filled=1;
					x+=ii;
					y+=jj;
					if(filled==1) break;
				}
				if(filled==1) break;
				}
				//if(filled==1) break;
				//}
				if(filled==0) break;
			}while(filled==1);

			//cout<<pe.outHit.size()<<" "<<i<<" ehc="<<ehc.num<<endl;
			//if(ehc.num>5)
				//ehc.print(fr);
			vehc.push_back(ehc);
		}
		//cout<<"vehc="<<vehc.size()<<endl;

		t.end();
		t.print("connect");



		int min=10000;
		int minPos=-1;
		int order=-1;
		int pow2;

		t.start();

		for(int i=vehc.size()-1;i>=0;i--){
		 min=10000;
		 minPos=-1;


		for(int j=0;j<i;j++){

			dx0=vehc[i].first->x-vehc[j].first->x;//尾首相减即算出长度
			dy0=vehc[i].first->y-vehc[j].first->y;
			dx1=vehc[i].first->x-vehc[j].last->x;
			dy1=vehc[i].first->y-vehc[j].last->y;

			if((abs(dx0)>mergeCut||abs(dy0)>mergeCut)&&(abs(dx1)>mergeCut||abs(dy1)>mergeCut)) {continue;}//mergeCut=4

			pow2=dx0*dx0+dy0*dy0;
			if(min>pow2) {min=pow2;minPos=j;order=0;}
			pow2=dx1*dx1+dy1*dy1;
			if(min>pow2) {min=pow2;minPos=j;order=1;}


		}

		if(sqrt(double(min))<mergeCut&&vehc[i].loopDis()>min){
			//cout<<i<<" ---  "<<minPos<<"  "<<vehc[i].num<<"  "<<vehc[minPos].num<<" "<<order<<" "<<min<<endl;
			if(order==0){

	CEdgeHit* pn;
	for(CEdgeHit*a=vehc[i].first;a!=NULL;a=pn){
	pn=a->next;
	vehc[minPos].fillBack(a);
	}			
			}else{
				vehc[minPos].last->next=vehc[i].first;
				vehc[i].first->prev=vehc[minPos].last;
				vehc[minPos].last=vehc[i].last;
				vehc[minPos].num+=vehc[i].num;
			}
	vehc[i].clear();
	continue;
		};



		 min=10000;
		 minPos=-1;
		for(int j=0;j<i;j++){
			dx0=vehc[i].last->x-vehc[j].first->x;
			dy0=vehc[i].last->y-vehc[j].first->y;
			dx1=vehc[i].last->x-vehc[j].last->x;
			dy1=vehc[i].last->y-vehc[j].last->y;

			if((abs(dx0)>mergeCut||abs(dy0)>mergeCut)&&(abs(dx1)>mergeCut||abs(dy1)>mergeCut)) {continue;}

			pow2=dx0*dx0+dy0*dy0;
			if(min>pow2) {min=pow2;minPos=j;order=0;}
			pow2=dx1*dx1+dy1*dy1;
			if(min>pow2) {min=pow2;minPos=j;order=1;}
		}

		if(sqrt(double(min))<mergeCut&&vehc[i].loopDis()>min){
			//if(cut==1)cout<<"what****"<<endl;
			//cout<<i<<" ---  "<<minPos<<"  "<<vehc[i].num<<"  "<<vehc[minPos].num<<" "<<order<<" "<<min<<endl;
			if(order==0){

				vehc[minPos].first->prev=vehc[i].last;
				vehc[i].last->next=vehc[minPos].first;
				vehc[minPos].first=vehc[i].first;
				vehc[minPos].num+=vehc[i].num;

	
			}else{
	CEdgeHit* pn;
	for(CEdgeHit*a=vehc[i].last;a!=NULL;a=pn){
	pn=a->prev;
	vehc[minPos].fill(a);
	}	
			}
	vehc[i].clear();	
	continue;
		};


	}		

		t.end();
		t.print("merge");






		v.clear();
		int sum=0;
		for(unsigned int i=0;i<vehc.size();i++){
			if(vehc[i].first==NULL) continue;
			if(vehc[i].num>minHits){
				v.push_back(vehc[i]);
			}
		}
		CLinePara lp;
		vlp.clear();
		for(unsigned int i=0;i<v.size();i++){
			v[i].checkLine(dw,pre,calLength,goodLength);
			lp.setup(&(v[i]));
			vlp.push_back(lp);
		}
		  cvReleaseImage(&newEdge);

}


void transformMask(vector <CvPoint> & mask){


		//int nCheck=2;
		double mergeCut=4;
		int minHits=20;

        double dw=2;
        int pre=3;
        int calLength=4;
        int goodLength=10;
//////////////////////////////////////////


		//cvCvtColor(pF,gray,CV_BGR2GRAY);
		vehc.clear();
		unsigned int ind;
		int filled;
		int x;
		int y;
		//int ii;
		//int jj;
		int dx0;
		int dy0;
		int dx1;
		int dy1;
		t.debug=0;
		t.start();
		//cvZero(pe.dst);
		pe.transformHitMask(newEdge,mask);

		t.end();
		t.print("transformHit");
		//cout<<pe.outHit.size()<<endl;
		//for(int i=0;i<pe.outHit.size();i++){
		//	x=pe.outHit[i].x;
		//	y=pe.outHit[i].y;
		//	if(i!=((int*)pe.posToIndex.pb(x,y))[0])cout<<i<<"=="<<((int*)pe.posToIndex.pb(x,y))[0]<<endl;
		//}
		t.start();
		for(unsigned int i=0;i<pe.outHit.size();i++){
			if(pe.outHit[i].used==1) continue;

			ehc.clear();
			ehc.fill(&(pe.outHit[i]));

			x=pe.outHit[i].x;
			y=pe.outHit[i].y;

			do{
				filled=0;

				//for(int k=1;k<nCheck;k++){
				//for(int m=0;m<4;m++){
				//for(int n=-k;n<k;n++){
				//	if(m==0) {ii=-k;jj=n;}
				//	if(m==1) {ii=k;jj=-n;}
				//	if(m==2) {ii=-n;jj=-k;}
				//	if(m==3) {ii=n;jj=k;}

				for(int ii=-1;ii<=1;ii++){
				for(int jj=-1;jj<=1;jj++){

					if(ii==0&&jj==0) continue;
					if(x+ii<0)continue;
					if(x+ii>=pe.frDst.width)continue;
					if(y+jj<0)continue;
					if(y+jj>=pe.frDst.height)continue;

					if(pe.frDst.b(x+ii,y+jj)==0) continue;

					ind=((int*)pe.posToIndex.pb(x+ii,y+jj))[0];
					//if(ind<0||ind>=pe.outHit.size())cout<<pe.outHit.size()<<" "<<x<<"  "<<y<<" "<<ii<<" "<<jj<<" ind="<<ind<<endl;
					if(pe.outHit[ind].used==1) continue;

					ehc.fill(&(pe.outHit[ind]));

					filled=1;
					x+=ii;
					y+=jj;
					if(filled==1) break;
				}
				if(filled==1) break;
				}
				//if(filled==1) break;
				//}
				if(filled==0) break;

			}while(filled==1);

	
			x=pe.outHit[i].x;
			y=pe.outHit[i].y;
			do{
				filled=0;


				//for(int k=1;k<nCheck;k++){
				//for(int m=0;m<4;m++){
				//for(int n=-k;n<k;n++){
				//	if(m==0) {ii=-k;jj=n;}
				//	if(m==1) {ii=k;jj=-n;}
				//	if(m==2) {ii=-n;jj=-k;}
				//	if(m==3) {ii=n;jj=k;}

				for(int ii=-1;ii<=1;ii++){
				for(int jj=-1;jj<=1;jj++){

					if(ii==0&&jj==0) continue;
					if(x+ii<0)continue;
					if(x+ii>=pe.frDst.width)continue;
					if(y+jj<0)continue;
					if(y+jj>=pe.frDst.height)continue;

					if(pe.frDst.b(x+ii,y+jj)==0) continue;
					ind=((int*)pe.posToIndex.pb(x+ii,y+jj))[0];
					//if(ind<0||ind>=pe.outHit.size())cout<<pe.outHit.size()<<" back "<<x<<"  "<<y<<" "<<ii<<" "<<jj<<" ind="<<ind<<endl;
					if(pe.outHit[ind].used==1) continue;

					ehc.fillBack(&(pe.outHit[ind]));
					filled=1;
					x+=ii;
					y+=jj;
					if(filled==1) break;
				}
				if(filled==1) break;
				}
				//if(filled==1) break;
				//}
				if(filled==0) break;
			}while(filled==1);

			//cout<<pe.outHit.size()<<" "<<i<<" ehc="<<ehc.num<<endl;
			//if(ehc.num>5)
				//ehc.print(fr);
			vehc.push_back(ehc);
		}
		//cout<<"vehc="<<vehc.size()<<endl;

		t.end();
		t.print("connect");



		int min=10000;
		int minPos=-1;
		int order=-1;
		int pow2;

		t.start();

		for(int i=vehc.size()-1;i>=0;i--){
		 min=10000;
		 minPos=-1;


		for(int j=0;j<i;j++){

			dx0=vehc[i].first->x-vehc[j].first->x;
			dy0=vehc[i].first->y-vehc[j].first->y;
			dx1=vehc[i].first->x-vehc[j].last->x;
			dy1=vehc[i].first->y-vehc[j].last->y;

			if((abs(dx0)>mergeCut||abs(dy0)>mergeCut)&&(abs(dx1)>mergeCut||abs(dy1)>mergeCut)) {continue;}

			pow2=dx0*dx0+dy0*dy0;
			if(min>pow2) {min=pow2;minPos=j;order=0;}
			pow2=dx1*dx1+dy1*dy1;
			if(min>pow2) {min=pow2;minPos=j;order=1;}


		}

		if(sqrt(double(min))<mergeCut&&vehc[i].loopDis()>min){
			//cout<<i<<" ---  "<<minPos<<"  "<<vehc[i].num<<"  "<<vehc[minPos].num<<" "<<order<<" "<<min<<endl;
			if(order==0){

	CEdgeHit* pn;
	for(CEdgeHit*a=vehc[i].first;a!=NULL;a=pn){
	pn=a->next;
	vehc[minPos].fillBack(a);
	}			
			}else{
				vehc[minPos].last->next=vehc[i].first;
				vehc[i].first->prev=vehc[minPos].last;
				vehc[minPos].last=vehc[i].last;
				vehc[minPos].num+=vehc[i].num;
			}
	vehc[i].clear();
	continue;
		};



		 min=10000;
		 minPos=-1;
		for(int j=0;j<i;j++){
			dx0=vehc[i].last->x-vehc[j].first->x;
			dy0=vehc[i].last->y-vehc[j].first->y;
			dx1=vehc[i].last->x-vehc[j].last->x;
			dy1=vehc[i].last->y-vehc[j].last->y;

			if((abs(dx0)>mergeCut||abs(dy0)>mergeCut)&&(abs(dx1)>mergeCut||abs(dy1)>mergeCut)) {continue;}

			pow2=dx0*dx0+dy0*dy0;
			if(min>pow2) {min=pow2;minPos=j;order=0;}
			pow2=dx1*dx1+dy1*dy1;
			if(min>pow2) {min=pow2;minPos=j;order=1;}
		}

		if(sqrt(double(min))<mergeCut&&vehc[i].loopDis()>min){
			//if(cut==1)cout<<"what****"<<endl;
			//cout<<i<<" ---  "<<minPos<<"  "<<vehc[i].num<<"  "<<vehc[minPos].num<<" "<<order<<" "<<min<<endl;
			if(order==0){

				vehc[minPos].first->prev=vehc[i].last;
				vehc[i].last->next=vehc[minPos].first;
				vehc[minPos].first=vehc[i].first;
				vehc[minPos].num+=vehc[i].num;

	
			}else{
	CEdgeHit* pn;
	for(CEdgeHit*a=vehc[i].last;a!=NULL;a=pn){
	pn=a->prev;
	vehc[minPos].fill(a);
	}	
			}
	vehc[i].clear();	
	continue;
		};


	}		

		t.end();
		t.print("merge");






		v.clear();
		int sum=0;
		for(unsigned int i=0;i<vehc.size();i++){
			if(vehc[i].first==NULL) continue;
			if(vehc[i].num>minHits){
				v.push_back(vehc[i]);
			}
		}
		CLinePara lp;
		vlp.clear();
		for(unsigned int i=0;i<v.size();i++){
			v[i].checkLine(dw,pre,calLength,goodLength);
			lp.setup(&(v[i]));
			vlp.push_back(lp);
		}
		

}

//
//void fxy(int i){
//	if(v[i].num<20) return;
//	fx.clear();
//	fy.clear();
//	double posx;
//	double posy;
//	double peakx;
//	double peaky;
//
//	for(CEdgeHit*a=v[i].first;a!=NULL;a=a->next){
//		peakx=pp.peakX(a->x,a->y,posx);
//		peaky=pp.peakY(a->x,a->y,posy);
//
//		if(peakx>peaky){
//			fx.push_back(posx);
//			fy.push_back(a->y);
//		}else{
//			fx.push_back(a->x);
//			fy.push_back(posy);		
//		}
//	}
//}
//
//void printFxy(CFrame *fr=NULL){
//
//	for(unsigned int i=0;i<fx.size();i++){
//		if(fr==NULL){
//		cout<<"index= "<<i<<"   fx="<<fx[i]<<"  fy="<<fy[i]<<"  "<<endl;
//		}else{
//	*(fr->pb(fx[i],fy[i]))=0;
//	*(fr->pg(fx[i],fy[i]))=255;
//	*(fr->pr(fx[i],fy[i]))=0;
//		}
//	}
//
//}

void ccs(CEdgeHit *a,CEdgeHit *b,CEdgeHit* pa,CEdgeHit* pb){

	double s1=atan2(double(a->y-b->y),double(a->x-b->x));
	double s2=atan2(double(pa->y-pb->y),double(pa->x-pb->x));
	c=cos(s2-s1);
	s=sin(s2-s1);

}

void expect(CEdgeHit * si,CEdgeHit * sk, CEdgeHit * mj){
	x=(sk->x-si->x)*c-(sk->y-si->y)*s+mj->x;//c=cos(s2-s1);
	y=(sk->x-si->x)*s+(sk->y-si->y)*c+mj->y;//s=sin(s2-s1);
}


CPicToEdge pe;
CPointToNorm pn;
CNormToPeak np;
CFrame fr;
CFrame frNewEdge;
IplImage* newEdge;
vector <CEdgeHitChain> vehc;
vector <CEdgeHitChain> v;
vector <CLinePara> vlp;
CEdgeHitChain ehc;
CTime t;
IplImage* f;

CPointToPeak pp;
double  c;
double  s;
double  x;
double  y;
};

