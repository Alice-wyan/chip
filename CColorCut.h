#pragma once

#include "lib.h"

class CColorCut
{
	private:
		int num;

	protected:
		double mean[3];
		double rms[3];	
		double scaleU;
		double scaleD;

	public:
		CFrame fsrc;
		CFrame fref;

		CColorCut();
		void clear();		
		void setsigma(double *MEAN,double *RMS);		
		void setscale(double scaleup,double scaledown);	
		void setcolor(IplImage *psrc,IplImage *pref);		
		void print();	
		int cut(unsigned char* a);
		void RefreshColorCut(CfChannel refresh_mean);
		int InitalColor(char *src_name,char *ref_name,double scaleup,double scaledown);
};



	CColorCut::CColorCut(){};

	void CColorCut::clear()
	{
		mean[0]=0;mean[1]=0;mean[2]=0;
		rms[0]=0;rms[1]=0;rms[2]=0;
		num=0;
	}
	
	void CColorCut::setsigma(double *MEAN,double *RMS)
	{
		mean[0]=MEAN[0];mean[1]=MEAN[1];mean[2]=MEAN[2];
		rms[0]=RMS[0];rms[1]=RMS[1];rms[2]=RMS[2];
	}		
	
	void CColorCut::setscale(double scaleup,double scaledown)
	{
		scaleU=scaleup;
		scaleD=scaledown;	
	}		

	void CColorCut::setcolor(IplImage *psrc,IplImage *pref)
	{
		fsrc.set(psrc);
		fref.set(pref);
		for(int i=0;i<fsrc.width;i++){
			for(int j=0;j<fsrc.height;j++){
				if(fsrc.b(i,j)==fref.b(i,j)&&
					fsrc.g(i,j)==fref.g(i,j)&&
					fsrc.r(i,j)==fref.r(i,j)) continue;
				num++;
				mean[0]+=int(fsrc.b(i,j));
				mean[1]+=int(fsrc.g(i,j));
				mean[2]+=int(fsrc.r(i,j));
				rms[0]+=int(fsrc.b(i,j))*int(fsrc.b(i,j));
				rms[1]+=int(fsrc.g(i,j))*int(fsrc.g(i,j));
				rms[2]+=int(fsrc.r(i,j))*int(fsrc.r(i,j));
			}	
		}			
						
		mean[0]/=num;
		mean[1]/=num;
		mean[2]/=num;
		rms[0]=sqrt(rms[0]/num-mean[0]*mean[0]);
		rms[1]=sqrt(rms[1]/num-mean[1]*mean[1]);
		rms[2]=sqrt(rms[2]/num-mean[2]*mean[2]);
		scaleU=1;
		scaleD=1;
	};
	
	void CColorCut::print()
	{
		cout<<"mean "<<mean[0]<<" "<<mean[1]<<" "<<mean[2]<<endl;
		cout<<"rms "<<rms[0]<<" "<<rms[1]<<" "<<rms[2]<<endl;
		cout<<"scale "<<scaleU<<" "<<scaleD<<endl;
		cout<<"color num "<<num<<endl;	
	};
	
	int CColorCut::cut(unsigned char* a)//传入的是像素点Blue那个通道的指针，依次是a[0],a[1],a[2]对应BGR通道
	{
		if(int(a[0])>mean[0]+scaleU*rms[0]) return 0;if(int(a[0])<mean[0]-scaleD*rms[0]) return 0;
		if(int(a[1])>mean[1]+scaleU*rms[1]) return 0;if(int(a[1])<mean[1]-scaleD*rms[1]) return 0;
		if(int(a[2])>mean[2]+scaleU*rms[2]) return 0;if(int(a[2])<mean[2]-scaleD*rms[2]) return 0;
		return 1;
	}

	void CColorCut::RefreshColorCut(CfChannel refresh_mean)
	{
		mean[0] = refresh_mean.b;
		mean[1] = refresh_mean.g;
		mean[2] = refresh_mean.r;
	}

	int CColorCut::InitalColor(char *src_name,char *ref_name,double scaleup,double scaledown)
	{
		//char ref_name[200]="ref.bmp";	
		//char src_name[200]="src.bmp";
		char pic_name[200]="pic.bmp";	
	
		IplImage *ppic = cvLoadImage(pic_name, 1);
		IplImage *psrc = cvLoadImage(src_name, 1);
		IplImage *pref = cvLoadImage(ref_name, 1);

		if(!ppic)
		{
			cout<<"Take Picture!"<<endl;
			return 0;
		}

		if(!psrc || !pref)
		{
			cout<<"Please Mark The Picture!"<<endl;
			exit(1);//异常退出
		}

		clear();
		setcolor(psrc,pref);
		setscale(scaleup,scaledown);	
		//print();

		cvReleaseImage(&psrc);
		cvReleaseImage(&pref);
		return 1;
	}