#include"getcorner.h"
#include"lib.h"
#include"getrectframe.h"
#include"GetCameraNoise.h"
#include"testdetect.h"
#include"testdetectpoint.h"
#include"testcorner.h"
#include"testnoise.h"
#include"getdecorner.h"
#include"linetodetect.h"




///
///������
///
int  main()
{
	///
	///������Ҫ�����������
	///
	GetCameraNoise gn;
	Getcorner gc;
	TestDetect  td;
	TestDetectPoint tp;
	DetectLine dl;
	Testcorner tc;
	CommomFunction cm;
	TestDetecteCorner tdc;
	TestNoise  tdt;

	LineToDetect  lt;


	CvCapture* capture1 = cvCreateCameraCapture( 0 ); ///����CvCapture ��������ͷ����˿�

	IplImage *src = cvQueryFrame(capture1 );///������ͷ����ͼƬ����IplImage ָ�� src

	//IplImage *src=cvLoadImage("3.png");

	///
	///������GetCameraNoise��ʱʹ������ռ� 
	///
	//int row=src->width;//ͼƬ����
	//int com=src->height;
	//double ** a = new double *[row];//���嶯̬����洢ͼƬ��ά�������ص�garyֵ
	//for(int i = 0;i < row;i++)
	//{a[i] = new double[com];}//��ʼ������
	//	double ** b = new double *[row];//���嶯̬����洢ͼƬ��ά�������ص�garyֵ
	//for(int i = 0;i < row;i++)
	//{b[i] = new double[com];}//��ʼ������
	//for(int i = 0;i < row;i++){  memset(a[i], 0, sizeof(double)*com);}
	//for(int i = 0;i < row;i++){  memset(b[i], 0, sizeof(double)*com);}



	//cvNamedWindow( "src", CV_WINDOW_AUTOSIZE );	 ///define a window
	//cvResizeWindow( "src", 320, 320 );///���ô��ڴ�С
	///
	///����i�����òɼ�100֡ͼƬ���˳�����
    ///
	int i=0;
	///
	///����sec��dst dst������Ϊ�˻���̽���߽�ֱ�۹۲�
	///
	//IplImage *dst=cvCloneImage(src);

	///
	///����ѭ��һֱ�ɼ�����ͷͼƬ
	///
	while(1)
	{
		   ///
		  ///����CvCapture ��������ͷ����˿�
		  ///
		CvCapture* capture1 = cvCreateCameraCapture( 0 );
		
		//������ͷ����ͼƬ����IplImage ָ�� src
		
		IplImage *src = cvQueryFrame( capture1 );

		if((src!=0))//���src��Ϊ����ִ�����³���'
		{  
			
			IplImage *dst=cvCloneImage(src);
		//if(i>=10)
		{	
			CvPoint p1,p2,p3;
			p1.x=281;p1.y=245;
			p2.x=317;p2.y=315;
			p3.x=333;p3.y=347;
			//tdc.getcorner(src,dst);
			//tc.getcorner(p1,p2,p3,src,dst);

		      gc.getcornerpoint(src,dst,8,9,30,10);//��ȡ�����Ͼ�ȷ�ǵ�����
		    // gn.getallgary(src,a,b);//��ȡ����ͷrgb����
		    ///
			///����TestDetect���е�getdetectedge��������ȡ��Ե����
			///
		    //tp.getdetectedge(src,dst);
			//tdt.getnoise();

		}
		cvShowImage( "src", src );///  show  the  IplImage *dst
		cvShowImage( "dst", dst );
		i++;
		cout<<"i="<<i<<endl;
		if(i==110){break;}
		char c=cvWaitKey(33);
		 cvReleaseImage(&dst);
		 cvReleaseImage(&src);
		 if(c==27)break;
		}
	}

	cvWaitKey(27);

	//double l1x=cm.getrmse(tdc.l1x);
	//cout<<"l1x=="<<l1x<<endl;
	//double l1y=cm.getrmse(tdc.l1y);
	//cout<<"l1y=="<<l1y<<endl;

	//double l2x=cm.getrmse(tdc.l2x);
	//cout<<"l2x=="<<l2x<<endl;
	//double l2y=cm.getrmse(tdc.l2y);
	//cout<<"l2y=="<<l2y<<endl;


	//double l3x=cm.getrmse(tdc.l3x);
	//cout<<"l3x=="<<l3x<<endl;
	//double l3y=cm.getrmse(tdc.l3y);
	//cout<<"l3y=="<<l3y<<endl;

	//double l4x=cm.getrmse(tdc.l4x);
	//cout<<"l4x=="<<l4x<<endl;
	//double l4y=cm.getrmse(tdc.l4y);
	//cout<<"l4y=="<<l4y<<endl;

	//double x=cm.getrmse(tdc.X);
	//cout<<"x=="<<x<<endl;
	//double y=cm.getrmse(tdc.Y);
	//cout<<"y=="<<y<<endl;




	//gn.getgarynoise(dst,a,b);
	//double xx=cm.getrmse(tdt.xpoint);
	//cout<<"XXXXXX.size()="<<tc.X.size()<<endl;
		//double xx=cm.getrmse(tc.X);
		//cout<<"rmsx="<<xx<<endl;
		//double yy=cm.getrmse(tdt.ypoint);
		//cout<<"rmsy="<<yy<<endl;
     // cout<<"YYYYYY.size()="<<tc.Y.size()<<endl;
		//double yy=cm.getrmse(tc.Y);
		//cout<<"rmsy="<<yy<<endl;
	/////
	/////����TestDetect���е�testdetectpointx��������ȡ�����rmseƽ��ֵ
	/////
	//cout<<"xxx==="<<endl;
	//td.calcpointrms(tc.px,tc.xpoint);//����̽�������ܲ���
	//td.calcpointrms(tp.px,tp.xdetectpoint);//����̽�������ܲ���
	//cout<<"yyy==="<<endl;
	/////
	/////����TestDetect���е�testdetectpointx��������ȡ�����rmseƽ��ֵ
	/////
	//td.calcpointrms(tc.py,tc.ypoint);
	//td.calcpointrms(tp.py,tp.ydetectpoint);
	//	delete a;//�ͷŶ�̬�����ڴ�
	//	delete b;//�ͷŶ�̬�����ڴ�

		return 0;

}
		