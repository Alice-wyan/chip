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
///主函数
///
int  main()
{
	///
	///定义需要被调用类对象
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


	CvCapture* capture1 = cvCreateCameraCapture( 0 ); ///定义CvCapture 对象即摄像头输入端口

	IplImage *src = cvQueryFrame(capture1 );///把摄像头输入图片传给IplImage 指针 src

	//IplImage *src=cvLoadImage("3.png");

	///
	///给调用GetCameraNoise类时使用申请空间 
	///
	//int row=src->width;//图片长宽
	//int com=src->height;
	//double ** a = new double *[row];//定义动态数组存储图片二维数组像素点gary值
	//for(int i = 0;i < row;i++)
	//{a[i] = new double[com];}//初始化数组
	//	double ** b = new double *[row];//定义动态数组存储图片二维数组像素点gary值
	//for(int i = 0;i < row;i++)
	//{b[i] = new double[com];}//初始化数组
	//for(int i = 0;i < row;i++){  memset(a[i], 0, sizeof(double)*com);}
	//for(int i = 0;i < row;i++){  memset(b[i], 0, sizeof(double)*com);}



	//cvNamedWindow( "src", CV_WINDOW_AUTOSIZE );	 ///define a window
	//cvResizeWindow( "src", 320, 320 );///设置窗口大小
	///
	///定义i有利用采集100帧图片后退出程序
    ///
	int i=0;
	///
	///复制sec给dst dst作用是为了画出探测线较直观观测
	///
	//IplImage *dst=cvCloneImage(src);

	///
	///进入循环一直采集摄像头图片
	///
	while(1)
	{
		   ///
		  ///定义CvCapture 对象即摄像头输入端口
		  ///
		CvCapture* capture1 = cvCreateCameraCapture( 0 );
		
		//把摄像头输入图片传给IplImage 指针 src
		
		IplImage *src = cvQueryFrame( capture1 );

		if((src!=0))//如果src不为空则执行以下程序'
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

		      gc.getcornerpoint(src,dst,8,9,30,10);//获取轮廓上精确角点坐标
		    // gn.getallgary(src,a,b);//获取摄像头rgb噪声
		    ///
			///调用TestDetect类中的getdetectedge函数来获取边缘坐标
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
	/////调用TestDetect类中的testdetectpointx函数来获取坐标的rmse平均值
	/////
	//cout<<"xxx==="<<endl;
	//td.calcpointrms(tc.px,tc.xpoint);//进行探测线性能测试
	//td.calcpointrms(tp.px,tp.xdetectpoint);//进行探测线性能测试
	//cout<<"yyy==="<<endl;
	/////
	/////调用TestDetect类中的testdetectpointx函数来获取坐标的rmse平均值
	/////
	//td.calcpointrms(tc.py,tc.ypoint);
	//td.calcpointrms(tp.py,tp.ydetectpoint);
	//	delete a;//释放动态数组内存
	//	delete b;//释放动态数组内存

		return 0;

}
		