#include"lib.h"
//初略查找出轮廓的角点坐标


class FindCorner
{
//!  a constructor of FindCorner
public:FindCorner(){};

///@brief   the setup function
/** the setup function is to get all points of contour of CvSeq* tempContour
*/
///@param   input the contour of CvSeq* tempContour 
///@return  vector<CvPoint> points
void setup(CvSeq* tempContour)
{
	points.clear();
	//cout<<"tempContour"<<tempContour->total<<endl;
	CvPoint* point = new CvPoint[tempContour->total ];///定义一个 CvPoint*的数组  tempContour->total  为轮廓上的所有点数
    
	for (int i = 0; i < tempContour->total ; i++) ///提取出轮廓上的所有像素点坐标
	{
       point[i]=*CV_GET_SEQ_ELEM(CvPoint,tempContour,i);///在当前contour下一个一个的读取数据
	   points.push_back(point[i]);///轮廓上所有像素点坐标都存储进 vector<CvPoint> points
	}

		delete point;
	//	return points;///return all points of contour

}


/** 利用欧氏距离来确定左右像素点
*具有旋转稳定性
*/
///\@brief   getcorner函数获取图片全部角点坐标
///\@param   vector<CvPoint> point 输入全部像素点坐标参数
///\@param   int T1 输入为两点像素之间的模的距离需要大于的阈值---单位欧式距离
///\@param   int T2 输入为两点像素之间的模的距离需要大于的阈值---单位欧式距离
///\@param   int T3 输入为两点像素之间的模的距离需要小于等于的阈值---单位欧式距离
void getcorner(vector<CvPoint> point,int T1,int T2,int T3)	
{
	
	///清空 corner 集合
	corner.clear();
	accos.clear();
	/// 双重for 循环为查找出当前遍历的点point[i]和距离L距离的point[i-left]  和point[i-right] 
	for(int i=0;i<point.size();i++)
	{		
		 vector<CvPoint> dispoint;///   vector<CvPoint> dispoint 用来储存进行余弦定理求夹角的三个像素点坐标		
		 vector<double> disL;///   vector<double> disL 用来储存进行余弦定理求夹角的三个像素点之间的距离		 
		 dispoint.clear();///清空vector<CvPoint>  dispoint 集合		
		 disL.clear(); ///清空vector<double>  dispoint 集合
		 double max1=0;double max2=0;double k1=0;double k2=0;///定义变量


		for (int j=point.size()-1; j>=0; j--)///for循环 像素点的末端来进行与正在遍历的像素点的两者之间的欧氏距离
			{		 			
				double diss=((point[i].x-point[j].x)*(point[i].x-point[j].x)+(point[i].y-point[j].y)*(point[i].y-point[j].y));///double diss 为两像素点之间的距离  x*x+y*y
			
				double L=sqrt(diss);///double L 为两像素之间的欧式距离 sqrt(x*x+y*y)

					if((L>T1)&&(L<=T2))///if 条件来进行筛选  选出两像素之间的欧式距离在 T1<L<=T2的两像素点 
					{
						disL.push_back(L);///符合条件的欧氏距离存储在vector<double> disL
						dispoint.push_back(point[j]);///符合条件的像素点坐标存储在 vector<CvPoint> dispoint		
					}
			 }
	
	    

		 if(dispoint.size()!=0)///(dispoint.size()!=0  说明有角点
		{
			for(int k=0;k<dispoint.size();k++)///遍历候选像素点求取出每个像素点的夹角
			{  
				if((abs(point[point.size()-1].x-point[0].x)<2)||(abs(point[point.size()-1].y-point[0].y)<2))///表示闭合轮廓
				{
					
					if(((dispoint[k].x<point[i].x)||(dispoint[k].y<point[i].y)))///选取出左边距离正在遍历像素点point[i]为T1-T2区域内L值最大的点 			
					{
						if((disL[k]>=max1))
								{max1=disL[k];k1=k;}///K1存储着左边或上方L最大值的角标
			
					 }
					
					if(((dispoint[k].x>point[i].x)||(dispoint[k].y>point[i].y)))///选取出右边距离正在遍历像素点point[i]为T1-T2区域内L值最大的点 
					{
						if((disL[k]>=max2))
						        {max2=disL[k];k2=k;}///K2存储着右边或下方L最大值的角标				
					}
				
				}
			else//非闭合轮廓
			{
				if((dispoint[k].x<point[i].x)||(dispoint[k].y<point[i].y))///选取出左边距离正在遍历像素点point[i]为T1-T2区域内L值最大的点 
				{
					if(disL[k]>=max1)
							{max1=disL[k];k1=k;}
			
				 }
				if((dispoint[k].x>=point[i].x)||(dispoint[k].y>=point[i].y))///选取出右边距离正在遍历像素点point[i]为T1-T2区域内L值最大的点 
				{
					if(disL[k]>=max2)
					{max2=disL[k];k2=k;}				
				}
			}
		}
					
		     ///调用getcos函数求取出像素点的夹角
			  double ccos=getcos(point[i],dispoint[k2],dispoint[k1]);
			  
	         if(abs((double)ccos-90)<=T3){accos.push_back(ccos);corner.push_back(point[i]);}///如果像素点的夹角在与90度阈值范围内就判定为候选角点
	 }
	}
}

///\brief   getcos函数获取图片角点坐标的返余弦值--即夹角大小
///\param   CvPoint point1 输入像素点坐标参数--正在遍历的像素点point[i]坐标
///\param   CvPoint point2 输入像素点坐标参数--正在遍历的像素点point[i+L]坐标 
///\param   CvPoint point3 输入像素点坐标参数--正在遍历的像素点point[i-L]坐标 
///\return     double sita  返回正在遍历的点point[i]的反余弦值
double   getcos(CvPoint point1,CvPoint point2,CvPoint point3)// point1--i     point2--point[i+L]    point3--point[i-L]
{
	/// double distance1 是point[i]与point[i+L] 之间的距离  x*x+y*y
	double distance1=(double)((point1.x-point2.x)*(point1.x-point2.x)+(point1.y-point2.y)*(point1.y-point2.y));
	/// double distance2 是point[i]与point[i-L] 之间的距离  x*x+y*y
	double distance2=(double)((point1.x-point3.x)*(point1.x-point3.x)+(point1.y-point3.y)*(point1.y-point3.y));
	/// double distance3 是point[i+L]与point[i-L] 之间的距离  x*x+y*y
	double distance3=(double)((point3.x-point2.x)*(point3.x-point2.x)+(point3.y-point2.y)*(point3.y-point2.y));
	/// double line1 是point[i]与point[i-L] 之间的欧式距离  sqrt(x*x+y*y)
	double line1=(double)sqrt(distance1);
	/// double line2 是point[i]与point[i+L] 之间的欧式距离  sqrt(x*x+y*y)
	double line2=(double)sqrt(distance2);
	/// double line3 是point[i-L]与point[i+L] 之间的欧式距离  sqrt(x*x+y*y)
	double line3=(double)sqrt(distance3);
	///double sita   是point[i]的反余弦值
	double sita=(double)acos((double)((line1*line1)+(line2*line2)-(line3*line3))/(double)(2*line1*line2))/CV_PI*180;
    return  sita;///返回正在遍历的点point[i]的反余弦值
}



///\brief   findcorner函数获取图片角点坐标
///\param CvSeq* tempContour 输入轮廓参数
///\param int T1 输入距离阈值最小值即 L>T1
///\param int T2 输入距离阈值最大值即 L<=T2
///\param int T3 输入角度阈值--若输入为0则找出为90度的角点，角度阈值设置为  abs(90-acos)<T3,角度阈值的设定
void findcorner(IplImage* dst,CvSeq* tempContour,int T1,int T2,int T3)
{    
	 setup( tempContour);///vector<CvPoint> point 集合用来存储轮廓所有的点坐标
	 getcorner(points,T1,T2,T3);///输入point点坐标集合   int T1，T2为两点像素之间的模的距离的阈值---单位欧式距离
	  	 
}

public:
	vector<CvPoint> corner;///存储着所有角点的坐标
	vector<double> accos;///存储着所有角点的反余弦值
    vector<CvPoint>  points;//轮廓上所有像素点



};