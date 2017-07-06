#include"lib.h"
//���Բ��ҳ������Ľǵ�����


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
	CvPoint* point = new CvPoint[tempContour->total ];///����һ�� CvPoint*������  tempContour->total  Ϊ�����ϵ����е���
    
	for (int i = 0; i < tempContour->total ; i++) ///��ȡ�������ϵ��������ص�����
	{
       point[i]=*CV_GET_SEQ_ELEM(CvPoint,tempContour,i);///�ڵ�ǰcontour��һ��һ���Ķ�ȡ����
	   points.push_back(point[i]);///�������������ص����궼�洢�� vector<CvPoint> points
	}

		delete point;
	//	return points;///return all points of contour

}


/** ����ŷ�Ͼ�����ȷ���������ص�
*������ת�ȶ���
*/
///\@brief   getcorner������ȡͼƬȫ���ǵ�����
///\@param   vector<CvPoint> point ����ȫ�����ص��������
///\@param   int T1 ����Ϊ��������֮���ģ�ľ�����Ҫ���ڵ���ֵ---��λŷʽ����
///\@param   int T2 ����Ϊ��������֮���ģ�ľ�����Ҫ���ڵ���ֵ---��λŷʽ����
///\@param   int T3 ����Ϊ��������֮���ģ�ľ�����ҪС�ڵ��ڵ���ֵ---��λŷʽ����
void getcorner(vector<CvPoint> point,int T1,int T2,int T3)	
{
	
	///��� corner ����
	corner.clear();
	accos.clear();
	/// ˫��for ѭ��Ϊ���ҳ���ǰ�����ĵ�point[i]�;���L�����point[i-left]  ��point[i-right] 
	for(int i=0;i<point.size();i++)
	{		
		 vector<CvPoint> dispoint;///   vector<CvPoint> dispoint ��������������Ҷ�����нǵ��������ص�����		
		 vector<double> disL;///   vector<double> disL ��������������Ҷ�����нǵ��������ص�֮��ľ���		 
		 dispoint.clear();///���vector<CvPoint>  dispoint ����		
		 disL.clear(); ///���vector<double>  dispoint ����
		 double max1=0;double max2=0;double k1=0;double k2=0;///�������


		for (int j=point.size()-1; j>=0; j--)///forѭ�� ���ص��ĩ�������������ڱ��������ص������֮���ŷ�Ͼ���
			{		 			
				double diss=((point[i].x-point[j].x)*(point[i].x-point[j].x)+(point[i].y-point[j].y)*(point[i].y-point[j].y));///double diss Ϊ�����ص�֮��ľ���  x*x+y*y
			
				double L=sqrt(diss);///double L Ϊ������֮���ŷʽ���� sqrt(x*x+y*y)

					if((L>T1)&&(L<=T2))///if ����������ɸѡ  ѡ��������֮���ŷʽ������ T1<L<=T2�������ص� 
					{
						disL.push_back(L);///����������ŷ�Ͼ���洢��vector<double> disL
						dispoint.push_back(point[j]);///�������������ص�����洢�� vector<CvPoint> dispoint		
					}
			 }
	
	    

		 if(dispoint.size()!=0)///(dispoint.size()!=0  ˵���нǵ�
		{
			for(int k=0;k<dispoint.size();k++)///������ѡ���ص���ȡ��ÿ�����ص�ļн�
			{  
				if((abs(point[point.size()-1].x-point[0].x)<2)||(abs(point[point.size()-1].y-point[0].y)<2))///��ʾ�պ�����
				{
					
					if(((dispoint[k].x<point[i].x)||(dispoint[k].y<point[i].y)))///ѡȡ����߾������ڱ������ص�point[i]ΪT1-T2������Lֵ���ĵ� 			
					{
						if((disL[k]>=max1))
								{max1=disL[k];k1=k;}///K1�洢����߻��Ϸ�L���ֵ�ĽǱ�
			
					 }
					
					if(((dispoint[k].x>point[i].x)||(dispoint[k].y>point[i].y)))///ѡȡ���ұ߾������ڱ������ص�point[i]ΪT1-T2������Lֵ���ĵ� 
					{
						if((disL[k]>=max2))
						        {max2=disL[k];k2=k;}///K2�洢���ұ߻��·�L���ֵ�ĽǱ�				
					}
				
				}
			else//�Ǳպ�����
			{
				if((dispoint[k].x<point[i].x)||(dispoint[k].y<point[i].y))///ѡȡ����߾������ڱ������ص�point[i]ΪT1-T2������Lֵ���ĵ� 
				{
					if(disL[k]>=max1)
							{max1=disL[k];k1=k;}
			
				 }
				if((dispoint[k].x>=point[i].x)||(dispoint[k].y>=point[i].y))///ѡȡ���ұ߾������ڱ������ص�point[i]ΪT1-T2������Lֵ���ĵ� 
				{
					if(disL[k]>=max2)
					{max2=disL[k];k2=k;}				
				}
			}
		}
					
		     ///����getcos������ȡ�����ص�ļн�
			  double ccos=getcos(point[i],dispoint[k2],dispoint[k1]);
			  
	         if(abs((double)ccos-90)<=T3){accos.push_back(ccos);corner.push_back(point[i]);}///������ص�ļн�����90����ֵ��Χ�ھ��ж�Ϊ��ѡ�ǵ�
	 }
	}
}

///\brief   getcos������ȡͼƬ�ǵ�����ķ�����ֵ--���нǴ�С
///\param   CvPoint point1 �������ص��������--���ڱ��������ص�point[i]����
///\param   CvPoint point2 �������ص��������--���ڱ��������ص�point[i+L]���� 
///\param   CvPoint point3 �������ص��������--���ڱ��������ص�point[i-L]���� 
///\return     double sita  �������ڱ����ĵ�point[i]�ķ�����ֵ
double   getcos(CvPoint point1,CvPoint point2,CvPoint point3)// point1--i     point2--point[i+L]    point3--point[i-L]
{
	/// double distance1 ��point[i]��point[i+L] ֮��ľ���  x*x+y*y
	double distance1=(double)((point1.x-point2.x)*(point1.x-point2.x)+(point1.y-point2.y)*(point1.y-point2.y));
	/// double distance2 ��point[i]��point[i-L] ֮��ľ���  x*x+y*y
	double distance2=(double)((point1.x-point3.x)*(point1.x-point3.x)+(point1.y-point3.y)*(point1.y-point3.y));
	/// double distance3 ��point[i+L]��point[i-L] ֮��ľ���  x*x+y*y
	double distance3=(double)((point3.x-point2.x)*(point3.x-point2.x)+(point3.y-point2.y)*(point3.y-point2.y));
	/// double line1 ��point[i]��point[i-L] ֮���ŷʽ����  sqrt(x*x+y*y)
	double line1=(double)sqrt(distance1);
	/// double line2 ��point[i]��point[i+L] ֮���ŷʽ����  sqrt(x*x+y*y)
	double line2=(double)sqrt(distance2);
	/// double line3 ��point[i-L]��point[i+L] ֮���ŷʽ����  sqrt(x*x+y*y)
	double line3=(double)sqrt(distance3);
	///double sita   ��point[i]�ķ�����ֵ
	double sita=(double)acos((double)((line1*line1)+(line2*line2)-(line3*line3))/(double)(2*line1*line2))/CV_PI*180;
    return  sita;///�������ڱ����ĵ�point[i]�ķ�����ֵ
}



///\brief   findcorner������ȡͼƬ�ǵ�����
///\param CvSeq* tempContour ������������
///\param int T1 ���������ֵ��Сֵ�� L>T1
///\param int T2 ���������ֵ���ֵ�� L<=T2
///\param int T3 ����Ƕ���ֵ--������Ϊ0���ҳ�Ϊ90�ȵĽǵ㣬�Ƕ���ֵ����Ϊ  abs(90-acos)<T3,�Ƕ���ֵ���趨
void findcorner(IplImage* dst,CvSeq* tempContour,int T1,int T2,int T3)
{    
	 setup( tempContour);///vector<CvPoint> point ���������洢�������еĵ�����
	 getcorner(points,T1,T2,T3);///����point�����꼯��   int T1��T2Ϊ��������֮���ģ�ľ������ֵ---��λŷʽ����
	  	 
}

public:
	vector<CvPoint> corner;///�洢�����нǵ������
	vector<double> accos;///�洢�����нǵ�ķ�����ֵ
    vector<CvPoint>  points;//�������������ص�



};