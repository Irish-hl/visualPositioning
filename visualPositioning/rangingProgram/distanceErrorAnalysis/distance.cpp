#include<iostream>
#include<math.h>
#include<fstream>
#include<cstdlib>
double PI=3.1415926;
using namespace std;
int main()
{
	double H=17;//摄像机高度					
	double y1=13;//摄像机垂直视角投影在地面的最近距离					
	double y12=159;//摄像机垂直视角投影在地面的最远距离					
	double a=52.6/180*PI;//摄像机垂直视角与地平面y轴的夹角的最大值					
	double b=6.1/180*PI;//摄像机垂直视角与地平面y轴的夹角的最小值
	double c=47.5/180*PI;//摄像机水平视角在地面上的投影与y轴的夹角					
	double Sx=320;//成像平面横向和纵向像素数量x值
	double Sy=240;//成像平面横向和纵向像素数量y值

	double L;//所求目标点与摄像机在地面投影之间的距离					
	double angle;//p点与小车前进方向所成夹角					
	double u,v;//像素点在图片中的坐标
	double ce_x,ce_y;//特征点在摄像机坐标系中的坐标
	double real_x,real_y;//特征点在摄像机坐标系中的实际坐标
	double ce_a,ce_b;//摄像机在世界坐标系中的坐标					
	double real_a,real_b;//摄像机在世界坐标系中的实际坐标
	double L1,num=0,jun=0;

	ifstream fin("数据1.txt",ios::in);
	  if(!fin){cerr<<"无法打开"<<endl;exit(1);}
	   while(fin>>real_a>>real_b>>real_x>>real_y>>u>>v)
	   {
		//cout<<"real_a  real_b  real_x  real_y  u  v:"<<endl;
		//cin>>real_a>>real_b>>real_x>>real_y>>u>>v;
		ce_y=H*tan(PI/2-a/2-b/2-atan((Sy-2*v)/Sy*tan(a/2-b/2)));
		double ce_x_1=ce_y*tan((Sx/2-u)/(Sx/2)*c);
		double ce_x_2=ce_y*tan((Sx/2-u)/v);
		double ce_x_3=u*sqrt(y1*y1+H*H)/sqrt(H*H+y1*y1+u*u-cos(PI/2-a/2+b/2)*2*u*sqrt(H*H+y1*y1));
		ce_x=(ce_x_1+ce_x_2+ce_x_3)/3;
		ce_a=real_x-ce_x;
		ce_b=real_y-ce_y;
		//cout<<"实际("<<real_a<<","<<real_b<<")  "<<"测量("<<ce_a<<","<<ce_b<<")"<<endl;
		L=sqrt(ce_x*ce_x+ce_y*ce_y);
		L1=sqrt((real_a-real_x)*(real_a-real_x)+(real_b-real_y)*(real_b-real_y));
		cout<<"L="<<L<<"    ,实际"<<L1<<"     误差："<<L1-L<<"     误差率："<<(L1-L)/L1<<endl;
		if(((L1-L)/L1)<0.3&&((L1-L)/L1)>0)
		{	jun+=(L1-L)/L1;
			num++;
		}
			if(((L1-L)/L1)<0)
		{	jun-=(L1-L)/L1;
			num++;
		}
	   }
	   cout<<"平均误差率："<<jun/num<<endl;
	return 0;
}
