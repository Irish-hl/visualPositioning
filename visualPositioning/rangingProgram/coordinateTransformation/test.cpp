#include<iostream>
#include<cstdio>
#include<iomanip>
#include<windows.h>
#include<math.h>
#include<fstream>
#include<cstdlib>
double PI=3.1415926;
double highRate=0.20;
double lowRate=0.05;
using namespace std;
int main()
{
	double H=17;//摄像机高度					
	double y1=13;//摄像机垂直视角投影在地面的最近距离					
	double y12=159;//摄像机垂直视角投影在地面的最远距离					
	double a=52.6/180*PI;//摄像机垂直视角与地平面y轴的夹角的最大值					
	double b=6.1/180*PI;//摄像机垂直视角与地平面y轴的夹角的最小值
	double c=85.0/180*PI;//摄像机水平视角在地面上的投影与y轴的夹角					
	double Sx=320;//成像平面横向和纵向像素数量x值
	double Sy=240;//成像平面横向和纵向像素数量y值

	double L;//所求目标点与摄像机在地面投影之间的距离					
	//double angle;//p点与小车前进方向所成夹角					
	double u,v;//像素点在图片中的坐标
	double ce_x,ce_y;//特征点在摄像机坐标系中的坐标
	double real_x,real_y;//特征点在摄像机坐标系中的实际坐标
	//double ce_a,ce_b;//摄像机在世界坐标系中的坐标					
	double real_a,real_b;//摄像机在世界坐标系中的实际坐标
	double L1,num=0,jun=0;

	ifstream fin("数据1.txt",ios::in);
	if(!fin){cerr<<"无法打开"<<endl;exit(1);}
		while(fin>>real_a>>real_b>>real_x>>real_y>>u>>v)
		{
			//ce_y=H*tan((PI/2-a+(a-b)/2+(v-Sy/2)*(a-b)/2/(Sy/2)));
			//修正比例模型的ce_y
			ce_y=H*tan((PI/2-a+(a-b)/2+atan((v-Sy/2)/(Sy/2)*tan((a-b)/2))));
			double ce_x_1=ce_y*tan((Sx/2-u)/(Sx/2)*c);
			double ce_x_2=ce_y*tan((Sx/2-u)/v);
			double ce_x_3=u*sqrt(y1*y1+H*H)/sqrt(H*H+y1*y1+u*u-cos(PI/2-a/2+b/2)*2*u*sqrt(H*H+y1*y1));
			double ce_x_4=ce_y*tan((u-Sx/2)*c/(Sx/2));
			//修正比例模型的ce_x
			double ce_x_5=ce_y*tan(atan((u-Sx/2)/(Sx/2)*tan(c/2)));
			ce_x=(ce_x_1+ce_x_2+ce_x_3+ce_x_4+ce_x_5)/3;
			//摄像机坐标系中目标特征点到原点的距离
			L=sqrt(ce_x_5*ce_x_5+ce_y*ce_y);
			//世界坐标系中目标特征点到摄像机坐标系原点的距离
			L1=sqrt((real_a-real_x)*(real_a-real_x)+(real_b-real_y)*(real_b-real_y));
			//将误差率（-rate，rate）的数据红色显示
			if(fabs((L1-L)/L1)>=highRate)
			{
				SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),FOREGROUND_RED);
				cout<<"摄像机("<<setprecision(2)<<real_a<<","<<setprecision(2)<<real_b<<"),  特征点（"<<setprecision(2)<<real_x<<","<<setprecision(2)<<real_y<<"),  像素点("<<setprecision(3)<<u<<","<<setprecision(3)<<v<<"),  ""测得距离："<<setprecision(5)<<L<<",  实际距离："<<setprecision(5)<<L1<<",  误差："<<setprecision(5)<<L1-L<<",  误差率："<<setprecision(5)<<(L1-L)/L1<<endl;	
				SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),7);
			}
			else if(fabs((L1-L)/L1)<=lowRate)
			{
				SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),FOREGROUND_GREEN);
				cout<<"摄像机("<<setprecision(2)<<real_a<<","<<setprecision(2)<<real_b<<"),  特征点（"<<setprecision(2)<<real_x<<","<<setprecision(2)<<real_y<<"),  像素点("<<setprecision(3)<<u<<","<<setprecision(3)<<v<<"),  ""测得距离："<<setprecision(5)<<L<<",  实际距离："<<setprecision(5)<<L1<<",  误差："<<setprecision(5)<<L1-L<<",  误差率："<<setprecision(5)<<(L1-L)/L1<<endl;
				SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),7);
			}
			else
				cout<<"摄像机("<<setprecision(2)<<real_a<<","<<setprecision(2)<<real_b<<"),  特征点（"<<setprecision(2)<<real_x<<","<<setprecision(2)<<real_y<<"),  像素点("<<setprecision(3)<<u<<","<<setprecision(3)<<v<<"),  ""测得距离："<<setprecision(5)<<L<<",  实际距离："<<setprecision(5)<<L1<<",  误差："<<setprecision(5)<<L1-L<<",  误差率："<<setprecision(5)<<(L1-L)/L1<<endl;
			//统计误差率位于（-rate，rate）之间的误差率绝对值平均值
			if(fabs((L1-L)/L1)<highRate)
			{	
				jun+=fabs(L1-L)/L1;
				num++;
			}
		}
	    cout<<"平均误差率："<<jun/num<<endl;
		fin.close();
	return 0;
}