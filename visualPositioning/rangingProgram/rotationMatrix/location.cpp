#include <iostream>
#include <iomanip>
#include<math.h>
using namespace std;

int main()
{
	double p_cx,p_cy;
	double p_wx,p_wy;
	double q_cx,q_cy;
	double q_wx,q_wy;
	double 
	double cos,sin,tx,ty;

	cout<<"请输入第一个点的摄像机坐标系下的坐标"<<endl;
	cin>>p_cx>>p_cy;
	cout<<"请输入第一个点的世界坐标系下的坐标"<<endl;
	cin>>p_wx>>p_wy;
	cout<<"请输入第二个点的摄像机坐标系下的坐标"<<endl;
	cin>>q_cx>>q_cy;
	cout<<"请输入第二个点的世界坐标系下的坐标"<<endl;
	cin>>q_wx>>q_wy;
	//(p_wx-q_wx)*cos-(p_wy-q_wy)*sin==(p_wx-q_wx)*sin+(p_wy-q_wy)*cos;
	sin=sqrt(1-cos*cos);
	cos=(p_wx-q_wx+p_wy-q_wy)/(p_wx-q_wx-p_wy+q_wy)*sin;
	cout<<"旋转矩阵为:"<<endl<<cos<<" "<<-sin<<endl<<sin<<" "<<cos<<endl;
}





