#ifndef _matrix_h
#define _matrix_h
#include <iostream>
#include <iomanip>
#include <cmath>
#define eps (1e-06)
#define KMAX 60  //最大迭代次数
using namespace std;
int cifang(int i,int j); //全局函数：求-1的(i+j)次方
int chartoint(const char* p);//字符转为整数型
int getRank(double *Array,int row,int col);//全局函数：求秩
int IsZero(double x){ return fabs(x) <= 0.00001?1:0;} //全局函数：判断是否为零
int jacobieigen(double *a,double *u,int jie);//a 和u 都是jie*jie的方阵
void getEigen(double *Array,int dim);//全局函数：获得特征值及特征向量。Array是dim*dim方阵
void show(double array[],int jie); //全局函数：显示矩阵
void transpose(double array[],int jie); //全局函数：转置
bool IsSymmetricMatrix(double *Array, int dim);//判断是对称矩阵
bool notfindchar(const char* p,int index);
bool Is0Array( double *Array,int row,int col); //全局函数：判断是否是0矩阵
bool chardouble(const char *p);//全局函数：判断字符是否满足转换要求
double det(double array[],int Jie); //全局函数：求行列式
double chartodouble(const char* p);//全局函数：字符型转double型
double* doInverse(double Array[],int row);//全局函数：求逆

int cifang(int i,int j) //全局函数：求-1的(i+j)次方
{ 
	if( i < 0 || j < 0 )
	{
		cerr<<"i,j非法!" << endl;
		exit(-1);
	}
	int temp = i+j;
	if( temp%2 == 0)
		return 1;
	else
		return -1;
	}
int chartoint(const char* p)//字符转为整数型
{
	if( p == NULL)
		return 0;
	int length = strlen(p);
	int INT = 0,temp =0;
	if( p[0] != '-' )
	{
		for( int i = 0; i < length; i++)
		{
			temp = p[i] - '0';
			INT = INT*10 + temp;
		}
		return INT;
	}
	else
	{
		for( int j = 1; j < length; j++)
		{
			temp = p[j] - '0';
			INT = INT*10 + temp;
		}
		return (-INT);
	}
}
int getRank(double *Array,int row,int col) //全局函数：求秩
{
	int i,j,k,l,i1,j1,main_row,main_col,rank;
	double main_element,temp;
	for( i = 0; i < row; i++ ) //行循环
	{ 
		main_element = *(Array + i*col + i );//对角元素
		main_row = i;//主元坐标：主对角元(i,i)
		main_col=i;
		for( j = i; j < row; j++ ) //循环：寻找当前(i,i)元之后的最大元，并记为坐标(main_row,mai_col)
			for( k = i; k < col; k++ )
			{ 
				if( fabs(*( Array + j*col + k )) >= fabs(main_element) && (i != j) )
				{ 
					main_element = *(Array + j*col + k);
					main_row=j;
					main_col=k;
				}
			}
		for(l = 0; l < col; l++ )
		{ 
			temp = *(Array + main_row*col + l);
			*(Array + main_row*col + l) = *(Array +i*col + l);
			*(Array +i*col + l) = temp;
		}
		for(l = 0; l < row; l++ )
		{ 
			temp = *(Array + l*col + main_col);
			*(Array + l*col + main_col) = *(Array + l*col + i);
			*(Array + l*col + i) = temp;
		}
		if( IsZero(main_element) == 1)
			break;
		for( i1 = i+1; i1 < row; i1++ )
		{ 
			temp = *(Array + i1*col + i);
			for( j1 = 0; j1 < col; j1++ )
				*(Array + i1*col + j1) = *(Array + i1*col +j1) - *(Array + i*col + j1)*temp/main_element;
		}
	}
	rank = 0;
	for( i = 0; i < row; i++ )
	{ 
		for( j = 0; j < col; j++ )
		{
			if( IsZero( *(Array + i*col + j)) == 0 )
			{
				rank++;
				break;
			}
			else
				continue;
		}
	}
	return rank;
}

int jacobieigen(double *a,double *u,int jie)//a 和u 都是jie*jie的方阵
{ 
	int i,j,k,p,q;
	double d,m,x,y,sn,cn,w;
	for( i = 0; i < jie; i++ ) //产生单位矩阵
	{ 
		(*(u+i*jie+i)) = 1;
		for( j = 0; j < jie; j++ )
		{ 
			if( i!=j )
				(*(u+i*jie+j)) = 0;
		}
	}
	k = 1;
	while(1)
	{
		m = 0;
		for( i = 0; i <= jie-1; i++ ) //选取绝对值最大的对角线元素
		{ 
			for( j = 0; j <= i-1; j++ )
			{ 
				d = fabs((*(a+i*jie+j)));
				if( (i!=j) && (d>m) )
				{ 
					m = d;
					p = i;
					q = j;
				}
			}
		}
		if( m < eps )  //满足精度要求，正常返回
			return(1);
		if( k > KMAX )  //超过最大迭代次数返回
			return(-1);
		k = k + 1;
		x = -(*(a+p*jie+q));
		y = ( (*(a+q*jie+q)) - (*(a+p*jie+p)) )/2.0;
		w = x/sqrt( x*x + y*y );
		if( y < 0 )
			w = -w;
		sn = 1 + sqrt( 1 - w*w );
		sn = w/sqrt( 2*sn );
		cn = sqrt( 1 - sn*sn );
		m = (*(a+p*jie+p));  //计算矩阵A的新元素
		(*(a+p*jie+p)) = m*cn*cn + (*(a+q*jie+q))*sn*sn + (*(a+p*jie+q))*w;
		(*(a+q*jie+q)) = m*sn*sn + (*(a+q*jie+q))*cn*cn - (*(a+p*jie+q))*w;
		(*(a+p*jie+q)) = 0;
		(*(a+q*jie+p)) = 0;
		for( j = 0; j < jie; j++ )
		{ 
			if( (j!=p) && (j!=q) )
			{ 
				m = (*(a+p*jie+j));
				(*(a+p*jie+j)) = m*cn + (*(a+q*jie+j))*sn;
				(*(a+q*jie+j)) = -m*sn + (*(a+q*jie+j))*cn;
			}
		}
		for( i = 0; i < jie; i++ )
		{
			if( (i!=p)&&(i!=q) )
			{ 
				m = (*(a+i*jie+p));
				(*(a+i*jie+p)) = m*cn + (*(a+i*jie+q))*sn;
				(*(a+i*jie+q)) = -m*sn + (*(a+i*jie+q))*cn;
			}
		}
		for( i = 0; i < jie; i++ )
		{ 
			m = (*(u+i*jie+p));
			(*(u+i*jie+p)) = m*cn + (*(u+i*jie+q))*sn;
			(*(u+i*jie+q)) = -m*sn + (*(u+i*jie+q))*cn;
		}
	}
}

void getEigen(double *Array,int dim)//全局函数：获得特征值及特征向量。Array是dim*dim方阵
{ 
	double *a = new double[dim*dim],*b = new double[dim*dim],*v = new double[dim*dim];
	int i,j,k;
	for( i = 0; i < dim; i++ )
		for( j = 0; j < dim; j++ )
			*(a+i*dim+j) = *(Array+i*dim+j);
	if( !IsSymmetricMatrix(a,dim) )
	{
		cerr << "不求特征值及特征向量\n";
		return; 
	}
	k = jacobieigen(a,v,dim);
	if( k == 1 )
	{ 
		cout << "\n特征值依次是\n";
		for( i = 0; i < dim; i++ )
		{ 
			cout << i+1 << ":";
			for( j = 0; j < dim; j++ )
				if( i == j )
					printf("%11f",(*(a+i*dim+j)));
			cout << endl;
		}
		cout << endl;
		for( i = 0; i < dim; i++ )
			for( j = 0; j < dim; j++ )
				*(b+i*dim+j) = *(v+i*dim+j);
		cout << "对应的特征向量\n";
		for( i = 0; i < dim; i++ )
		{ 
			printf(" %d:",i+1);
			for( j = 0; j < dim; j++ )
				printf(" %11f",*(b+i*dim+j));
			cout << endl;
		}
	}
	else
		cout << "dimoEigenValue\n";
}

void transpose(double array[],int jie) //全局函数：转置
{ 
	if( jie <= 0){cerr<< "数组非法!" << endl;exit(-1);}
	double temp;
	for( int i = 0; i < jie; i++ )
		for( int j = 0; j < i; j++ )
		{ 
			temp = *(array + i*jie + j);
			*(array + i*jie + j) = *(array + j*jie + i );
			*(array + j*jie + i) = temp;
		}
}

void show(double array[],int jie) //全局函数：方阵输出
{ 
	if(jie <= 0)
	{ 
		cerr << " 非法" << endl;
		exit(-1);
	}
	cout << "矩阵输出为: " << endl;
	int i, j;
	for( i = 0; i < jie; i++ )
	{ 
		for( j = 0; j < jie; j++ )
			cout << *(array + i*jie + j) << '\t';
		cout << endl;
	}
}
bool IsSymmetricMatrix(double *Array, int dim)//判断是对称矩阵
{ 
	if( dim <= 0 )
	{
		cerr << "参数非法,退出!" << endl;
		exit(-1);
	}
	int i,j;
	for( i = 0; i < dim; i++ )
		for( j = 0; j <= i; j++ )
		{ 
			if( *(Array+i*dim+j) != *(Array+j*dim+i) )
			{ 
				cout << "不是对称矩阵" << endl;
				return false;
			}
		}
	cout << "是对称矩阵!" << endl;
	return true;
}

bool Is0Array( double *Array,int row,int col)//全局函数：判断是否是0矩阵
{ 
	if( row <= 0 || col<= 0 )
	{
		cerr <<"行列值非法,退出!\n";
		exit(-1);
	}
	int i,j;
	double temp1 = 0,temp2 = 0;
	for( i = 0; i < row; i++ )
		for( j = 0; j < col; j++ )
		{
			temp2 = fabs(*(Array + i*col +j));
			temp1 += temp2;
		}
	if(temp1 < 0.00001)
	{
		cout << "是" << row << '*' << col << "的0矩阵\n";
		return true;
	}
	return false;
}

bool notfindchar(const char* p,int index)
{
	if( p == NULL )
		return true;
	int length = strlen(p);
	for( int i = index; i <= length-1; i++ )
		if(p[i]<'0' || p[i] >'9')
			return false;
	return true;
}

double det(double array[],int Jie) //全局函数：求行列式
{ 
	if( Jie <= 0 )
	{
		cerr << "阶小于0或等于0!" << endl;
		return 0;
	}
	else if( Jie == 1)
		return array[0];
	else
	{ 
		int i,j,k,tag;
		double *subArray[500];  //须填固定值?
		for( i = 0; i < Jie; i++ )
			subArray[i] = new double[(Jie-1)*(Jie-1)];
		for( i = 0; i < Jie; i++ )
			for( j = 0; j < Jie-1; j++ )
				for( k = 0; k < Jie-1; k++ )
					*(subArray[i] + j*(Jie-1) + k) = 0;
		for( i = 0; i < Jie; i++ )
			for( j = 0; j < Jie-1; j++ )
				for( k = 0; k < Jie-1; k++ )
				{ 
					if( k < i )
						*(subArray[i] + j*(Jie-1) + k) = *(array + (j+1)*Jie + k );
					else
						*(subArray[i] + j*(Jie-1) + k) = *(array + (j+1)*Jie + k+1 );     
				}
		double temp= 0;
		tag = 1;
		for( i = 0 ; i < Jie; i++)
		{ 
			temp += tag * det(subArray[i],Jie-1) * array[i];
			tag *= -1;
		}
		return temp;
	}
}

double *doInverse(double Array[],int row)//全局函数：求逆  
{ 
	double d_det = det(Array,row);
	if(d_det == 0)
	{
		cerr << "行列式为0,不存在逆矩阵!" << endl;
		exit(-1);
	} 
	int i,j,k,h,subRow = row-1;
	double *subArray[1000];      //须填固定值?
	for( i = 0; i < 500; i++ )
	{ 
		if( i < row*row )
		{
			subArray[i] = new double[subRow*subRow];
			for( j = 0; j < subRow; j++ )
				for( k = 0; k < subRow; k++ )
					*(subArray[i] + j*subRow + k) = 0;
		}
		else
			subArray[i] = NULL;
	}
	for( i= 0; i < row; i++ )
		for( h = 0; h < row; h++ )
			for( j = 0; j < subRow; j++ )
				for( k = 0; k < subRow; k++ )
				{
					if( j < i && k < h )
						*(subArray[i*row+h] + j*subRow + k) = *(Array + j*row + k );
					if( j < i && k >= h )
						*(subArray[i*row+h] + j*subRow + k) = *(Array + j*row + k+1 );
					if( j >= i && k < h )
						*(subArray[i*row+h] + j*subRow + k) = *(Array + (j+1)*row + k );
					if( j >= i && k >= h)
						*(subArray[i*row+h] + j*subRow + k) = *(Array + (j+1)*row + k+1 );
				}
	double *tempArray = new double[row*row]; //创建临时数组，赋值后并返回
	for( i = 0; i < row*row; i++ )    //初始化
		tempArray[i] = 0; 
	for( i = 0; i < row; i++ )
		for( j = 0; j < row; j++ )
			*(tempArray + i*row + j) = cifang(i,j)*det(subArray[i*row+j],subRow)/d_det;//求逆矩阵的(i,j)元
	transpose(tempArray,row); //调用转置后的伴随
	return tempArray;
}
bool chardouble(const char *p)
{
	if( p == NULL )
		return true;//风险?
	int length=strlen(p),count1=0,count2=0,i;
	for( i = 0; i < length; i++ )
	{
		if( p[i] =='.')
			count1++;
		if( p[i] < '0' || p[i] > '9' )
			count2++;
	}
	if( count2 == 0)
		return true;
	if( count2 == 1  && p[0] == '-' && p[1]!='\0' )
		return true;
	if( count2 == 1 && count1 == 1 && p[0]!='.' )
		return true;
	if( count2 == 2 && count1 == 1 && p[0] == '-' && p[1] !='.' )
		return true;
	return false;
}

double chartodouble(const char* p)
{
	if( p == NULL )
		return 0;
	double DOUBLE = 0,DOUBLE1 = 0;
	int length=strlen(p),count1=0,count2=0,temp = 0,i,j=-1;
	for( i = 0; i < length; i++ )
	{
		if( p[i] =='.')
		{
			count1++;
			j=i;
		}
		if( p[i] < '0' || p[i] > '9' )
			count2++;
	}
	if( count2 == 0)
	{
		for( i = 0;i< length;i++)
		{
			temp = p[i]-'0';
			DOUBLE = DOUBLE*10 + temp;
		}
		return DOUBLE;
	}
	if( count2 == 1  && p[0] == '-'&& p[1]!='\0' )
	{
		for( i = 1;i < length;i++)
		{
			temp = p[i] - '0';
			DOUBLE = DOUBLE*10 + temp;
		}
		return -DOUBLE;
	}
	if( count2 == 1 && count1 == 1 && p[0]!='.' )
	{
		for( i = 0;i < j;i++)
		{
			temp = p[i] - '0';
			DOUBLE = DOUBLE*10 + temp;
		}
		for( i = length-1; i > j;i--)
		{
			temp = p[i] - '0';
			DOUBLE1 = DOUBLE1*0.1+temp;
		}
		return (DOUBLE+0.1*DOUBLE1);
	}
	if( count2 == 2 && count1 == 1 && p[0] == '-' && p[1] != '.' )
	{
		for( i = 1;i < j;i++)
		{
			temp = p[i] - '0';
			DOUBLE = DOUBLE*10 + temp;
		}
		for( i = length-1; i > j;i--)
		{
			temp = p[i] - '0';
			DOUBLE1 = DOUBLE1*0.1+temp;
		}
		return -(DOUBLE+0.1*DOUBLE1);
	}
	return 0;
}
#endif