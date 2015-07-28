#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/opencv.hpp"
#include <cmath>
#include <iostream>

using namespace cv;
using namespace std;

const double PI = 3.1415926535897932384626433832795;

double step_len = 3 ;
double min_r = 5;
double max_r = 20;
Mat res_tmp;



double CrossDot(const Point2d &a,const Point2d &b,const Point2d &c)
{
	return (b-a).cross(c-a);
}

bool calcCrossPoint(const Point2d & a,const Point2d & b,const Point2d & c,const Point2d & d,Point2d &ret)
{
	double s1 = CrossDot(c,d,a),s2 = CrossDot(c,d,b);
	double s3 = CrossDot(a,b,c),s4 = CrossDot(a,b,d);
	if (s1*s2 > 0 || s3 * s4 >0)
		return false;

	ret = a + (b - a) * (s1 / (s1-s2)); 

	return true;
}

double calcAngle(const Point2d & a,const Point2d & b,const Point2d & c,const Point2d & d)
{
	Point2d x = b - a, y = d - c;
	return acos(x.dot(y) / (norm(x) * norm(y)));
}

double CalcR2(const Point2d & x,const Point2d & y,const Point2d &z)
{
	Point2d a = y-x, b = z-x;
	double cos_alpha = abs(a.dot(b) / (norm(a) * norm(b)));
	return norm(b) / (2.0 * cos_alpha);
}

//在多边形轮廓中，计算圆心在x->y射线上，且x点在圆上的圆，最大直径。
//同时返回与该圆恰好相交的线段。
double calcR(const vector<vector<Point> >& poly, const Point2d & x,const Point2d & y,int a,int b,int &ret_x,int& ret_y)
{
	double ret_R = max_r,tmpR=0;
	for( size_t i = 0; i < poly.size(); ++i )
	{
		Point2d c,d,cross;
		size_t len = poly[i].size();
		for (size_t j = 0; j<len; ++j)
			if (i!=a||(j!=b))
			{
				c = poly[i][j];
				d = poly[i][(j+1)%len];
				Point2d tmp = d - c;
				if (calcCrossPoint(x,y,c,d,cross))
				{
					double alpha = 	calcAngle(x,y,c,d);
					tmpR = norm(cross-x) / (1 + 1.0/sin(alpha));
					if (tmpR < ret_R)
					{
						ret_x = i;
						ret_y = j;
						ret_R = tmpR;
					}
				}
				tmpR = min( CalcR2(x,y,c),CalcR2(x,y,d));
				if (tmpR < ret_R)
				{
					ret_x = i;
					ret_y = j;
					ret_R = tmpR;
				}
			}
	}
	return ret_R;
}

Point2d ComplexMul(const Point2d& a,const Point2d& b)
{
	return Point2d(a.x*b.x-a.y*b.y,a.y*b.x+a.x*b.y); 
}

map< pair<int,int>, vector<Point2d> > crosspoint;
void addCirclePoint(const vector<vector<Point> >& poly,int a,int b,const Point2d& x)
{
	if (CrossDot(poly[a][b],poly[a][(b+1)%poly[a].size()],x) < 1e-3)
		crosspoint[make_pair(a,b)].push_back(x);
}

bool check(int a,int b,const Point2d& x,double eps)
{
	auto it = crosspoint.find(make_pair(a,b));
	if (it!=crosspoint.end())
	{
		size_t len = it->second.size();
		for (size_t i=0;i<len;++i)
			if (norm(it->second.at(i)-x)<eps)
				return false;
	}
	return true;
}

double atan2(const Point2d a)
{
	return atan2(a.y,a.x);
}

vector<pair<Point2d,double> > center; //圆心 <x,y,r>
int addCircle(const Point2d o,double R)
{
	center.push_back(make_pair(o,R));
	return center.size()-1;
}

map<int,vector<int> > lines;

void addLine(int c1,int c2)
{
	lines[c1].push_back(c2);
	lines[c2].push_back(c1);
}

int main()
{
	Mat src,threshold_output,res;
	vector<vector<Point> > contours;
	vector<Vec4i> hierarchy;
	

	src = imread("thu.png",IMREAD_GRAYSCALE);
	erode(src,src,Mat());

	cvtColor(src,res,COLOR_GRAY2BGR);
	//res.zeros(res.size());

	threshold(src,threshold_output,70,255,THRESH_BINARY_INV);
	
	findContours(threshold_output,contours,hierarchy,RETR_TREE,CHAIN_APPROX_SIMPLE);
	
	vector<vector<Point> > contours_poly( contours.size() );
	for( size_t i = 0; i < contours.size(); ++i )
		approxPolyDP( Mat(contours[i]), contours_poly[i], 1, true);

	for( size_t i = 0; i < contours_poly.size(); ++i )
	{
		size_t len = contours_poly[i].size();
		cout << "Poly ("<< len  << ") : ";
		for (size_t j = 0; j < len; ++j)
		{
			//line(res,contours_poly[i][j],contours_poly[i][(j+1)%len],Scalar(0,255,0));
			cout<<contours_poly[i][j]<<" ";
		}
		cout<<endl;
	}
	
	int p,q;
	for( size_t i = 0; i < contours_poly.size(); ++i )
	if ( contours_poly[i].size() > 3)
	{

		Point2d a,b,c,c1,x,y;
		size_t len = contours_poly[i].size();
		int lastCid = -1;
		double lastR = max_r;
		for (size_t j = 0; j < len; ++j)
		{
			
			res.copyTo(res_tmp);
			a = contours_poly[i][j];
			b = contours_poly[i][(j+1)%len];

			//直线上逐一模拟
			c = b - a; 
			if (norm(c) >= step_len)
			{
				int step_max = int( norm(c) / step_len);
				c1 = c / norm(c);
				for (int k=0;k<step_max;++k)
				{
					x = c * ( (k + 0.5) / step_max ) + a ;
					if (!check(i,j,x,step_len))
					{
						lastCid = -1;
						continue;
					}
					y = x + Point2d(c1.y,-c1.x) * 1000 ;
					double R = calcR(contours_poly,x,y,i,j,p,q);
					R = min(R,lastR*1.3);
					lastR = R;
					addCirclePoint(contours_poly,p,q,x+Point2d(c1.y,-c1.x) * (2 * R));
					Point2d o = x+Point2d(c1.y,-c1.x) * R;

					int cid = addCircle(o,R);

					if (lastCid != -1)
						addLine(lastCid,cid);
					lastCid = cid;

					//circle(res,o,R,Scalar(255,255,0),1);

					//line(res,x,y,Scalar(0,0,255));
					//imshow("test",res);
					//waitKey(0);
				}
			}

			//线段终点拐角，逐一模拟
/*			c = contours_poly[i][(j+2)%len];
			double alpha = atan2(a-b);
			double beta = atan2(c-b);
			double phi = 
			x = b;
			y = x + ComplexMul((a-b),Point2d(cos(phi),sin(phi)))*1000;
			line(res,x,y,Scalar(0,0,255));
	*/		
			

			//imshow("test",res);
			//waitKey(0);
		}

	}


//	for (size_t i=0;i<center.size();++i)
//		if (lines[i].size()==1)
//			for (size_t j=0;j<center.size();++j)
//				if (i!=j && norm(center[i].first-center[j].first) < center[i].second+center[j].second+step_len )
//					addLine(i,j);

	for (size_t i=0;i<lines.size();++i)
	{
		for (size_t j=0;j<lines[i].size();++j)
			if (lines[i][j]>i)
				line(res,center[i].first,center[lines[i][j]].first,Scalar(255,255,0),1);
		//circle(res,center[i].first,center[i].second,Scalar(255,0,0),1);
	}
	
	imshow("test",res);
	waitKey(0);

	return 0;
}