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
double calcR(const vector<vector<Point> >& poly, const Point2d & x,const Point2d & y,int a,int b)
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
					ret_R =min(ret_R, norm(cross-x) / (1 + 1.0/sin(alpha)));
				}
				ret_R = min(ret_R, min(CalcR2(x, y, c), CalcR2(x, y, d)));
			}
	}
	return ret_R;
}

Point2d ComplexMul(const Point2d& a,const Point2d& b)
{
	return Point2d(a.x*b.x-a.y*b.y,a.y*b.x+a.x*b.y); 
}

double atan2(const Point2d a)
{	return atan2(a.y,a.x);	}

vector<pair<Point2d,double> > center; //圆心 <x,y,r>

int addCircle(const Point2d o,double R)
{
	center.push_back(make_pair(o,R));
	return center.size()-1;
}

double check(const Point2d& x)
{
	double ret = max_r;
	for (size_t i = 0; i < center.size(); ++i)
		ret = min(ret, norm(center[i].first - x));
	return ret;
}

map<int,vector<int> > lines;
map<int, int> color;
void addLine(int c1,int c2)
{
	lines[c1].push_back(c2);
	lines[c2].push_back(c1);
}

#define sqr(x) ((x)*(x))

double fillColor(int x, int c)
{
	double ret = sqr(center[x].second);
	color[x] = c;
	for (size_t i = 0; i < lines[x].size(); ++i)
		if (color.find(lines[x][i]) == color.end())
			ret += fillColor(lines[x][i], c);
	return ret;
}

int main()
{
	Mat src,threshold_output,res;
	vector<vector<Point> > contours;
	vector<Vec4i> hierarchy;
	
	src = imread("mew.png",IMREAD_GRAYSCALE);

	cvtColor(src,res,COLOR_GRAY2BGR);

	threshold(src,threshold_output,70,255,THRESH_BINARY_INV);
	
	findContours(threshold_output,contours,hierarchy,RETR_TREE,CHAIN_APPROX_SIMPLE);
	
	vector<vector<Point> > contours_poly( contours.size() );
	for( size_t i = 0; i < contours.size(); ++i )
		approxPolyDP( Mat(contours[i]), contours_poly[i], 1, true);

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
					double R = check(x);
					if ( R <= step_len )
						continue;
					y = x + Point2d(c1.y,-c1.x) * 1000 ;
					R = calcR(contours_poly,x,y,i,j);
					R = min(R, lastR * 1.1);

					Point2d o = x+Point2d(c1.y,-c1.x) * R;

					int cid = addCircle(o,R);

					if (lastCid != -1 && norm(center[lastCid].first - o) <= step_len * 3)
						addLine(lastCid, cid);
					lastR = R;
					lastCid = cid;
				}
			}

			//线段终点拐角逐一模拟
/*			c = contours_poly[i][(j+2)%len];
			double alpha = atan2(a-b);
			double beta = atan2(b-c);
			double phi = 
			x = b;
			y = x + ComplexMul((a-b),Point2d(cos(phi),sin(phi)))*1000;
			line(res,x,y,Scalar(0,0,255));
	*/		
			
		}
	}


	int tot_color = 1;
	for (size_t i = 0; i < center.size(); ++i)
		if (color.find(i) == color.end())
			if (fillColor(i, tot_color) > 5)
				++tot_color;
			else
				fillColor(i, -1);

	for (size_t i = 0; i<center.size(); ++i)
		if (lines[i].size() == 1 && color[i] != -1)
		{
			circle(res, center[i].first, 1, Scalar(0,255, 0), 2);
			int cp = -1;
			double tmp = max_r + 1;
			for (size_t j = 0; j < center.size(); ++j)
				if ( i != j && norm(center[i].first - center[j].first) < tmp && ( color[i] != color[j] && color[j] != -1 || color[i] == color[j] && lines[j].size() == 1 ) )
				{
					tmp = norm(center[i].first - center[j].first);
					cp = j;
				}
			if (cp != -1 && tmp <= step_len * 3)
				addLine(i, cp);
		}


	
	for (size_t i = 0; i<center.size(); ++i)
	{
		for (size_t j = 0; j<lines[i].size(); ++j)
			line(res, center[i].first, center[lines[i][j]].first, Scalar(0, 255, 0 ), 2);
		//circle(res, center[i].first, 2, Scalar(255, 0, 0), 2);
	}





	imshow("test",res);
	waitKey(0);

	return 0;
}