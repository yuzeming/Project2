#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/opencv.hpp"
#include <cmath>
#include <iostream>
#include <fstream>

#define sqr(x) ((x)*(x))

using namespace cv;
using namespace std;

const double PI = 3.1415926535897932384626433832795;

double step_len = 2 ;
double min_r = 4;
double max_r = 10;
double eps = 4;
Mat src, threshold_output, res;

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
double calcR(const vector<vector<Point> >& poly, const Point2d & x,const Point2d & y)
{
	double ret_R = max_r,tmpR=0;
	for( size_t i = 0; i < poly.size(); ++i )
	{
		Point2d c,d,cross;
		size_t len = poly[i].size();
		for (size_t j = 0; j<len; ++j)
			{
				c = poly[i][j];
				d = poly[i][(j+1)%len];
				Point2d tmp = d - c;
				if (calcCrossPoint(x,y,c,d,cross)&& norm(cross-x) > 0.1 )//
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


#define CalcCircle(a)  ((a).o) = ((a).x + (a).p * (a).r)
class TCircle
{
public:
	Point2d x, p; //垂点，法线
	double r; // 半径
	Point2d o; //圆心 o = x + p * r
	double s; //与之前的点不重合的面积
};

vector< vector<TCircle> > center; //圆心 <x,y,r>
#define VecLast(x) (x[x.size()-1])

bool check(Point2d o,double r)
{
	for (size_t i = 0; i < center.size(); ++i)
	{
		int size = center[i].size();
		for (size_t j = 0; j < size; ++j)
		{
			if (!(i + 1 == center.size() && j + 1 == size) && norm(center[i][j].o - o) < 3)
				return false;
			if (norm(center[i][j].o - o) < 1)
				return false;
		}
	}
	return true;
}

bool CheckLine(const vector<vector<Point> >& poly,Point2d x, Point2d y)
{
	Point2d a, b, ret;
	for (size_t i = 0; i < poly.size(); ++i)
	{
		size_t len = poly[i].size();
		for (size_t j = 0; j < len; ++j)
		{
			a = poly[i][j];
			b = poly[i][(j + 1) % len];
			if (calcCrossPoint(a, b, x, y, ret))
				return false;
		}
	}
	return true;
}

void addCircle(const vector<vector<Point> >& poly,const Point2d x, const Point2d p,double r)
{
	TCircle tmp;
	tmp.x = x;
	tmp.p = p;
	tmp.r = r;
	tmp.o = x + p * r;

	if (!check(tmp.o,tmp.r))
		return;

	if (center.size() == 0 || ( VecLast(center).size() != 0 && (!CheckLine(poly,VecLast(VecLast(center)).o, tmp.o) || norm(VecLast(VecLast(center)).o- tmp.o) > eps*2) ))
		center.push_back(vector<TCircle>());
	VecLast(center).push_back(tmp);
	
}

void CalcS()
{
	for (size_t i = 0; i < center.size(); ++i)
	{
		int size = center[i].size();
		for (size_t j = 0; j < size; ++j)
		{

			double rx = 0.1;
			double r_step = 0.2;
			center[i][j].s = sqr(rx) * PI;
			while (rx < center[i][j].r)
			{
				int k = int(rx * PI * 2 / r_step);
				double unit_s = (sqr(rx + r_step) - sqr(rx)) * PI / k;
				double tmp_r = rx + r_step / 2;
				for (int phi = 0; phi < k; ++phi)
				{
					Point2d p = center[i][j].o + Point2d(cos(2 * PI*phi / k), sin(2 * PI*phi / k)) * tmp_r;
					bool flag = true;
					for (size_t q = max(int(j) - 10, 0); q < j && flag; ++q)
						if (norm(p - center[i][q].o) < center[i][q].r)
							flag = false;
					if (flag)
						center[i][j].s += unit_s;
				}
				rx += r_step;
			}
			//printf("%lf %lf\n", center[i][j].r, center[i][j].s);
		}
	}
}

int main()
{
	vector<vector<Point> > contours;
	vector<Vec4i> hierarchy;
	
	src = imread("mew.png",IMREAD_GRAYSCALE);

	cvtColor(src,res,COLOR_GRAY2BGR);

	threshold(src,threshold_output,70,255,THRESH_BINARY_INV);
	
	findContours(threshold_output,contours,hierarchy,RETR_TREE,CHAIN_APPROX_SIMPLE);
	
	vector<vector<Point> > contours_poly( contours.size() );
	for( size_t i = 0; i < contours.size(); ++i )
		approxPolyDP( Mat(contours[i]), contours_poly[i], 2, true);

	for( size_t i = 0; i < contours_poly.size(); ++i )
	{

		Point2d a,b,c,c1,x,y;
		size_t len = contours_poly[i].size();
		double R;
		for (size_t j = 0; j < len; ++j)
		{
			
			a = contours_poly[i][j];
			b = contours_poly[i][(j+1)%len];

			//直线上逐一模拟
			c = b - a; 

			int step_max = max(int( norm(c) / step_len  ),1);
			c1 = Point2d(c.y, -c.x)/ norm(c); //单位法线方向向量
			for (int k=0;k<step_max;++k)
			{
				x = c * ( (k + 0.5) / step_max ) + a ; // 取样点
				y = x + c1 * 1000 ; // 无限长射线
				R = calcR(contours_poly,x,y);
				addCircle(contours_poly,x,c1,R);
			}
			
			c = contours_poly[i][(j+2)%len];
			double alpha = atan2(a - b);
			double beta = atan2(c - b);

			if (beta > alpha)
				alpha += 2 * PI;
			double turn = (alpha - beta) / 2.0;
			x = b;
			c1 = ComplexMul(b - c, Point2d(cos(turn), sin(turn)));
			c1 = c1 / norm(c1);
			y = x + c1 * 1000;
			R = calcR(contours_poly, x, y);

			addCircle(contours_poly, x, c1, R);

		}
	}

	//平滑半径
	int KX = 5;
	for (size_t i = 0; i < center.size(); ++i)
	{
		vector <double> R;
		size_t len = center[i].size();
		R.push_back(0.0);
		for (size_t j = 0; j < len; ++j)
			R.push_back(center[i][j].r + VecLast(R));
		for (int j = 0; j < len; ++j)
		{
			center[i][j].r = (R[j + 1] - R[max(j + 1 - KX, 0)]) / min(j + 1, KX);
			CalcCircle(center[i][j]);
		}
	}

	for (size_t i = 0; i < center.size(); ++i)
		for (size_t j = i + 1; j < center.size(); ++j)
			for (size_t k = 0; k < 4; ++k)
			{
				Point2d a, b;
				a = (k & 1) ? center[i][0].o : VecLast(center[i]).o;
				b = (k & 2) ? center[j][0].o : VecLast(center[j]).o;
				if (norm(a - b) < eps * 4 && CheckLine(contours_poly, a, b))  //
				{
					if (k & 1)	reverse(center[i].begin(), center[i].end());
					if (~k & 2) reverse(center[j].begin(), center[j].end());
					for (size_t p = 0; p < center[j].size(); ++p)
						center[i].push_back(center[j][p]);
					center.erase(center.begin() + j);
					--j;
					break;
				}
			}
	
	for (int i = 0; i < center.size(); ++i)
	{
		double len = 0;
		for (int j = 1; j < center[i].size(); ++j)
			len += norm(center[i][j].o - center[i][j - 1].o);
		if (len < 10)
		{
			center.erase(center.begin() + i);
			--i;
		}
	}
	
	//排序
	for (int i = 1; i < center.size(); ++i)
	{
		bool reverse_flag = false;
		double tmp = 1e99;
		Point2d last = VecLast(center[i - 1]).o;
		int next = -1;

		for (int j = i; j < center.size(); ++j)
		{
			if (norm(last - center[j][0].o) < tmp)
			{
				tmp = norm(last - center[j][0].o);
				reverse_flag = false;
				next = j;
			}
			if (norm(last - VecLast(center[j]).o) < tmp)
			{
				tmp = norm(last - VecLast(center[j]).o);
				reverse_flag = true;
				next = j;
			}
		}

		if (next != -1)
		{
			if (reverse_flag)
				reverse(center[next].begin(), center[next].end());
			swap(center[i], center[next]);
		}
	}


	CalcS();
	
	double X_MAX = 200;
	double Y_MAX = 300;
	double resize_k = min(X_MAX / src.cols,   src.rows/ Y_MAX);
	
	int power_on= 10;
	int power_off = -10;
	int power_stop = -100;
	double expextrude = 2500;

	ofstream gcode("test.gcode");


	for (size_t i = 0; i < center.size(); ++i)
	if(center[i].size() >= 5)
	{
		size_t len = center[i].size();
		size_t endp = len - 1;
		double end_len = 0;
		while (end_len < 5 && endp > 0)
		{
			end_len += resize_k*norm(center[i][endp - 1].o - center[i][endp].o);
			endp--;
		}
		double last_r = -1;
		//移动到起始点
		for (size_t j = 0; j < len; ++j)
		{
			gcode << "G1 X" << center[i][j].o.x*resize_k << " Y" << center[i][j].o.y*resize_k << " F"<< int(expextrude / pow( center[i][j].s * sqr(resize_k),1.2) )<< endl;
			if (j == 0)
			{
				double power = power_on;
				gcode << "G90" << endl << "M400" << endl;
				if (power_on < 0)
					gcode << "M804" << endl; // 回抽
				else
					gcode << "M803" << endl; // 挤出
				gcode << "M801 S" << abs(power_on) << endl; //设定气压

				gcode << "G4 P300" << endl; //等待300ms

				last_r = center[i][j].r;
			}
			if (j == endp)
			{
				gcode << "G90" << endl << "M400" << endl;
				gcode << "M804" << endl;
				gcode << "M801 S"<< abs(power_off) << endl; // 回抽
			}
		}
		gcode << "G4 P300" << endl;
	}
	gcode << "M804" << endl;
	gcode << "M801 S" << abs(power_stop) << endl; // 回抽
	gcode << "G4 P400" << endl;
	gcode << "G28" << endl; //回原点

	//debug
	for (int i = 0; i < center.size(); ++i)
	{
		size_t len = center[i].size();
		for (int j = 1; j < len; ++j)
		{
			line(res, center[i][j].o, center[i][j - 1].o, Scalar(255, 255, 0), 2);
			//circle(res, center[i][j].o, center[i][j].r+0.9, Scalar(0, 0, 255), 1);
		}
		if (i + 1 < center.size())
			line(res, center[i + 1][0].o, center[i][len - 1].o, Scalar(0, 0, 255),1);
		circle(res, center[i][0].o, 2, Scalar(0, 0, 255), 2);
		circle(res, VecLast(center[i]).o, 2, Scalar(0, 0, 255), 2);
		//imshow("test", res);
		//waitKey(0);
	}

	imshow("test",res);
	waitKey(0);
	 
	return 0;
}
