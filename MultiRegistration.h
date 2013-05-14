#pragma once


#include <iostream>
#include <vector>
#include <Eigen/Dense>

#include "Pair.h"

typedef Eigen::Matrix<double, 3, 3>  MyMatrix3f;
typedef Eigen::Vector3d MyVector3f;
typedef Eigen::RowVector3d MyRowVector3f;

typedef Eigen::MatrixXd MyMatrixXf;

const int MAX_ITERATION_R = 50;

class Pair;


class PairVertex
{
public:
	MyVector3f a;
	MyVector3f b;
};


typedef std::vector<PairVertex> PointList;


class MultiRegistration
{
public:
	MultiRegistration(int m, int p, std::vector<Pair> *pair, std::vector<PointList> *pointlist);
	~MultiRegistration(void);

	void init();

	void block_index(int &b_x, int &b_y, int &g_x, int &g_y, int block_size_x=3, int block_size_y=3);

	void test();

	void get_R_and_T(std::vector<MyMatrix3f> &R, std::vector<MyVector3f> &T);

private:
	void pinv( MyMatrixXf &input, MyMatrixXf &pinv_matrix);


private:
	int M; // M¸öviews
	int P; // P¸öpairs

	MyMatrixXf R;
	MyMatrixXf T;

	MyMatrixXf Ca;
	MyMatrixXf Cb;

	MyMatrixXf H_xx;
	MyMatrixXf H_xy;
	MyMatrixXf H_yx;
	MyMatrixXf H_yy;
	MyMatrixXf Q_r;
	MyMatrixXf N;
	MyMatrixXf gamma;
	MyMatrixXf Q_xx;
	MyMatrixXf Q_xy;
	MyMatrixXf Q_yx;
	MyMatrixXf Q_yy;
	MyMatrixXf Q_rt;
	MyMatrixXf Q;

	std::vector<PointList> *pl;
	std::vector<Pair> *pair_index;
};

