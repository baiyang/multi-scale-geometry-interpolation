#include "MultiRegistration.h"
#include <Eigen/SVD>

MultiRegistration::MultiRegistration( int m, int p, std::vector<Pair> *pair, std::vector<PointList> *pointlist )
{
	M = m;
	P = p;

	pair_index = pair;
	pl = pointlist;

}


MultiRegistration::~MultiRegistration(void)
{
}

void MultiRegistration::init()
{
	R.resize(3, 3 * M);
	T.resize(1, 3 * M);

	Ca.resize(3 * M, 3 * P);
	Cb.resize(3 * M, 3 * P);

	H_xx.resize( 3 * P, 3 * P);
	H_xy.resize( 3 * P, 3 * P);
	H_yx.resize( 3 * P, 3 * P);
	H_yy.resize( 3 * P, 3 * P);
	Q_r.resize( 3 * M, 3 * M);
	N.resize( 3 * P, 3 * P);

	int b_x, b_y; // block的index
	int g_x, g_y; // 原矩阵的index
	int a, b;


	/*** begin:构建Ca和Cb矩阵 ***/
	for(int i = 0; i != M; i++)
	{
		for(int j = 0; j != P; j++)
		{
			a = (*pair_index)[j].a;
			b = (*pair_index)[j].b;

			b_x = i;
			b_y = j;
			block_index(b_x, b_y, g_x, g_y);

			if( i == a )
			{
				Ca.block<3, 3>(g_x, g_y) = MyMatrix3f::Identity();
			}else
			{
				Ca.block<3, 3>(g_x, g_y) = MyMatrix3f::Zero();
			}

			if( i == b )
			{
				Cb.block<3, 3>(g_x, g_y) = MyMatrix3f::Identity();
			}else
			{
				Cb.block<3, 3>(g_x, g_y) = MyMatrix3f::Zero();
			}
		}
	}

	/*** end ***/

	/*** begin:构建H_xx,N等矩阵***/
	MyVector3f x, y;
	MyRowVector3f x_r, y_r;
	int Nu;
	MyMatrix3f xx, xy, yx, yy;

	H_xx.setZero();
	H_yy.setZero();
	H_xy.setZero();
	H_yx.setZero();
	N.setZero();

	for (int u = 0; u != P; u++)
	{

		Nu = (*pl)[u].size();

		xx = MyMatrix3f::Zero();
		xy = MyMatrix3f::Zero();
		yx = MyMatrix3f::Zero();
		yy = MyMatrix3f::Zero();

		for (int i = 0; i != Nu; i++)
		{
			x = (*pl)[u][i].a;
			y = (*pl)[u][i].b;
			x_r = x.transpose();
			y_r = y.transpose();

			xx += x * x_r;
			xy += x * y_r;
			yx += y * x_r;
			yy += y * y_r;
		}

		b_x = u;
		b_y = u;
		block_index(b_x, b_y, g_x, g_y);

		N.block<3, 3>(g_x, g_y) = Nu * MyMatrix3f::Identity();

		H_xx.block<3, 3>(g_x, g_y) = xx;
		H_xy.block<3, 3>(g_x, g_y) = xy;
		H_yx.block<3, 3>(g_x, g_y) = yx;
		H_yy.block<3, 3>(g_x, g_y) = yy;
	}

	/*** end ***/

	/*** begin: 计算Q_r  ***/
	Q_r = Ca * H_xx * Ca.transpose() + Cb * H_yy * Cb.transpose() - Ca * H_xy * Cb.transpose() - Cb * H_yx * Ca.transpose();

	MyMatrixXf Ca_Cb = Ca - Cb;
	MyMatrixXf Ca_Cb_t = Ca_Cb.transpose();

	MyMatrixXf tmp_pinv = Ca_Cb * N * Ca_Cb_t;
	pinv(tmp_pinv, tmp_pinv);

	MyMatrixXf ZK = N * Ca_Cb_t * tmp_pinv * Ca_Cb * N;

	gamma.resize( P, P);
	
	for (int i = 0; i != P; i++)
	{
		for (int j = 0; j != P; j++)
		{
			b_x = i;
			b_y = j;
			block_index(b_x, b_y, g_x, g_y);

			gamma(i, j) = ZK(g_x, g_y);
		}
	}

	std::vector<MyVector3f> x_mean;
	std::vector<MyVector3f> y_mean;

	MyVector3f x_sum, y_sum;
	for (int u = 0; u != P; u++)
	{
		x_sum = y_sum = MyVector3f::Zero();

		Nu = (*pl)[u].size();
		for (int i = 0; i != Nu; i++)
		{
			x_sum += (*pl)[u][i].a;
			y_sum += (*pl)[u][i].b;
		}

		x_mean.push_back( x_sum / Nu);
		y_mean.push_back( y_sum / Nu);

	}


	MyMatrixXf G_xx;
	MyMatrixXf G_xy;
	MyMatrixXf G_yx;
	MyMatrixXf G_yy;

	G_xx.resize( 3 * P, 3 * P);
	G_xy.resize( 3 * P, 3 * P);
	G_yx.resize( 3 * P, 3 * P);
	G_yy.resize( 3 * P, 3 * P);
	
	for (int u = 0; u != P; u++)
	{
		for (int n = 0; n != P; n++)
		{
			b_x = u;
			b_y = n;

			block_index(b_x, b_y, g_x, g_y);
			xx = gamma(u, n) * x_mean[u] * x_mean[n].transpose();
			xy = gamma(u, n) * x_mean[u] * y_mean[n].transpose();
			yx = gamma(u, n) * y_mean[u] * x_mean[n].transpose();
			yy = gamma(u, n) * y_mean[u] * y_mean[n].transpose();

			G_xx.block<3, 3>(g_x, g_y) = xx;
			G_xy.block<3, 3>(g_x, g_y) = xy;
			G_yx.block<3, 3>(g_x, g_y) = yx;
			G_yy.block<3, 3>(g_x, g_y) = yy;

		}
	}

	Q_xx = Ca * G_xx * Ca.transpose();
	Q_xy = Ca * G_xy * Cb.transpose();
	Q_yx = Cb * G_yx * Ca.transpose();
	Q_yy = Cb * G_yy * Cb.transpose();

	Q_rt = Q_xy + Q_yx - Q_yy - Q_xx;

	Q = Q_r + Q_rt;

	MyMatrixXf Sj;
	MyMatrixXf V, U;
	MyMatrixXf last_R;

	last_R.resize( 3, M * 3);

	last_R.setZero();

	int count = MAX_ITERATION_R;
	while( count-- )
	{
		for (int j= 0; j != M; j++)
		{
			Sj = Q.block(j * 3, 0, 3, j * 3) * R.block(0, 0, 3, j * 3).transpose() + Q.block( 3 * j , j * 3 + 3, 3, (M - j - 1) * 3 ) * R.block(0, j*3 + 3, 3, (M-j-1)*3).transpose();
			Sj = -1 * Sj;

			Eigen::JacobiSVD<MyMatrixXf> svd( Sj, Eigen::ComputeThinU | Eigen::ComputeThinV);
			V = svd.matrixV();
			U = svd.matrixU();

			if( V.determinant() < 0 )
			{
				for (int i = 0; i != 3; i++)
				{
					V(i, 2) *= -1;
				}
			}


			if( U.determinant() < 0 )
			{
				for (int i = 0; i != 3; i++)
				{
					U(i, 2) *= -1;
				}
			}

			R.block<3, 3>(0, j * 3) = V * U.transpose();

		}

		MyMatrixXf diff = R - last_R;
		if( fabs( diff.maxCoeff() ) < 1 * 1e-2 && fabs( diff.minCoeff() ) < 1 * 1e-2 )
		{
			break;
		}
		
		last_R = R;
	}

	std::cout<<count<<std::endl;

	MyMatrixXf Z;
	MyMatrixXf RCa;
	MyMatrixXf RCb;
	MyMatrixXf RCa_x;
	MyMatrixXf RCb_y;

	RCa = R * Ca;
	RCb = R * Cb;

	RCa_x.resize( 3 * P, 1);
	RCb_y.resize( 3*P, 1);

	for (int i = 0; i != P; i++)
	{
		b_x = i;
		b_y = 0;
		block_index(b_x, b_y, g_x, g_y);

		RCa_x.block<3, 1>(g_x, 0) = RCa.block<3, 3>(0, g_x) * x_mean[i];
		RCb_y.block<3, 1>(g_x, 0) = RCb.block<3, 3>(0, g_x) * y_mean[i];

	}

	Z = RCa_x - RCb_y;

	T = (-1) * tmp_pinv * Ca_Cb * N * Z;
}

void MultiRegistration::block_index( int &b_x, int &b_y, int &g_x, int &g_y, int block_size_x/*=3*/, int block_size_y/*=3*/ )
{
	g_x = b_x * block_size_x;
	g_y = b_y * block_size_y;
}

void MultiRegistration::pinv( MyMatrixXf &input, MyMatrixXf &pinv_matrix )
{
	double epsilon = 1e-5;
	int x, y;

	x = input.rows();
	y = input.cols();

	Eigen::JacobiSVD<MyMatrixXf> svd(input, Eigen::ComputeThinU | Eigen::ComputeThinV);

	Eigen::Matrix<double, Eigen::Dynamic, 1> inv_singular(y, 1);

	for(int i = 0; i < y; i++)
	{
		if ( fabs( svd.singularValues()[i] ) < epsilon )
		{
			inv_singular[i] = 0;
		}else
		{
			inv_singular[i] = 1 / svd.singularValues()[i]; 
		}
	}

	pinv_matrix = svd.matrixV() * inv_singular.asDiagonal() * svd.matrixU().transpose();
}

void MultiRegistration::test()
{
	MyMatrixXf m = MyMatrixXf::Random(3, 2);

	MyMatrixXf output;

	std::cout<<m.block(0, 0, 1,0)<<std::endl;

	pinv(m , output);
}

void MultiRegistration::get_R_and_T( std::vector<MyMatrix3f> &_R, std::vector<MyVector3f> &_T )
{
	int b_x, b_y; // block的index
	int g_x, g_y; // 原矩阵的index

	for (int i = 0; i != M; i++)
	{
		b_x = 0;
		b_y = i;

		block_index(b_x, b_y, g_x, g_y);

		_R.push_back( R.block<3, 3>(g_x, g_y) );

		_T.push_back( T.block<3, 1>(g_y, 0) );

	}
}
