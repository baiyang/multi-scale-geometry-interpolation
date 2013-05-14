#include "MsInterpolation.h"
#include "LeafNodeInterpolation.h"
#include <cmath>

#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <map>
using namespace std;


void get_beta_vector( MyVector3f &in_vect, double cos_theta, MyVector3f &normal, MyVector3f &out_vect, int sign )
{
	in_vect.normalize();
	normal.normalize();

	double sin_theta = sqrt( 1 - cos_theta * cos_theta) * sign;

	out_vect = cos_theta * in_vect + sin_theta * ( normal.cross( in_vect));
	out_vect.normalize();
}


void leaf_node_interpolation( double t, FaceNode *leaf_node, MyMesh &src_mesh, MyMesh &target_mesh )
{
	queue<MyFaceHandle> closed_list; // 保存已经处理过的triangles
	map<int, bool> face_done;        // 标记每个triangle的状态
	map<int, int> r_idx;             // 反向索引在pts_index里的值


	for (int i = 0; i != leaf_node->idx.size(); i++)
	{
		face_done.insert( pair<int, bool>(leaf_node->idx[i], false) );
	}

	/*** 分配空间  ***/
	leaf_node->pts.resize( leaf_node->pts_index.size() );


	for (int i = 0; i != leaf_node->pts_index.size(); i++)
	{
		r_idx.insert( pair<int, int>(leaf_node->pts_index[i], i) );
	}

	/*** 初始化closed_list ***/
	MyFaceHandle curr_face_handle = MyFaceHandle( leaf_node->idx[0] );
	vector<double> src_edge_len;// 表示三角形三边的长度
	vector<double> target_edge_len;
	vector<int> curr_vertex_idx;
	double a0, b0, c0, a1, b1, c1, a, b, c;
	int x, y, z;
	MyVector3f p_x, p_y, p_z;
	MyVector3f in_vect, out_vect, face_normal;
	MyMesh::FaceEdgeIter fe_iter;


	for (fe_iter = src_mesh.fe_begin( curr_face_handle); fe_iter; ++fe_iter)
	{
		src_edge_len.push_back( src_mesh.calc_edge_length( fe_iter.handle() ) );
		target_edge_len.push_back( target_mesh.calc_edge_length( fe_iter.handle() ));

		curr_vertex_idx.push_back( src_mesh.from_vertex_handle( fe_iter.current_halfedge_handle() ).idx() );
	}

	a0 = src_edge_len[0];
	b0 = src_edge_len[1];
	c0 = src_edge_len[2];

	a1 = target_edge_len[0];
	b1 = target_edge_len[1];
	c1 = target_edge_len[2];

	x = curr_vertex_idx[0];
	y = curr_vertex_idx[1];
	z = curr_vertex_idx[2];

	/*** 插值后的边长 ***/
	a = ( 1 - t) * a0 + t * a1;
	b = ( 1 - t) * b0 + t * b1;
	c = ( 1 - t) * c0 + t * c1;

	MyMesh::Point tmp_x = src_mesh.point( MyVertexHandle( x ) );
	MyMesh::Point tmp_y = src_mesh.point( MyVertexHandle( y ) );
	p_x = MyVector3f( tmp_x[0], tmp_x[1], tmp_x[2] );
	p_y = MyVector3f( tmp_y[0], tmp_y[1], tmp_y[2] );

	p_y = (p_y - p_x) / (p_y - p_x).norm() * a + p_x;
	double cos_x = ( a*a + c*c - b*b) / (2*a*c);
	in_vect = p_y - p_x;
	
	/*** 计算当前face的normal ***/
	MyMesh::Normal tmp_n = src_mesh.calc_face_normal( curr_face_handle );
	face_normal = MyVector3f(tmp_n[0], tmp_n[1], tmp_n[2]);

	get_beta_vector(in_vect, cos_x, face_normal, out_vect);

	p_z = c * out_vect + p_x;

	leaf_node->pts[ r_idx[x] ] = p_x;
	leaf_node->pts[ r_idx[y] ] = p_y;
	leaf_node->pts[ r_idx[z] ] = p_z;

	closed_list.push( curr_face_handle );
	face_done[curr_face_handle.idx() ] = true;

	/***初始化完毕***/

    while( !closed_list.empty() )
	{
		MyFaceHandle done_face = closed_list.front();
		closed_list.pop();

		MyMesh::FaceFaceIter ff_iter;
		for (ff_iter = src_mesh.ff_begin(done_face); ff_iter; ++ff_iter)
		{
			/*** overlap_halfedge 半边是属于done_face里面的 ***/
			MyHalfedgeHandle overlap_halfedge = ff_iter.current_halfedge_handle();
			curr_face_handle = ff_iter.handle();

			/***不存在，或则该face已经插值完都返回 ***/
			if( !face_done.count(curr_face_handle.idx())  || face_done[curr_face_handle.idx()])
			{
				continue;
			}

			double src_dihedral_angle = src_mesh.calc_dihedral_angle( overlap_halfedge );
			double target_dihedral_angle = target_mesh.calc_dihedral_angle( overlap_halfedge );

			double dihedral_angle = (1 - t) * src_dihedral_angle + t * target_dihedral_angle;

			x = src_mesh.from_vertex_handle( overlap_halfedge ).idx();
			y = src_mesh.from_vertex_handle( src_mesh.next_halfedge_handle( overlap_halfedge )).idx();
			z = src_mesh.from_vertex_handle( src_mesh.next_halfedge_handle(src_mesh.next_halfedge_handle( overlap_halfedge )) ).idx();

			int sign = 1;
			if( dihedral_angle < 0)
			{
				sign = -1;
			}

			/***注意：获取已插值好的三角形的face normal（即这里的done_face),保存在in_vect***/
			get_triangle_normal( leaf_node->pts[ r_idx[x]], leaf_node->pts[ r_idx[y] ], leaf_node->pts[r_idx[z] ], in_vect);
			
			/***旋转轴***/
			face_normal =  leaf_node->pts[ r_idx[y] ] - leaf_node->pts[ r_idx[x] ];

			/***out_vect保存了新插值后的三角形的法向量***/
			get_beta_vector(in_vect, cos( dihedral_angle ), face_normal, out_vect, sign);

			MyHalfedgeHandle curr_halfedge = src_mesh.opposite_halfedge_handle( overlap_halfedge);

			c = (leaf_node->pts[ r_idx[y] ] - leaf_node->pts[ r_idx[x] ]).norm();
			a0 = src_mesh.calc_edge_length( src_mesh.next_halfedge_handle( curr_halfedge ));
			b0 = src_mesh.calc_edge_length( src_mesh.next_halfedge_handle( src_mesh.next_halfedge_handle( curr_halfedge )));

			a1 = target_mesh.calc_edge_length( target_mesh.next_halfedge_handle( curr_halfedge ) );
			b1 = target_mesh.calc_edge_length( target_mesh.next_halfedge_handle( target_mesh.next_halfedge_handle( curr_halfedge )));

			a = (1 - t) * a0 + t * a1;
			b = (1 - t) * b0 + t * b1;
			double cos_z = (b*b + c*c - a*a) / (2*b*c);

			face_normal *= -1;
			get_beta_vector( face_normal, cos_z, out_vect, out_vect);

			z = src_mesh.to_vertex_handle( src_mesh.next_halfedge_handle( curr_halfedge )).idx();
			p_z = out_vect * b + leaf_node->pts[ r_idx[y] ];

			leaf_node->pts[ r_idx[z] ] = p_z;

			/*** 插值完成 ***/
			face_done[ curr_face_handle.idx()]  = true;
			closed_list.push( curr_face_handle );
		}
	}
}

void get_triangle_normal( MyVector3f &x, MyVector3f &y, MyVector3f &z, MyVector3f &out_normal )
{
	MyVector3f n0 = z - x;
	MyVector3f n1 = y - x;

	out_normal = n1.cross( n0 );

	out_normal.normalize();
}
