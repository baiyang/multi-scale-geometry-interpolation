#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Tools/Utils/getopt.h>
#include <iostream>
#include <string>
#include <ctime>
#include <cstdlib>
#include <map>
#include <algorithm>
#include <set>
#include <fstream>
#include <assert.h>
#include "MsInterpolation.h"
#include "LeafNodeInterpolation.h"
#include "MultiRegistration.h"


MsInterpolation::MsInterpolation(void)
{
	face_root = NULL;

	min_size_node = 0x7fffffff;
}


MsInterpolation::~MsInterpolation(void)
{
	release_facenode_help(face_root);
}


void MsInterpolation::read_mesh_data( std::string filename, std::string target_file)
{
	bool ret;
	
	ret = IO::read_mesh( this->m_mesh, filename);
	if( ret == false )
	{
		printf("Reading Mesh Data Error in %s", __FUNCTION__);
		exit( -1 );
	}

	std::cout<<"Source Mesh Data:\n";
	std::cout << "# Vertices: " << m_mesh.n_vertices() << std::endl;
	std::cout << "# Edges   : " << m_mesh.n_edges() << std::endl;
	std::cout << "# Faces   : " << m_mesh.n_faces() << std::endl;

	ret = IO::read_mesh( target_mesh, target_file );
	if( ret == false )
	{
		printf("Reading Mesh Data Error in %s", __FUNCTION__);
		exit( -1 );
	}

	std::cout<<"Target Mesh Data:\n";
	std::cout << "# Vertices: " << m_mesh.n_vertices() << std::endl;
	std::cout << "# Edges   : " << m_mesh.n_edges() << std::endl;
	std::cout << "# Faces   : " << m_mesh.n_faces() << std::endl;
}


void MsInterpolation::write_mesh_data( std::string output )
{
	bool ret;
	IO::Options wopt;

	ret = OpenMesh::IO::write_mesh( this->m_mesh, output, wopt);

	if( ret == false )
	{
		printf("Writing Mesh Data Error In %s\n", __FUNCTION__);
	}
	else
	{
		printf("Writing Mesh Data Success!\n");
	}

}



void MsInterpolation::test()
{
	while( 1 )
	{
		release_facenode_help(face_root);
		build_hierarchy_on_face();

		printf("最小的节点的点数是：%d\n", min_size_node);

		int min = 0x7fffffff, max = 0;
		max_min_height(face_root, min, max);
		printf("Face min:%d, max:%d\n", min, max);

        if( max < 15)
		{
			break;
		}
	}

	build_registration_pair();

	std::cout<<"Build Pair Done!\n";

	if( BLENDING )
	{
		pre_blending_process(face_root);
		std::cout<<"Pre_blending_process Done!\n";
	}

	build_interpolation(face_root, m_mesh, target_mesh, 0.75);
		
	MyMesh::VertexIter v_iter, v_end( m_mesh.vertices_end());

	for (v_iter = m_mesh.vertices_begin(); v_iter != v_end; ++v_iter)
	{
		MyVector3f v = face_root->pts[ v_iter.handle().idx() ];
		m_mesh.set_point( v_iter, MyMesh::Point( v[0], v[1], v[2]) );
	}

	write_mesh_data( "interpolation_result.off");
}


void MsInterpolation::random_seed_point_helper( int size, int *idx )
{
	srand( time(NULL) );

	int curr_idx = 0;
	bool ret;

	while ( 1 )
	{
		int random_value = rand() % size;

		ret = false;
		for( int i = 0; i != curr_idx; i++)
		{
			if( random_value == idx[ i ] )
			{
				ret = true;
				break;
			}
		}

		if( ! ret )
		{
			idx[ curr_idx++ ] = random_value;
		}

		if( curr_idx == NR_MIN_PATCH_PER_LEVEL )
		{
			break;
		}
	}
}


void MsInterpolation::build_hierarchy_on_face()
{
	if( face_root == NULL)
	{
		face_root = new FaceNode();
	}else
	{
		printf("face root已经被创建了!\n");
		return ;
	}

	MyMesh::FaceIter f_iter, f_end( m_mesh.faces_end() );
	MyFaceHandle fh;

	for(f_iter = m_mesh.faces_begin(); f_iter != f_end; ++f_iter)
	{
		fh = f_iter.handle();

		{
			face_root->idx.push_back( fh.idx() );
		}
	}

	build_hierarchy_based_on_face(face_root, 0);
}


void MsInterpolation::build_hierarchy_based_on_face( FaceNode *subroot, int level /*= 0*/ )
{
	FaceNode *curr_node;
	MyFaceHandle curr_face_handle;

	int nr_faces = subroot->idx.size();
	int last_pos[NR_MIN_PATCH_PER_LEVEL];
	
	std::vector<int> face_flag;
	face_flag.assign(nr_faces, -1);

	std::vector<bool> need_iter;
	need_iter.assign(nr_faces, true);

	std::map<int, int> ridx;
	for(int i = 0; i < nr_faces; i++)
	{
		ridx.insert( std::pair<int, int>(subroot->idx[i], i) );
	}

	int seed_idx[NR_MIN_PATCH_PER_LEVEL];
	random_seed_point_helper(nr_faces, seed_idx);

	for(int i = 0; i < NR_MIN_PATCH_PER_LEVEL; i++)
	{
		int selected_idx = seed_idx[ i ];
		curr_node = subroot->next[i] = new FaceNode();

		curr_node->idx.push_back( subroot->idx[ selected_idx ] );
		face_flag[ selected_idx ] = i;

		last_pos[ i ] = 0;//设置起点
	}

	MyMesh::FaceIter face_iter;
	MyMesh::FaceFaceIter faceface_iter;
	MyFaceHandle close_face_handle;

	bool done;

	int start_pos, end_pos;// 定义iterator区间

	while( 1 )
	{
		done = true;
		for(int i = 0; i < nr_faces; i++)
		{
			if( face_flag[ i ] == -1 )
			{
				done = false;
				break;
			}
		}

		if( done )
		{
			break;// 此level结束
		}

		for(int i = 0; i < NR_MIN_PATCH_PER_LEVEL; i++)
		{

			curr_node = subroot->next[ i ];

			start_pos = last_pos[ i ];
			end_pos = curr_node->idx.size();

			for(int j = start_pos; j != end_pos; j++)
			{
				curr_face_handle = MyFaceHandle( curr_node->idx[ j ] );

				if( ! need_iter[ ridx[ curr_face_handle.idx() ] ] )
				{
					continue;
				}

				for( faceface_iter = m_mesh.ff_begin( curr_face_handle ); faceface_iter; ++faceface_iter )
				{
					close_face_handle = faceface_iter.handle();
					int close_idx = close_face_handle.idx();

					int ret = binary_search( subroot->idx.begin(), subroot->idx.end(), close_idx );
					if ( ret == false )
					{
						continue;// 该face不在上一层的patch中
					}

					int used = face_flag[ ridx[ close_idx ] ];
					if( used == -1 )
					{
						// 还未被使用
						curr_node->idx.push_back( close_idx );
						face_flag[ ridx[ close_idx ] ] = i;

					}else
					{
						if( used != i)
						{
							// 邻近patch的三角形
							curr_node->boundray.push_back( close_idx );
						}
					}

				}
				
			}
			last_pos[ i ] = end_pos;

		}

	}


    /*** 添加边缘 ***/
	for (int i = 0; i != NR_MIN_PATCH_PER_LEVEL; i++)
	{

		FaceNode *node_i = subroot->next[i];
		node_i->idx.insert( node_i->idx.end(), node_i->boundray.begin(), node_i->boundray.end() );
	}


	/*** 去重 ***/
	for (int i = 0; i != NR_MIN_PATCH_PER_LEVEL; i++)
	{

		FaceNode *node_i = subroot->next[i];

		sort( node_i->idx.begin(), node_i->idx.end() ); // 排序

		std::vector<int>::iterator iter;

		iter = unique(node_i->idx.begin(), node_i->idx.end() );
		node_i->idx.erase( iter, node_i->idx.end() );
	}

	level++;
	
	for(int i = 0; i < NR_MIN_PATCH_PER_LEVEL; i++)
	{
		int next_idx_len = subroot->next[ i ]->idx.size();
		
		if( next_idx_len < min_size_node )
		{
			min_size_node = next_idx_len;
		}

		if( next_idx_len > LEAF_NODE_MIN_NR )
		{
			build_hierarchy_based_on_face( subroot->next[ i ], level );
		}
	}
}


void MsInterpolation::release_facenode_help( FaceNode * &subroot )
{
	if( subroot == NULL )
	{
		return;
	}

	for(int i = 0; i < NR_MIN_PATCH_PER_LEVEL; i++)
	{
		release_facenode_help( subroot->next[ i ] );
	}

	delete subroot;
	subroot = NULL;
}

void MsInterpolation::max_min_height( FaceNode *subroot, int &min, int &max )
{
	if( subroot == NULL)
	{
		min = 0;
		max = 0;
		return ;
	}

	int mins[NR_MIN_PATCH_PER_LEVEL];
	int maxs[NR_MIN_PATCH_PER_LEVEL];

	for(int i = 0; i < NR_MIN_PATCH_PER_LEVEL; i++)
	{	
		mins[i] = 0x7fffffff;
		maxs[i] = 0;
		max_min_height(subroot->next[i], mins[i], maxs[i]);
	}

	for(int i = 0; i < NR_MIN_PATCH_PER_LEVEL; i++)
	{
		if( mins[i] < min )
		{
			min = mins[i];
		}

		if( maxs[i] > max )
		{
			max = maxs[i];
		}
	}

	min += 1;
	max += 1;
}


void MsInterpolation::build_registration_pair( FaceNode *subroot, int level  )
{

	if( subroot == NULL )
	{
		return;
	}

	bool leaf_node = true;
	for (int i = 0; i != NR_MIN_PATCH_PER_LEVEL; i++)
	{
		if( subroot->next[i] )
		{
			leaf_node = false;
			break;
		}
	}

	if( leaf_node )
	{
		return;
	}

	/*** 计算该节点的点index  ***/
	std::set<int> set_idx[NR_MIN_PATCH_PER_LEVEL];
	for (int i = 0; i != NR_MIN_PATCH_PER_LEVEL; i++)
	{	
		FaceNode *node_i = subroot->next[i];

		MyMesh::FaceVertexIter fv_iter;
		MyFaceHandle fh;

		int size = node_i->idx.size();
		for (int k = 0; k != size; k++)
		{
			fh = MyFaceHandle( node_i->idx[k] );
			for (fv_iter = m_mesh.fv_begin(fh); fv_iter; ++fv_iter)
			{
				set_idx[i].insert( fv_iter.handle().idx() );
			}
		}


		node_i->pts_index.resize( set_idx[i].size() );
		std::copy( set_idx[i].begin(), set_idx[i].end(), node_i->pts_index.begin() );

		for (int j = 0;  j != set_idx[i].size(); j++)
		{
			node_i->r_idx.insert( std::pair<int, int>(node_i->pts_index[j], j) );
		}

	}

	/***计算交集***/
	for (int i = 0; i != NR_MIN_PATCH_PER_LEVEL; i++)
	{
		for (int j = i + 1; j != NR_MIN_PATCH_PER_LEVEL; j++)
		{
			std::vector<int> ret;
			std::vector<int>::iterator iter;

			ret.resize( set_idx[i].size() );

			iter = std::set_intersection(set_idx[i].begin(), set_idx[i].end(), set_idx[j].begin(), set_idx[j].end(), ret.begin() );

			ret.resize( iter - ret.begin() );

			if( ret.empty() )
			{
				continue;
			}


			Pair p = Pair(i, j);

			subroot->P.push_back( p );
			subroot->pl.push_back( ret );
		}

	}

	level++;

	for (int i = 0; i !=  NR_MIN_PATCH_PER_LEVEL; i++)
	{
		build_registration_pair(subroot->next[i], level);
	}

}

void MsInterpolation::build_registration_pair()
{
	MyMesh::VertexIter iter, v_end( m_mesh.vertices_end() );

	face_root->pts_index.reserve( m_mesh.n_vertices() );

	int k = 0;
	for (iter = m_mesh.vertices_begin(); iter != v_end; ++iter)
	{
		face_root->pts_index.push_back( iter.handle().idx());
		face_root->r_idx.insert( std::pair<int, int>(iter.handle().idx(), k++) );
	}

	build_registration_pair( face_root );
}


void MsInterpolation::build_interpolation( FaceNode *subroot, MyMesh &src_mesh, MyMesh &target_mesh, double t )
{
	bool leaf_node = true;
 
	for (int i = 0; i != NR_MIN_PATCH_PER_LEVEL; i++)
	{
		if( subroot->next[i] != NULL)
		{
			leaf_node = false;
			break;
		}
	}

	/***如果是叶节点***/
	if( leaf_node )
	{
		leaf_node_interpolation(t, subroot, src_mesh,target_mesh);
		return;
	}

	for (int i = 0; i != NR_MIN_PATCH_PER_LEVEL; i++)
	{
		build_interpolation(subroot->next[i], src_mesh, target_mesh, t);
	}

	/***进行配准***/
	int M = NR_MIN_PATCH_PER_LEVEL;
	int P = subroot->P.size();

	std::vector<PointList> pl;

	for (int i = 0; i != P; i++)
	{
		int a = subroot->P[i].a;
		int b = subroot->P[i].b;

		FaceNode *node_a = subroot->next[a];
		FaceNode *node_b = subroot->next[b];
		
		std::vector<int> *curr_p = &subroot->pl[i];
		int size = curr_p->size();
		
		PointList pts;

		for (int j = 0; j != size; j++)
		{
			int idx = (*curr_p)[j];//共有的点索引
			PairVertex pv;

			pv.a = node_a->pts[ node_a->r_idx[idx] ];
			pv.b = node_b->pts[ node_b->r_idx[idx] ];


			pts.push_back( pv );
		}

		pl.push_back( pts );
	}

	MultiRegistration mr(M, P, &subroot->P, &pl);
	mr.init();

	std::vector<MyMatrix3f> R;
	std::vector<MyVector3f> T;

	mr.get_R_and_T(R, T);

	/***转换到统一的坐标系中***/
	for (int i = 0; i != NR_MIN_PATCH_PER_LEVEL; i++)
	{
		int size = subroot->next[i]->pts.size();

		for (int j = 0; j != size; j++)
		{
			subroot->next[i]->pts[j] = R[i] * subroot->next[i]->pts[j] + T[i];
		}

	}

	/***分配空间***/
	int size = subroot->pts_index.size();
	subroot->pts.assign(size, MyVector3f(0, 0, 0) );

	if( BLENDING )
	{
		int m = 0;
		MyMatrixXf E;

		E.resize(subroot->M.rows(), 3);

		for (int i = 0; i != NR_MIN_PATCH_PER_LEVEL; i++)
		{
			FaceNode *node_i = subroot->next[i];
			int nr_edges = subroot->edges_vector[i].size();

			for (int j = 0; j !=  nr_edges; j++)
			{
				int g_from = subroot->edges_vector[i][j].a;
				int g_to = subroot->edges_vector[i][j].b;

				int n_from =node_i->r_idx[ g_from ];
				int n_to = node_i->r_idx[ g_to ];

				MyVector3f e = node_i->pts[n_from] - node_i->pts[n_to];

				int from = subroot->r_idx[ g_from ];
				int to = subroot->r_idx[ g_to ];

				assert( subroot->M(m, from) == 1);
				assert( subroot->M(m, to) == -1);

				E(m, 0) = e(0);
				E(m, 1) = e(1);
				E(m, 2) = e(2);
				m++;
			}
		}

		E(m, 0) = 0;
		E(m, 1) = 0;
		E(m, 2) = 0;

		Eigen::SimplicialCholesky<MySMatrixXf > solver;
		solver.compute( (subroot->M.transpose() * subroot->M) );

		MyMatrixXf result = solver.solve( subroot->M.transpose() * E );

		int pts_size = subroot->pts_index.size();

		for (int i = 0; i != pts_size; i++)
		{
			double x, y, z;

			if( fabs( result(i, 0) ) < 1e-10 )
			{
				x = 0;
			}else
			{
				x = result(i, 0);
			}

			if( fabs( result(i, 1) ) < 1e-10 )
			{
				y = 0;
			}else
			{
				y = result(i, 1);
			}

			if( fabs( result(i, 2) ) < 1e-10 )
			{
				z = 0;
			}else
			{
				z = result(i, 2);
			}

			MyVector3f vi = MyVector3f(x, y, z);
			subroot->pts[i] = vi;
		}

	}else
	{
		// 如果没有blending的过程，则平均
		std::vector<int> cnt;

		cnt.assign( size, 0);

		for (int i = 0; i != NR_MIN_PATCH_PER_LEVEL; i++)
		{
			FaceNode *node_i = subroot->next[i];
			int pts_size = node_i->pts_index.size();

			for (int j = 0; j != pts_size; j++)
			{
				int idx = node_i->pts_index[j];
				int r_idx = subroot->r_idx[ idx ];

				cnt[r_idx] += 1;

				subroot->pts[ r_idx ] += node_i->pts[j];
			}

		}

		/**** 对于重复的点的坐标取平均值 ***/
		for (int i = 0; i != size; i++)
		{

			subroot->pts[ i ]  /= cnt[i];

		}
	}
}

void MsInterpolation::pre_blending_process( FaceNode * subroot )
{
	bool leaf_node = true;

	for (int i = 0; i != NR_MIN_PATCH_PER_LEVEL; i++)
	{
		if( subroot->next[i])
		{
			leaf_node = false;
			break;
		}
	}

	if( leaf_node )
	{
		return ;
	}


	/*** 先统计子节点的边向量 ***/
	std::vector<Pair> *edges_vector = subroot->edges_vector;

	for (int i = 0; i != NR_MIN_PATCH_PER_LEVEL; i++)
	{
		FaceNode *node_i = subroot->next[i];
		int nr_face = node_i->idx.size();

		/*** 确保边向量起点索引小于终点索引 ***/
		for (int j = 0; j != nr_face; j++)
		{
			int face_idx = node_i->idx[j];

			MyMesh::FaceHalfedgeIter fe_iter;

			for (fe_iter = m_mesh.fh_begin( MyFaceHandle(face_idx) ); fe_iter; ++fe_iter  )
			{
				int from = m_mesh.from_vertex_handle( fe_iter ).idx();
				int to = m_mesh.to_vertex_handle( fe_iter ).idx();

				if( from > to )
				{
					std::swap( from, to);
				}

				edges_vector[i].push_back( Pair(from, to ) );
			}
		}
	}


	/*** 去重 ***/
	int nr_cond = 0; // 约束条件个数
	for (int i = 0; i != NR_MIN_PATCH_PER_LEVEL; i++)
	{
		std::sort(edges_vector[i].begin(), edges_vector[i].end());

		std::vector<Pair>::iterator end_iter;

		end_iter = std::unique(edges_vector[i].begin(), edges_vector[i].end() );

		edges_vector[i].erase( end_iter, edges_vector[i].end() );

		nr_cond += edges_vector[i].size();
	}

	int nr_vertex = subroot->pts_index.size();
	subroot->M.resize( nr_cond + 1, nr_vertex);

	subroot->M.setZero();// 必要的

	int m = 0; 
	for (int i = 0; i != NR_MIN_PATCH_PER_LEVEL; i++)
	{
		int pair_size = edges_vector[i].size();

		for (int j = 0; j != pair_size; j++)
		{
			int from = edges_vector[i][j].a;
			int to = edges_vector[i][j].b;

			//local index
			from = subroot->r_idx[from];
			to = subroot->r_idx[to];

	        subroot->M.insert(m, from) = 1;
			subroot->M.insert(m, to) = -1;

			m++;
		}
	}

	subroot->M.insert(m, 0) = 1;

	for (int i = 0; i != NR_MIN_PATCH_PER_LEVEL; i++)
	{
		pre_blending_process(subroot->next[i]);
	}
}
