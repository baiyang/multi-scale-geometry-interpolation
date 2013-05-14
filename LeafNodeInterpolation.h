#pragma once
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
using namespace OpenMesh;

#include <Eigen/Dense>


typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;
typedef OpenMesh::VertexHandle MyVertexHandle;
typedef OpenMesh::FaceHandle MyFaceHandle;
typedef Eigen::Vector3d MyVector3f;

typedef OpenMesh::HalfedgeHandle MyHalfedgeHandle;


class FaceNode;


/*** in_vect， normal 是单位长度是1的向量，cos_theta是cos值 ***/
void get_beta_vector(MyVector3f &in_vect, double cos_theta, MyVector3f &normal, MyVector3f &out_vect, int sign = 1.0);

/*** 对叶节点进行插值，插值后的数据存放在leaf_node的pts数据里 ***/
void leaf_node_interpolation(double t, FaceNode *leaf_node, MyMesh &src_mesh, MyMesh &target_mesh);

/*** x,y,z是逆时针 ***/
void get_triangle_normal(MyVector3f &x, MyVector3f &y, MyVector3f &z, MyVector3f &out_normal);