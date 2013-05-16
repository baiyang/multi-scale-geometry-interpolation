#include <OpenMesh/Core/IO/MeshIO.hh>
#include "MsInterpolation.h"

#include "LeafNodeInterpolation.h"

#include <iostream>
using namespace std;


int main()
{
	MsInterpolation ms;
	std::string src = "./horse-01.obj";
	std::string target = "./horse-02.obj";

	ms.read_mesh_data( src, target );
	ms.test();
	return 0;
}
