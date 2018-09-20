#include "ncm.h"
TMyMesh mesh;
igl::viewer::Viewer viewer;
Harmonic har(mesh,viewer);

int main()
{

	igl::readOFF("e:/off/face1.off",har.V,har.F);
	OpenMesh::IO::Options opt;
	if (!OpenMesh::IO::read_mesh(mesh, "e:/off/face1.off", opt))
	{
		printf_s("∂¡»° ß∞‹\r\n");
		return 0;
	}
	//TMyMesh::VertexHandle v = mesh.vertex_handle(0);
	
	
	/*TMyMesh::VertexHandle v=mesh.vertex_handle(0);
	printf_s("%d\r\n",mesh.from_vertex_handle(mesh.halfedge_handle(v)).idx());*/
	har.init();
	har.init_weight();
	har.construct_LL();
	//har.test();
	har.natural_confromal_boudary();
	//har.disk_boudary();
	har.map();
	Eigen::MatrixXd V_uv;
	har.set_uv(V_uv);
	viewer.data.set_uv(V_uv);
	viewer.core.show_texture = true;
    //har.draw();
	//har.draw_edge();
	viewer.data.set_mesh(har.V,har.F);
	viewer.launch();

	return 0;
}