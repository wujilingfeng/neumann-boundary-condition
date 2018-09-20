#pragma once
#include <iostream>
#include <list>
#include<stack>
#include<nanogui\screen.h>
#include<nanogui\formhelper.h>
#include<igl\readOFF.h>
#include<igl\viewer\Viewer.h>
#include<igl\jet.h>
//#include<Eigen/Core>
#include<Eigen/IterativeLinearSolvers>
//#include <Eigen/SparseLU>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include<OpenMesh\Core\Mesh\TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/System/config.h>
#include <OpenMesh/Core/Mesh/Status.hh>
#include<GenEigsSolver.h>
#include<MatOp/SparseGenMatProd.h>
#include "igl/cotmatrix.h"
#ifndef PI
#define PI 3.1415926
#endif 
struct MyTraits :public OpenMesh::DefaultTraits
{
	VertexTraits{ double ux = 0,uy = 0; int index = 0; };
	EdgeTraits{ double weight = 0; };
	HalfedgeTraits{  };
	FaceTraits{  };
};
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> TMyMesh;
typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits> PMesh;
class Harmonic {
public:
	Harmonic(TMyMesh &me,igl::viewer::Viewer &vie):mesh(me),viewer(vie)
	{
		printf_s("create harmonic\r\n");
	}
	~Harmonic()
	{}
	void init()
	{

		int i = 0, j = 0;
		for (TMyMesh::VertexIter vter=mesh.vertices_begin();vter!=mesh.vertices_end();vter++)
		{
		//	printf_s("kaishi\r\n");
			if (mesh.is_boundary(*vter))
			{
				mesh.data(*vter).index = j;
				j++;
			}
			else
			{
				mesh.data(*vter).index = i;
				i++;
			}
		}
		
		//printf_s("%d  %d",mesh.n_vertices(),j);
		H.resize(i,i);
		H.setZero();
	}
	void init_weight()
	{
		Eigen::SparseMatrix<double> L;
		igl::cotmatrix(V, F, L);
		TMyMesh::HalfedgeHandle he;
		printf_s("kaishi\r\n");
		for (TMyMesh::EdgeIter eter = mesh.edges_begin(); eter != mesh.edges_end(); eter++)
		{
			//printf_s("kaishi\r\n");
			he = mesh.s_halfedge_handle(*eter, 0);
			mesh.data(*eter).weight = L.coeff(mesh.to_vertex_handle(he).idx(), mesh.from_vertex_handle(he).idx());
			/*if (mesh.is_boundary(*eter))
			{
			printf_s(" %lf   %lf     %lf\r\n",mesh.data(*eter).weight,1/tan(mesh.calc_sector_angle(mesh.next_halfedge_handle(he))), 1 / tan(mesh.calc_sector_angle(mesh.prev_halfedge_handle(he))));
			}*/
			//printf_s("%f,  %f\r\n", mesh.data(*eter).weight, weight(*eter));
			if (L.coeff(mesh.to_vertex_handle(he).idx(), mesh.from_vertex_handle(he).idx()) != L.coeff(mesh.from_vertex_handle(he).idx(), mesh.to_vertex_handle(he).idx()))
			{
				printf_s("´íÎó\r\n");
			}
		}
	
	}
	void construct_LL()
	{
		Eigen::MatrixXd v(H.rows(),mesh.n_vertices()-H.rows());
		
		//int i = 0, j = 0;
		TMyMesh::HalfedgeHandle he; TMyMesh::VertexHandle v1, v2;
		v.setZero();
		LL = v;
		for (TMyMesh::EdgeIter eter=mesh.edges_begin();eter!=mesh.edges_end();eter++)
		{
			//printf_s("he\r\n");
			if (mesh.is_boundary(*eter))
			{
			}
			else
			{
				
				he = mesh.s_halfedge_handle(*eter,0);
				v1 = mesh.to_vertex_handle(he);
				v2 = mesh.from_vertex_handle(he);
				if (!mesh.is_boundary(v1)&&!mesh.is_boundary(v2))
				{
					H.coeffRef(mesh.data(v1).index, mesh.data(v2).index) =- mesh.data(*eter).weight;
					H.coeffRef(mesh.data(v2).index, mesh.data(v1).index) =- mesh.data(*eter).weight;
					H.coeffRef(mesh.data(v2).index, mesh.data(v2).index) += mesh.data(*eter).weight;
					H.coeffRef(mesh.data(v1).index, mesh.data(v1).index) += mesh.data(*eter).weight;
				}
				else if (mesh.is_boundary(v1)&&!mesh.is_boundary(v2))
				{
					//printf_s("cunzai\r\n");
					H.coeffRef(mesh.data(v2).index, mesh.data(v2).index) += mesh.data(*eter).weight;
					v.coeffRef(mesh.data(v2).index,mesh.data(v1).index)=mesh.data(*eter).weight;
				}
				else if (!mesh.is_boundary(v1)&&mesh.is_boundary(v2))
				{
					H.coeffRef(mesh.data(v1).index, mesh.data(v1).index) += mesh.data(*eter).weight;
					v.coeffRef(mesh.data(v1).index, mesh.data(v2).index) = mesh.data(*eter).weight;
				}
				else
				{
					
					printf_s("cuowullll\r\n");
				}
				
				//H.coeffRef()
			}
		}
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>llt;
		llt.compute(H);
		for (int j=0;j<v.cols();j++)
		{
			LL.col(j) = llt.solve(v.col(j));
		}
	}
	void disk_boudary()
	{
		TMyMesh::VertexHandle v; TMyMesh::VertexIter vter = mesh.vertices_begin();
		for (;vter!=mesh.vertices_end();vter++)
		{
			if (mesh.is_boundary(*vter))
			{
				v = *vter;
				break;
			}
		}
		Eigen::VectorXd temp_v(LL.cols());
		temp_v.setZero(); int i = 0;
		TMyMesh::HalfedgeHandle he;
		do {
			he = mesh.halfedge_handle(v);
			temp_v.coeffRef(i)=mesh.calc_edge_length(he)+temp_v.coeff((i-1+LL.cols())%LL.cols());
			v = mesh.to_vertex_handle(he);
			i++;
		//	printf_s("%d\r\n",i);
		} while (v!=(*vter));
		double sum = temp_v.coeff(i - 1);
		i = 0; double theta = 0;
		do {
			 theta = 2 * PI*temp_v.coeff(i) / sum;
			mesh.data(v).ux = cos(theta);
			mesh.data(v).uy = sin(theta);
			he = mesh.halfedge_handle(v);
			//temp_v.coeffRef(i) = mesh.calc_edge_length(he) + temp_v.coeff((i - 1 + LL.cols()) % LL.cols());
			v = mesh.to_vertex_handle(he);
			i++;
			//printf_s("%d\r\n", i);
		} while (v != (*vter));
		

	}
	void natural_confromal_boudary()
	{
		Eigen::MatrixXd temp_H(2*LL.cols(),2*mesh.n_vertices()),temp_L(2*mesh.n_vertices(),2*LL.cols());
		Eigen::VectorXd vx(2 * LL.cols()), x;
		vx.setZero();
		vx.coeffRef(0) = 0; vx.coeffRef(0 + LL.cols()) = 0;
		vx.coeffRef(1) = 0; vx.coeffRef(1 + LL.cols()) = 0.2;
		temp_H.setZero(); temp_L.setZero();	
		int n = mesh.n_vertices(),int_v=LL.rows();
		printf_s("%d  %d\r\n",n,LL.rows()+LL.cols());
		for (TMyMesh::VertexIter vter=mesh.vertices_begin();vter!=mesh.vertices_end();vter++)
		{
			
			if (mesh.is_boundary(*vter))
			{
				if (mesh.data(*vter).index == 0 || mesh.data(*vter).index == 1)
				{
					temp_H.coeffRef(mesh.data(*vter).index, int_v + mesh.data(*vter).index) = 1;
					temp_H.coeffRef(LL.cols() + mesh.data(*vter).index, n + int_v + mesh.data(*vter).index) = 1;
					
				}
				else
				{
					/*for (TMyMesh::VertexVertexCWIter vvter=mesh.vv_cwbegin(*vter);vvter.is_valid();vvter++)
					{
						temp_H.coeffRef(mesh.data(*vter).index,mesh.data(*vter).index)+=
					}*/
					//int temp_d = 0;
					for (TMyMesh::VertexOHalfedgeCWIter vhter=mesh.voh_cwbegin(*vter);vhter.is_valid();vhter++)
					{
						temp_H.coeffRef(mesh.data(*vter).index, int_v+mesh.data(*vter).index) += 2*mesh.data(mesh.edge_handle(*vhter)).weight;
					    temp_H.coeffRef(LL.cols()+mesh.data(*vter).index, n+int_v + mesh.data(*vter).index)+= 2*mesh.data(mesh.edge_handle(*vhter)).weight;
						TMyMesh::VertexHandle v=mesh.to_vertex_handle(*vhter);
						if (mesh.is_boundary(v))
						{
							temp_H.coeffRef(mesh.data(*vter).index, int_v + mesh.data(v).index) = -2*mesh.data(mesh.edge_handle(*vhter)).weight;
							temp_H.coeffRef(LL.cols()+mesh.data(*vter).index, n+int_v + mesh.data(v).index) = -2*mesh.data(mesh.edge_handle(*vhter)).weight;
						
							/*temp_H.coeffRef(mesh.data(*vter).index, n + int_v + mesh.data(v).index) = pow(-1,temp_d);
							temp_H.coeffRef(LL.cols()+mesh.data(*vter).index, int_v + mesh.data(v).index) =-pow(-1,temp_d);
							temp_d++;*/
						}
						else
						{
							temp_H.coeffRef(mesh.data(*vter).index,  mesh.data(v).index) = -2*mesh.data(mesh.edge_handle(*vhter)).weight;
							temp_H.coeffRef(LL.cols() + mesh.data(*vter).index, n + mesh.data(v).index) = -2*mesh.data(mesh.edge_handle(*vhter)).weight;

						}
					
						
					}
					TMyMesh::VertexOHalfedgeCCWIter vhter = mesh.voh_ccwbegin(*vter);
					TMyMesh::VertexHandle v = mesh.to_vertex_handle(*vhter);
					if (!mesh.is_boundary(v))
					{
						printf_s("´íÎó\r\n");
					}
					temp_H.coeffRef(mesh.data(*vter).index, n + int_v + mesh.data(v).index) = 1;
					temp_H.coeffRef(LL.cols() + mesh.data(*vter).index, int_v + mesh.data(v).index) = -1;
					vhter++;
					v = mesh.to_vertex_handle(*vhter);
					if (!mesh.is_boundary(v))
					{
						printf_s("´íÎó\r\n");
					}
					temp_H.coeffRef(mesh.data(*vter).index, n + int_v + mesh.data(v).index) = -1;
					temp_H.coeffRef(LL.cols() + mesh.data(*vter).index, int_v + mesh.data(v).index) = 1;
					//printf_s("%d\r\n",temp_d);
				}


			}
		}
		temp_L.block(0,0,LL.rows(),LL.cols())=LL;
		temp_L.block(LL.rows(),0,LL.cols(),LL.cols())=Eigen::MatrixXd::Identity(LL.cols(), LL.cols());
		temp_L.block(n, LL.cols(), LL.rows(), LL.cols())=LL;
		temp_L.block(n+LL.rows(),LL.cols(), LL.cols(), LL.cols()) = Eigen::MatrixXd::Identity(LL.cols(), LL.cols());
		Eigen::MatrixXd temp_M = temp_H * temp_L;
		Eigen::SparseMatrix<double >temp_MM(2 * LL.cols(), 2 * LL.cols());
		for (int i=0;i<2 * LL.cols();i++)
		{
			for(int j=0;j<2 * LL.cols();j++)
			{
				temp_MM.coeffRef(i, j) = temp_M.coeff(i, j);
				//printf_s("%lf\r\n",temp_M.coeff(i,j));
			}
		}
		for (int i = 0; i<2 * LL.cols(); i++)
		{
			for (int j = 0; j<2 * LL.cols(); j++)
			{
				temp_MM.coeffRef(i, j) = temp_M.coeff(i, j);
				//printf_s("%lf\r\n",temp_M.coeff(i,j));
			}
		}
		Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver;
		//Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::AMDOrdering< int >> LU;
		//Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>llt;
		//LU.compute(temp_MM);
		//x = LU.solve(vx);
		solver.compute(temp_MM);
		x = solver.solve(vx);
		//x = temp_M.lu().solve(vx);
		printf_s("%lf\r\n", vx.coeffRef(1 + LL.cols()));
		for (TMyMesh::VertexIter vter = mesh.vertices_begin(); vter != mesh.vertices_end(); vter++)
		{
			if (mesh.is_boundary(*vter))
			{
				mesh.data(*vter).ux = x.coeff(mesh.data(*vter).index);
				mesh.data(*vter).uy = x.coeff(LL.cols()+mesh.data(*vter).index);
				printf_s("%lf    %lf\r\n", vx.coeff(mesh.data(*vter).index), vx.coeff(LL.cols()+mesh.data(*vter).index));
			}
		}
		//temp_M.sparseView();
		//Eigen::MatrixXd::Identity
	}
	void map()
	{
	
		Eigen::VectorXd ux(LL.cols()), uy(LL.cols()),x(LL.cols()),y(LL.cols()); int i = 0;
		for (TMyMesh::VertexIter vter=mesh.vertices_begin();vter!=mesh.vertices_end();vter++)
		{
			if (mesh.is_boundary(*vter))
			{
				ux.coeffRef(i) = mesh.data(*vter).ux;
				uy.coeffRef(i) = mesh.data(*vter).uy;
				i++;
			}
		}
		x = LL * ux;
		y = LL * uy;
		//i = 0;
		for (TMyMesh::VertexIter vter = mesh.vertices_begin(); vter != mesh.vertices_end(); vter++)
		{
			if (mesh.is_boundary(*vter))
			{
			}
			else
			{
				mesh.data(*vter).ux = x.coeff(mesh.data(*vter).index);
				mesh.data(*vter).uy = y.coeff(mesh.data(*vter).index);
				i++;
			}
		}
	}
	void draw_edge()
	{
		printf_s("begin drawedge\r\n");
		Eigen::RowVector3d v1, v2; TMyMesh::VertexHandle vv1, vv2;
		
		TMyMesh::HalfedgeHandle he;
		for (TMyMesh::EdgeIter eter = mesh.edges_begin(); eter != mesh.edges_end(); eter++)
		{
			he = mesh.s_halfedge_handle(*eter, 0);
			vv1 = mesh.to_vertex_handle(he);
			vv2 = mesh.from_vertex_handle(he);
			v1.coeffRef(0) = mesh.data(vv1).ux; v2.coeffRef(0) = mesh.data(vv2).ux;
			v1.coeffRef(1) = mesh.data(vv1).uy; v2.coeffRef(1) = mesh.data(vv2).uy;
			v1.coeffRef(2) = 0;                 v2.coeffRef(2) = 0;
			
			if (mesh.is_boundary(vv1)&&mesh.is_boundary(vv2))
			{
				if (mesh.data(vv1).index==0||mesh.data(vv2).index==0|| mesh.data(vv1).index == 1|| mesh.data(vv2).index == 1)
				{
					viewer.data.add_edges(1.011*v1, 1.011*v2, Eigen::RowVector3d(1, 0, 0));
				}
				
				viewer.data.add_edges(v1, v2, Eigen::RowVector3d(0.5, 1, 0));
			}
			
			
		}
	}
	void draw()
	{
		Eigen::MatrixXd temp_V = V;
		for (TMyMesh::VertexIter vter=mesh.vertices_begin();vter!=mesh.vertices_end();vter++)
		{
			temp_V.coeffRef((*vter).idx(), 0) = mesh.data(*vter).ux;
			temp_V.coeffRef((*vter).idx(), 1) = mesh.data(*vter).uy;
			temp_V.coeffRef((*vter).idx(), 2) = 0;
		}
	
		viewer.data.set_mesh(temp_V,F);
		viewer.core.align_camera_center(temp_V,F);
	}
	void set_uv(Eigen::MatrixXd &V_uv)
	{
		V_uv.resize(mesh.n_vertices(), 2);
		V_uv.setZero();
		for (TMyMesh::VertexIter vter = mesh.vertices_begin(); vter != mesh.vertices_end(); vter++)
		{
			//V_uv.coeffRef((*vter).idx(), 0) = fabs(sin(mesh.data(vter).v)) / 3.0;
			V_uv.coeffRef((*vter).idx(), 0) = mesh.data(vter).ux* 40;
			//V_uv.coeffRef((*vter).idx(), 1) = mesh.data(vter).v;
			//V_uv.coeffRef((*vter).idx(), 1) = fabs(sin(mesh.data(vter).u))/3.0;
			V_uv.coeffRef((*vter).idx(), 1) = mesh.data(vter).uy * 40;
			printf_s("u: %f\r\n", mesh.data(vter).uy);
		}

	}
	
	
	
	void mywrite(const char* filename)
	{
		//ofstream myfile(filename);
		FILE *myfile = fopen(filename, "w");
		if (myfile == NULL)
		{
			printf("Ê§°Ü\r\n");
			return;
		}
		/*	for (list<VertexHandle>::iterator iter = vlist.begin(); iter != vlist.end(); iter++)
		{
		fprintf(myfile, "%d ", (*iter).idx());
		}*/
		// myfile << "libo hello";
		fclose(myfile);
	}
	void myread(const char* filename)
	{
		//ifstream myfile(filename);
		// string t;
		//myfile >> t;
		// printf(t.c_str());
		FILE *myfile = fopen(filename, "r");
		if (myfile == NULL)
		{
			printf("Ê§°Ü\r\n");
			return;
		}
		int a;
		while (fscanf(myfile, "%d", &a) != EOF)
		{
			vlist.push_back(mesh.vertex_handle(a));
			printf("%d \r\n", a);
		}
		/* for (list<VertexHandle>::iterator iter = vlist.begin(); iter != vlist.end(); iter++)
		{
		fprintf(myfile, "%d ", (*iter).idx());
		}*/
		// myfile << "libo hello";
		fclose(myfile);

	}
	void test()
	{
		TMyMesh::VertexHandle v = mesh.vertex_handle(0);
		TMyMesh::VertexOHalfedgeCCWIter vhter = mesh.voh_ccwbegin(v);
		printf_s("%d\r\n",mesh.from_vertex_handle(mesh.prev_halfedge_handle(*vhter)));
		for (;vhter.is_valid();vhter++)
		{
			printf_s("%d \r\n",mesh.to_vertex_handle(*vhter).idx());
		}
	
	}
	
	Eigen::SparseMatrix<double > H;
	Eigen::MatrixXd V,LL;
	Eigen::MatrixXi F;
	std::list<TMyMesh::VertexHandle> vlist;
protected:
	igl::viewer::Viewer &viewer;
	TMyMesh & mesh;


};