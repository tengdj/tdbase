/*
 * himesh_3dpro.cpp
 *
 *  Created on: Apr 10, 2023
 *      Author: teng
 */


#include "himesh.h"

namespace tdbase{

bool HiMesh::isProtruding(const std::vector<Halfedge_const_handle> &polygon) const
{
	if(polygon.size()<2){
		return true;
	}
	Point top(polygon[0]->opposite()->vertex()->point());
	vector<Point> rings;

	for(Halfedge_const_handle h:polygon){
		Point p(h->vertex()->point());
		rings.push_back(p);
	}

	// find one possible triangulation way that fulfill protruding
	bool is_recessing = false;
	// t is the starting vertex on the ring
	for(int t=0;t<rings.size()-1;t++){
	//for(int t=0;t<1;t++){

		is_recessing = false;
		// evaluate all the tetrahedrons
		for(int i=1;i<rings.size()-1;i++){
			// the random vector pointing from the bottom triangle to the top vertex
			double r1 = top.x()-rings[t].x();
			double r2 = top.y()-rings[t].y();
			double r3 = top.z()-rings[t].z();

			// calculate the normal vector of the bottom triangle
			//printf("3 0 %d %d|",i+1,i+2);
			double a1 = rings[t].x()-rings[(i+t)%rings.size()].x();
			double a2 = rings[t].y()-rings[(i+t)%rings.size()].y();
			double a3 = rings[t].z()-rings[(i+t)%rings.size()].z();
			double b1 = rings[t].x()-rings[(i+t+1)%rings.size()].x();
			double b2 = rings[t].y()-rings[(i+t+1)%rings.size()].y();
			double b3 = rings[t].z()-rings[(i+t+1)%rings.size()].z();

			// n*a=0 n*b=0
			double n1 = a2*b3-a3*b2;
			double n2 = a3*b1-a1*b3;
			double n3 = a1*b2-a2*b1;

			if(!global_ctx.counter_clock){
				n1 = -n1;
				n2 = -n2;
				n3 = -n3;
			}

			// calculate the angle between the normal and vector t->0
			double cosvalue = (r1*n1+r2*n2+r3*n3)/(sqrt(r1*r1+r2*r2+r3*r3)*sqrt(n1*n1+n2*n2+n3*n3));
			double angle = acos(cosvalue)*180/M_PI;
			// avoid some corner case, such that bigger than 90.5 instead of 90, increase the tolerance.
			if(angle>90.5){
				is_recessing = true;
			}
			//printf("%d\tangle: %f %f\n",t,cosvalue,angle);
		}
		// this vertex can be protruding
		if(is_recessing == false){
			break;
		}
	}

	// print the removed part into a single polyhedron for visualization
	if(is_recessing && false){
		printf("%d\n",is_recessing);
		printf("OFF\n");
		printf("%ld %ld 0\n",1+rings.size(),1+rings.size());
		printf("%f %f %f\n",top.x(),top.y(),top.z());
		for(Point p:rings){
			printf("%f %f %f\n",p.x(),p.y(),p.z());
		}

		if(global_ctx.counter_clock){
			for(int i=0;i<rings.size()-1;i++){
				printf("3 0 %d %d 0 255 0\n",i+1,i+2);
			}
			printf("3 0 %ld 1 0 255 0\n",rings.size());

			printf("%ld",rings.size());
			for(int i=rings.size();i>=1;i--){
				printf(" %d",i);
			}
		}else{
			for(int i=0;i<rings.size()-1;i++){
				printf("3 %d %d 0 0 255 0\n",i+2,i+1);
			}
			printf("3 1 %ld 0 0 255 0\n",rings.size());

			printf("%ld",rings.size());
			for(int i=1;i<=rings.size();i++){
				printf(" %d",i);
			}
		}

		printf(" 255 0 0\n");

		vector<Point> norms;

		for(int i=1;i<rings.size()-1;i++){
			// the random vector pointing from the bottom triangle to the top vertex
			double r1 = top.x()-rings[0].x();
			double r2 = top.y()-rings[0].y();
			double r3 = top.z()-rings[0].z();

			// calculate the normal vector of the bottom triangle
			//printf("3 0 %d %d|",i+1,i+2);
			double a1 = rings[0].x()-rings[(i+0)%rings.size()].x();
			double a2 = rings[0].y()-rings[(i+0)%rings.size()].y();
			double a3 = rings[0].z()-rings[(i+0)%rings.size()].z();
			double b1 = rings[0].x()-rings[(i+0+1)%rings.size()].x();
			double b2 = rings[0].y()-rings[(i+0+1)%rings.size()].y();
			double b3 = rings[0].z()-rings[(i+0+1)%rings.size()].z();

			// n*a=0 n*b=0
			double n1 = a2*b3-a3*b2;
			double n2 = a3*b1-a1*b3;
			double n3 = a1*b2-a2*b1;

			if(!global_ctx.counter_clock){
				n1 = -n1;
				n2 = -n2;
				n3 = -n3;
			}
			Point p(n1*100,n2*100,n3*100);
			norms.push_back(p);
			//printf("%d\tangle: %f %f\n",t,cosvalue,angle);
		}

		// the vertical lines
		printf("OFF\n");
		printf("%ld %ld 0\n",2*(rings.size()-2)+1,rings.size()-2);
		for(int i=1;i<rings.size()-1;i++){
			printf("%f %f %f\n",rings[i].x(),rings[i].y(),rings[i].z());
		}
		for(int i=1;i<rings.size()-1;i++){
			printf("%f %f %f\n",rings[i].x()+norms[i-1].x(),rings[i].y()+norms[i-1].y(),rings[i].z()+norms[i-1].z());
		}
		printf("%f %f %f\n",rings[0].x(),rings[0].y(),rings[0].z());
		for(int i=1;i<rings.size()-1;i++){
			if(i==1){
				printf("3 %d %ld %ld 255 0 0\n",i-1,i-1+rings.size()-2,2*(rings.size()-2));
			}else if(i==2){
				printf("3 %d %ld %ld 0 255 0\n",i-1,i-1+rings.size()-2,2*(rings.size()-2));
			}else if(i==3){
				printf("3 %d %ld %ld 0 0 255\n",i-1,i-1+rings.size()-2,2*(rings.size()-2));
			}else{
				printf("3 %d %ld %ld\n",i-1,i-1+rings.size()-2,2*(rings.size()-2));
			}
		}
		exit(0);
	}
	// no recessing point
	return !is_recessing;
}

/*
 *
 * Test whether a vertex is protruding
 * Core function for 3DPro
 *
 * */

void HiMesh::profileProtruding(){
	int protruding = 0;
	int recessing = 0;
	for(HiMesh::Vertex_iterator vit = vertices_begin(); vit!=vertices_end(); ++vit){
		std::vector<Halfedge_const_handle> heh_oneRing;
		heh_oneRing.reserve(vit->vertex_degree());
		//vh_oneRing.push_back(v);
		Halfedge_around_vertex_const_circulator hit(vit->vertex_begin()), end(hit);
		do
		{
		  heh_oneRing.push_back(hit->opposite());
		}while(++hit != end);

		if(isProtruding(heh_oneRing)){
			protruding++;
		}else{
			recessing++;
		}
	}
	printf("%d %d %f\n",protruding,recessing,protruding*100.0/(protruding+recessing));
}

}
