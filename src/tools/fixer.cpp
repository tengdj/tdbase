
#include "../mymesh/mymesh.h"
using namespace std;
using namespace hispeed;

int main(int argc, char **argv){
//	Polyhedron *poly = new Polyhedron();
//	poly->load(argv[1]);
//	vector<Polyhedron *> ind = poly->depart();
//	char path[256];
//	for(int i=0;i<ind.size();i++){
//		if(ind[i]->points.size()>100){
//			ind[i]->evaluate();
//			ind[i]->remove_redundant();
//			ind[i]->merge_vertex();
//			ind[i]->fill_holes();
//			//ind[i]->evaluate();
//			//ind[i]->fill_holes();
//			if(argc>2){
//				sprintf(path,"%s/fixed_%ld_%ld.OFF",argv[2],ind[i]->points.size(),ind[i]->faces.size());
//			}else{
//				sprintf(path,"/gisdata/fixed_%ld_%ld.OFF",ind[i]->points.size(),ind[i]->faces.size());
//			}
//			ind[i]->dumpto(path);
//		}
//		delete ind[i];
//	}
//	if(argc>2){
//		printf("%ld components are dumped to %s\n",ind.size(),argv[2]);
//	}
//	ind.clear();
//	delete poly;
	string str = hispeed::read_file(argv[1]);
	Polyhedron *poly = new Polyhedron();
	poly->parse(str.c_str(), str.size());
	poly->remove_redundant();
	poly->print();

}
