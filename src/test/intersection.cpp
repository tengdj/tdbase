/*
 * intersection.cpp
 *
 *  Created on: Dec 2, 2019
 *      Author: teng
 */

#include "../geometry/geometry.h"

using namespace hispeed;
using namespace std;

int main(int argc, char **argv)
{
	struct timeval start = get_cur_time();
	float PS[9];
	float QS[9];
	int sum = 0;
	for(int i=0; i<10000; i++){
		for(int j=0; j<3; j++){
			for(int k=0; k<3; k++){
				PS[j*3+k] = drand48();
				QS[j*3+k] = drand48();
				if(j!=0){
					PS[j*3+k] -= PS[k];
					QS[j*3+k] -= QS[k];
				}
			}
		}
		sum += TriInt(PS, QS);
	}
	report_time("computation", start);
	cout<<sum<<endl;

}


