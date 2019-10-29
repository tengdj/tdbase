/*
 * test.cpp
 *
 *  Created on: Oct 24, 2019
 *      Author: teng
 */

#include <stdio.h>
#include <iostream>
using namespace std;

float TriDist(float S[3][3], float T[3][3]){
	return 0.0;
}

int main(int argc, char **argv){

	int *S = new int[9];
	for(int i=0;i<9;i++){
		S[i]=i;
	}

	for(int i=0;i<3;i++){
		int *curs = S+i*3;
		for(int j=0;j<3;j++){
			cout<<curs[j]<<endl;
		}
	}
}


