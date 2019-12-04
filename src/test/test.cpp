/*
 * test.cpp
 *
 *  Created on: Oct 24, 2019
 *      Author: teng
 */

#include <stdio.h>
#include <iostream>
#include <type_traits>
#include "util/util.h"
using namespace std;

int main(int argc, char **argv){


	int count = 0;
	for(int i =0;i<1000000;i++){
		if(hispeed::get_rand_sample(atoi(argv[1]))){
			count++;
		}
	}
	cout<<count<<endl;

}
