#include "util.h"
#include <iostream>
#include <map>

using namespace std;

int main(int argc, char **argv) {
    std::fstream file1(argv[1], std::ios_base::in);
    int total1 = 0;
    int matched = 0;
    int total2 = 0;

    int source;
    int target;
    map<int, map<int,int>> results;
    while (file1 >> source){
        file1 >> target;
        if (results.find(source)==results.end()) {
            map<int, int> lst;
            results[source] = lst;
        }
        results[source][target] = 1;
        total1++; 
    }
    file1.close();

    std::fstream file2(argv[2], std::ios_base::in);
    while (file2 >> source) {
        file2 >> target;
        assert(results.find(source) != results.end());
        if (results[source].find(target)!= results[source].end()) {
            matched++;
        }
        total2++;
    }
    file2.close();

    printf("%d,%d,%d\n", total1, total2, matched);
}