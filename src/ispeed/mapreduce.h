/* 
 * This header contains relevant function for interaction between Hadoop MapReduce and CPP code
 * */
#ifndef MAP_REDUCE_CPP_HPP
#define MAP_REDUCE_CPP_HPP

#include <cstring>
using namespace std;

// MapReduce constants
const string JAR_FILE_NAME = "hadoop-streaming.jar";
// relative path with respect to the lib path
const string CUSTOM_JAR_REL_PATH = "../lib/myCustomLibs.jar";
const string PARTITION_FILE_NAME = "partfile.idx";


string update_ld_lib_path();
string getHadoopJarPath();
string getHdfsCmdPath();
string getHadoopCmdPath();
string hdfs_str_result(string programpath, vector<string> &arr_args);
bool hdfs_check_data(string hadoopcmdpath, string input_path);
bool hdfs_cat(string hadoopcmdpath, string path, int outputfd = -1);
long hdfs_get_size(string programpath, string inputpath);
bool hdfs_delete(string hadoopcmdpath, string path);
bool hdfs_move(string programpath, string source, string destination);
bool hdfs_put(string programpath, string source, string destination);
pid_t execute_command(string programpath,vector<string> &strargs, int outputfd = -1);


#endif

