/* 
 * This is the framework controller program. It manages calls to MapReduce and holds results, parameters.
 * This serves as the basic front-end interface
 * */

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <map>
#include <cstdlib>
#include <getopt.h>
#include <time.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/fcntl.h>
#include <fcntl.h>

#include <boost/program_options.hpp>
#include "iutil.h"
#include "mapreduce.h"

namespace po = boost::program_options;
using namespace std;

// Combined framework variables
struct framework_vars {

	// MapReduce-related parameter
	std::string binary_prefix;  // absolute path to the folder contains the binary executables
	std::string streaming_path; // Absolute path to the hadoop streaming jar
	std::string hadoopcmdpath; // Absolute path to the command line hadoop executable
	std::string hdfscmdpath; // Absolute path to the command line hadoop executable

	int numreducers = 1; // number of reducer

	// Input/output data variables
	std::string input_path_1;
	std::string input_path_2;
	std::string output_path;

	//for resque
	std::string predicate="st_intersects";
	int decomp_lod = 100;
};

bool extract_params(int argc, char **argv, struct framework_vars &fr_vars){
	string binpath;
	string querytype;
	try {
		po::options_description desc("Options");
		desc.add_options()
			// common parameters
			("help,h", "This help message")
			("binpath",po::value<string>(&binpath)->required(), "path to the binary executables")
			("outputpath,o", po::value<string>(&fr_vars.output_path)->required(), "Output path")
			("numreducers,n", po::value<int>(&fr_vars.numreducers), "The number of reducers")
			("input1,a", po::value<string>(&fr_vars.input_path_1)->required(), "HDFS file path to data set 1")
			("input2,b", po::value<string>(&fr_vars.input_path_2)->required(), "HDFS file path to data set 2")
			// resque
			("lod,l", po::value<int>(&fr_vars.decomp_lod) , "Decompression LOD. (0, 100]. Default is 100.")
			("predicate,p", po::value<string>(&fr_vars.predicate)->required(), "Predicate for spatial join and nn queries "
					"[ st_intersects | st_touches | st_crosses | st_contains | st_adjacent | st_disjoint "
					"| st_equals | st_dwithin | st_within | st_overlaps | st_nn_voronoi | st_nn_rtree ] ")
			;
		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);

		if (vm.count("help")) {
			cerr << desc << endl;
			return false;
		}
		if(get_join_predicate(fr_vars.predicate.c_str())==ST_ERROR){
			cerr << desc << endl;
			cerr << "type of predict should be set properly"<<endl;
			return false;
		}

		remove_slash(fr_vars.output_path);
		/* Update environment variables with HADOOP_HOME defined*/
		fr_vars.hadoopcmdpath = getHadoopCmdPath();
		fr_vars.hdfscmdpath = getHdfsCmdPath();
		fr_vars.streaming_path = getHadoopJarPath();
		fr_vars.binary_prefix = binpath + SLASH;
	} catch (exception& e) {
		cerr << "error here: " << e.what() << "\n";
		return false;
	} catch (...) {
		cerr << "Exception of unknown type!\n";
		return false;
	}

	return true;
}

bool join_data(struct framework_vars &fr_vars) {

	hdfs_delete(fr_vars.hadoopcmdpath, fr_vars.output_path);
	vector<string> arr_args = {"hadoop", "jar", fr_vars.streaming_path};
	arr_args.push_back("-files");
	stringstream ss;
	ss << fr_vars.binary_prefix + MANIPULATE << ","
	   << fr_vars.binary_prefix + RESQUE;
	arr_args.push_back(ss.str());

	arr_args.push_back("-input");
	arr_args.push_back(fr_vars.input_path_1);
	arr_args.push_back("-input");
	arr_args.push_back(fr_vars.input_path_2);

	arr_args.push_back("-output");
	arr_args.push_back(fr_vars.output_path);

	// the mapper phase assign each object to
	// different tiles with its mbb
	arr_args.push_back("-mapper");
	arr_args.push_back(MANIPULATE);

	// the resque tool received mbbs for spatial join
	arr_args.push_back("-reducer");
	ss.str("");
	ss <<RESQUE
	   <<" -l "<<fr_vars.decomp_lod
	   <<" -p "<<fr_vars.predicate;
	arr_args.push_back(ss.str()); // Offset to account for tile id and join index

	arr_args.push_back("-numReduceTasks");
	arr_args.push_back(to_string(fr_vars.numreducers));

	arr_args.push_back("-jobconf");
	arr_args.push_back("mapreduce.task.timeout=36000000");

	int status = 0;
	pid_t childpid;
	if ((childpid = execute_command(fr_vars.hadoopcmdpath, arr_args))) {
		if (wait(&status)) {
			cerr << "Succeeded in sp join: " << status << endl;
		} else {
			cerr << "Failed in sp join: " << status << endl;
			exit(1);
		}
		return status == 0 ? true : false;
	}
	return false;
}

int main(int argc, char** argv) {
	struct framework_vars fr_vars;
	if (!extract_params(argc, argv, fr_vars)) {
		return 1;
	}
	join_data(fr_vars);
}





