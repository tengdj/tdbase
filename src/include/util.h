#ifndef TENG_UTIL_H_
#define TENG_UTIL_H_

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <sstream>
#include <vector>
#include <thread>
#include <iostream>
#include <stdarg.h>
#include <time.h>
#include <string>
#include <fstream>
#include <streambuf>
#include <random>
#include <sys/stat.h>
#include <sys/types.h>

#if __linux
#include <sys/syscall.h>
#include <sys/time.h>
#include <unistd.h>
#include <pthread.h>
#include <dirent.h>
#include <cstdint>
#elif defined(_WIN32) || defined(_WIN64)
#include <windows.h>       // Or something like it. 
#include <filesystem>
namespace fs = std::filesystem;

#endif

using namespace std;

namespace tdbase{

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
 * timer functions
 * */
#define TENG_RANDOM_NUMBER 0315
#define INIT_TIME struct timeval start = get_cur_time();

#if defined(_WIN32) || defined(_WIN64)
inline int gettimeofday(struct timeval* tp, struct timezone* tzp)
{
	// Note: some broken versions only have 8 trailing zero's, the correct epoch has 9 trailing zero's
	// This magic number is the number of 100 nanosecond intervals since January 1, 1601 (UTC)
	// until 00:00:00 January 1, 1970 
	static const uint64_t EPOCH = ((uint64_t)116444736000000000ULL);

	SYSTEMTIME  system_time;
	FILETIME    file_time;
	uint64_t    time;

	GetSystemTime(&system_time);
	SystemTimeToFileTime(&system_time, &file_time);
	time = ((uint64_t)file_time.dwLowDateTime);
	time += ((uint64_t)file_time.dwHighDateTime) << 32;

	tp->tv_sec = (long)((time - EPOCH) / 10000000L);
	tp->tv_usec = (long)(system_time.wMilliseconds * 1000);
	return 0;
}
#endif


inline struct timeval get_cur_time(){
	struct timeval t1;
	gettimeofday(&t1, NULL);
	return t1;
}
inline double get_time_elapsed(struct timeval &t1, bool update_start = false){
	struct timeval t2;
    double elapsedTime;
	gettimeofday(&t2, NULL);
	// compute and print the elapsed time in millisec
	elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
	elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
	if(update_start){
		t1 = get_cur_time();
	}
	return elapsedTime;
}

inline string time_string(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	struct tm *nowtm;
	char tmbuf[100];
	char buf[256];
#if defined(_WIN32) || defined(_WIN64)
	long long tvsec = tv.tv_sec;
	nowtm = localtime(&tvsec);
#else
	nowtm = localtime(&tv.tv_sec);
#endif

	strftime(tmbuf, sizeof tmbuf, "%H:%M:%S", nowtm);
	sprintf(buf,"%s.%04ld", tmbuf, tv.tv_usec/1000);
	return string(buf);
}

class Timer{
	struct timeval t;
public:
	Timer(){
		t = get_cur_time();
	}
	double time_elapsed(bool update_start = false){
		return get_time_elapsed(t, update_start);
	}
	void reset(){
		t = get_cur_time();
	}
};

/*
 * log functions
 * */

static pthread_mutex_t print_lock;
inline double logt(const char *format, struct timeval &start, ...){

	pthread_mutex_lock(&print_lock);
	va_list args;
	va_start(args, start);
	char sprint_buf[200];
	int n = vsprintf(sprint_buf, format, args);
	va_end(args);
	long tid = 0;
#if __linux
	tid = syscall(__NR_gettid);
#elif defined(_WIN32) || defined(_WIN64)
	tid = GetCurrentThreadId();
#endif

	fprintf(stderr,"%s thread %ld:\t%s", time_string().c_str(), tid,sprint_buf);

	double mstime = get_time_elapsed(start, true);
	if(mstime>1000){
		fprintf(stderr," takes %f s\n", mstime/1000);
	}else{
		fprintf(stderr," takes %f ms\n", mstime);
	}
	fflush(stderr);

	pthread_mutex_unlock(&print_lock);
	return mstime;
}

inline void log(const char *format, ...){
	pthread_mutex_lock(&print_lock);
	va_list args;
	va_start(args, format);
	char sprint_buf[200];
	int n = vsprintf(sprint_buf, format, args);
	va_end(args);
	long tid = 0;
#if __linux
	tid = syscall(__NR_gettid);
#elif defined(_WIN32) || defined(_WIN64)
	tid = GetCurrentThreadId();
#endif
	fprintf(stderr,"%s thread %ld:\t%s\n", time_string().c_str(), tid,sprint_buf);
	fflush(stderr);
	pthread_mutex_unlock(&print_lock);
}


inline int get_num_threads(){
	return std::thread::hardware_concurrency();
}

static pthread_mutex_t plock;
inline void process_lock(){
	pthread_mutex_lock(&plock);
}

inline void process_unlock(){
	pthread_mutex_unlock(&plock);
}

/*
 * some random number based functions
 *
 * */

// get an integer between 0 and max_value
inline int get_rand_number(int max_value){
	std::random_device dev;
	std::mt19937 rng(dev());
	std::uniform_int_distribution<std::mt19937::result_type> dist(0, max_value);
	return dist(rng);
}

// get a double number between 0 and 1
inline double get_rand_double(){
	return get_rand_number(RAND_MAX)*1.0/RAND_MAX;
}

// try luck with a given possibility
inline bool tryluck(float possibility){
	assert(possibility>=0);
	return possibility>=1.0||get_rand_double()<=possibility;
}

inline bool flip_coin(){
	return tryluck(0.5);
}

/*
 * file or folder operations
 * */
inline bool is_dir(const char* path) {
    struct stat buf;
    stat(path, &buf);
    return S_ISDIR(buf.st_mode);
}

inline bool is_file(const char* path) {
    struct stat buf;
    stat(path, &buf);
    return S_ISREG(buf.st_mode);
}

inline void list_files(const char *path, std::vector<string> &f_list){
	if(is_file(path)){
		f_list.push_back(std::string(path));
		return;
	}
	DIR *dir;
	struct dirent *ent;
	if ((dir = opendir (path)) != NULL) {
		/* print all the files and directories within directory */
		while ((ent = readdir (dir)) != NULL) {
			if(strcmp(ent->d_name,"..")==0||
					strcmp(ent->d_name,".")==0){
				continue;
			}
			stringstream ss;
			ss<<path<<"/"<<ent->d_name;
			list_files(ss.str().c_str(), f_list);
		}
		closedir (dir);
	}
}


inline long file_size(const char *file){
	struct stat stat_buf;
	int rc = stat(file, &stat_buf);
	return rc == 0 ? stat_buf.st_size : -1;
}

inline long file_size(std::vector<string> &f_list){
	long size = 0;
	for(string s:f_list){
		long ls = file_size(s.c_str());
		if(ls>0){
			size += ls;
		}
	}
	return size;
}

inline bool file_exist(const char *path) {
  struct stat buffer;
  return (stat(path, &buffer) == 0);
}

inline string read_line(){
	string input_line;
	getline(std::cin, input_line);
	return input_line;
}

inline string read_file(const char *path){
	std::ifstream t(path);
	std::string str((std::istreambuf_iterator<char>(t)),
					 std::istreambuf_iterator<char>());
	t.close();
	return str;
}

inline void write_file(const string &content, const char *path){
	std::filebuf fb;
	fb.open(path, std::ios::out | std::ios::trunc);
	if(fb.is_open())
	{
		std::ostream os(&fb);
		os << content.c_str();
	}else{
		std::cerr<<"cannot find path "<<path<<std::endl;
	}
}


inline void tokenize( const std::string& str, std::vector<std::string>& result,
	const std::string& delimiters = " ,;:\t",
	const bool keepBlankFields=false,
	const std::string& quote="\"\'"
	){
    // clear the vector
    if (!result.empty()){
    	result.clear();
    }

    // you must be kidding
    if (delimiters.empty())
	return ;

    std::string::size_type pos = 0; // the current position (char) in the string
    char ch = 0; // buffer for the current character

    char current_quote = 0; // the char of the current open quote
    bool quoted = false; // indicator if there is an open quote
    std::string token;  // string buffer for the token
    bool token_complete = false; // indicates if the current token is
    // read to be added to the result vector
    std::string::size_type len = str.length();  // length of the input-string

    // for every char in the input-string
	while(len > pos){
		// get the character of the string and reset the delimiter buffer
		ch = str.at(pos);

		bool add_char = true;
		if ( false == quote.empty()){
			// if quote chars are provided and the char isn't protected
			if (std::string::npos != quote.find_first_of(ch)){
				if (!quoted){
					quoted = true;
					current_quote = ch;
					add_char = false;
				} else {
					if (current_quote == ch){
						quoted = false;
						current_quote = 0;
						add_char = false;
					}
				}
			}
		}

		if (!delimiters.empty()&&!quoted){
			// if ch is delemiter
			if (std::string::npos != delimiters.find_first_of(ch)){
				token_complete = true;
				// don't add the delimiter to the token
				add_char = false;
			}
		}

		// add the character to the token
		if (add_char){
			token.push_back(ch);
		}

		// add the token if it is complete
		// if ( true == token_complete && false == token.empty() )
		if (token_complete){
			if (token.empty())
			{
			if (keepBlankFields)
				result.push_back("");
			}
			else
			result.push_back( token );
			token.clear();
			token_complete = false;
		}
		++pos;
    } // while
    // add the final token
    if ( false == token.empty() ) {
    	result.push_back( token );
    } else if(keepBlankFields && std::string::npos != delimiters.find_first_of(ch) ){
    	result.push_back("");
    }
}


inline void remove_slash(string &str){
	if(str.at(str.size() - 1) == '/'){
		str = str.substr(0, str.size() - 1);
	}
}

}
#endif
