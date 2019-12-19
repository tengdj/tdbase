#ifndef TOKENIZER_H
#define TOKENIZER_H

#include <string>
#include <vector>
#include <chrono>
using namespace std;

// name of binary tools
const std::string MANIPULATE = "manipulate_3d";
const std::string RESQUE = "resque_3d";

#define SPACE " "
#define SLASH "/"

enum Jointype{
	ST_ERROR = 0,
	ST_INTERSECTS = 1,
	ST_TOUCHES = 2,
	ST_CROSSES = 3,
	ST_CONTAINS = 4,
	ST_ADJACENT = 5,
	ST_DISJOINT = 6,
	ST_EQUALS = 7,
	ST_DWITHIN = 8,
	ST_WITHIN = 9,
	ST_OVERLAPS = 10,
	ST_NEAREST = 11,
	ST_NEAREST_2 = 12,
	ST_NN_VORONOI = 13,
	ST_NN_RTREE = 14
};

const std::string join_type_str[15] = {
		"st_error", "st_intersects", "st_touches", "st_crosses", "st_contains",
		"st_adjacent", "st_disjoint", "st_equals", "st_dwithin", "st_within",
		"st_overlaps", "st_nearest", "st_nearest2", "st_nn_voronoi", "st_nn_rtree"};

inline Jointype get_join_predicate(const char * predicate_str){
	for(int i=1;i<15;i++){
		if (strcmp(predicate_str, join_type_str[i].c_str()) == 0) {
			return (Jointype)i ;
		}
	}
	std::cerr << "unrecognized join predicate " <<predicate_str<< std::endl;
	return ST_ERROR;
}


class Timer
{
	public:
	Timer() : beg_(clock_::now()) {}
	void restart() { beg_ = clock_::now(); }

	double elapsed() const
	{
		return std::chrono::duration_cast<second_>(clock_::now() - beg_).count();
	}

	private:
	typedef std::chrono::high_resolution_clock clock_;
	typedef std::chrono::duration<double, std::ratio<1> > second_;
	std::chrono::time_point<clock_> beg_;
};


inline void tokenize ( const std::string& str, std::vector<std::string>& result,
	const std::string& delimiters = " ,;:\t", 
	const bool keepBlankFields=false,
	const std::string& quote="\"\'"
	)
{
    // clear the vector
    if ( false == result.empty() )
    {
	result.clear();
    }

    // you must be kidding
    if (delimiters.empty())
	return ;

    std::string::size_type pos = 0; // the current position (char) in the string
    char ch = 0; // buffer for the current character
    char delimiter = 0;	// the buffer for the delimiter char which
    // will be added to the tokens if the delimiter
    // is preserved
    char current_quote = 0; // the char of the current open quote
    bool quoted = false; // indicator if there is an open quote
    std::string token;  // string buffer for the token
    bool token_complete = false; // indicates if the current token is
    // read to be added to the result vector
    std::string::size_type len = str.length();  // length of the input-string

    // for every char in the input-string
    while ( len > pos )
    {
	// get the character of the string and reset the delimiter buffer
	ch = str.at(pos);
	delimiter = 0;

	bool add_char = true;

	// check ...

	// ... if the delimiter is a quote
	if ( false == quote.empty())
	{
	    // if quote chars are provided and the char isn't protected
	    if ( std::string::npos != quote.find_first_of(ch) )
	    {
		// if not quoted, set state to open quote and set
		// the quote character
		if ( false == quoted )
		{
		    quoted = true;
		    current_quote = ch;

		    // don't add the quote-char to the token
		    add_char = false;
		}
		else // if quote is open already
		{
		    // check if it is the matching character to close it
		    if ( current_quote == ch )
		    {
			// close quote and reset the quote character
			quoted = false;
			current_quote = 0;

			// don't add the quote-char to the token
			add_char = false;
		    }
		} // else
	    }
	}

	if ( false == delimiters.empty() && false == quoted )
	{
	    // if ch is delemiter 
	    if ( std::string::npos != delimiters.find_first_of(ch) )
	    {
		token_complete = true;
		// don't add the delimiter to the token
		add_char = false;
	    }
	}

	// add the character to the token
	if ( true == add_char )
	{
	    // add the current char
	    token.push_back( ch );
	}

	// add the token if it is complete
	// if ( true == token_complete && false == token.empty() )
	if ( true == token_complete )
	{
	    if (token.empty())
	    {
		if (keepBlankFields)
		    result.push_back("");
	    }
	    else 
		// add the token string
		result.push_back( token );

	    // clear the contents
	    token.clear();

	    // build the next token
	    token_complete = false;

	}
	// repeat for the next character
	++pos;
    } // while
    
    /* 
    cout << "ch: " << (int) ch << endl;
    cout << "token_complete: " << token_complete << endl;
    cout << "token: " << token<< endl;
    */
    // add the final token
    if ( false == token.empty() ) {
	result.push_back( token );
    }
    else if(keepBlankFields && std::string::npos != delimiters.find_first_of(ch) ){
	result.push_back("");
    }
}


inline void remove_slash(string &str){
	if (str.at(str.size() - 1) == '/') {
		str = str.substr(0, str.size() - 1);
	}
}

#endif

