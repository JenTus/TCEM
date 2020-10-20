#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <unistd.h>
#include <sys/types.h>
#include <sstream>
#include <climits>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>

namespace _Cide{
	extern std::vector< std::vector<int> > graphT;
}


typedef long long int64;
typedef unsigned long long uint64;

#ifdef _WIN32
#define OS_SEP '\\'
#else
#define OS_SEP '/'
#endif

using namespace std;

float getCurrentMemoryUsage();
float getRunningTime(time_t startTime);

void stringTokenizer(string& str, float *tokens, int size, string& delimiters);
void stringTokenizer(string& str, double *tokens, int size, string& delimiters);

void rtrim(char *str);
void ltrim(char *str);
void trim(char *str);

inline unsigned int strToInt(string s) {
  unsigned int i;
  istringstream myStream(s);

  if(myStream >> i) 
    return i;
  
  else {
    std::cout << "String " << s << " is not a number." << endl;
    return atoi(s.c_str());
  }
  
  return i;
}

inline unsigned int strToUInt(const char* s) {

    unsigned int i;
    istringstream myStream(s);
    
    if (myStream >> i)
        return i;
    
    else {
        cout << "String " << s << " is not a number." << endl;
        return stoi(s);
    }
    
    return i;
}


inline int strToInt(const char* s) {

  int i;
  istringstream myStream(s);

  if (myStream >> i)
    return i;
  
  else {
    cout << "String " << s << " is not a number." << endl;
    return atoi(s);
    }
  
  return i;
}

inline float strToFloat(const char* s) {
  return atof(s);  
}

inline float strToFloat(string s) {
  return atof(s.c_str());
}

inline string floatToStr(float f) {
  stringstream ss;
  ss << f;
  return ss.str();
}

inline int64_t strToInt64(string s) {
  int64_t i;
  istringstream myStream(s);

  if (myStream >> i)
    return i;
  
  else {
    cout << "String " << s << " is not a number." << endl;
    exit(1);    
  }
  
  return i;
}

inline string intToStr(int i) {
  stringstream ss;
  ss << i;
  return ss.str();  
}

inline double strToDouble(string s) {	
// 	return std::stod(s);
	double a = 0; 
	stringstream ss;
	ss << s;
	ss >> a;
	return a;
}

inline double logcnk(int n, int k) {
    double ans = 0;
    for (int i = n - k + 1; i <= n; i++)
    {
        ans += std::log(i);
    }
    for (int i = 1; i <= k; i++)
    {
        ans -= std::log(i);
    }
    return ans;
}

//
inline double logP(int n, int k_r, int tau) {
    double ans = 0;
    
    // computes log of perm(n, k_r * (tau+1))
    for (int i = n - k_r *(tau+1) + 1; i <= n; i++)
    {
        ans += std::log(i);
    }
    
    // computes log of k_r!
    for (int i = 1; i <= k_r; i++)
    {
        ans -= std::log(i);
    }
    
    // computes log of tau!
    double t = 0;
    for (int i = 1; i <= tau; i++)
    {
        t += std::log(i);
    }
    ans -= k_r * t;
    return ans;
}


//
// add some calculations
inline double fact(int n){
    if(n == 1)
        return 1;
    else
        return n*fact(n-1);
}

//
inline double Permutation(int n, int k){
    return (double) fact(n)/ (double) fact(n - k);
//    return (double) fact(n)/ (double) fact(k);
}

inline double Combination(int n, int k){
    return Permutation(n, k) / fact(k);
}

inline double sqr(double t)
{
    return t * t;
}

template <class Container>
void split1(const std::string& str, Container& cont)
{
    std::istringstream iss(str);
    std::copy(std::istream_iterator<std::string>(iss),
              std::istream_iterator<std::string>(),
              std::back_inserter(cont));
}

// commented out below - using from the c++ math library instead
//static double log2(int n){
//    return std::log(n) / std::log(2);
//}

std::vector<int> intersection(std::vector<int> &v1, std::vector<int> &v2);

struct Counter
{
    struct value_type { template<typename T> value_type(const T&) { } };
    void push_back(const value_type&) { ++count; }
    size_t count = 0;
};

template<typename T1, typename T2>
size_t intersection_size(const T1& s1, const T2& s2)
{
    Counter c;
    set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter(c));
    return c.count;
}


#endif
