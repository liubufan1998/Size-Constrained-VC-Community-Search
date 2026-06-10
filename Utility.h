#ifndef _UTILITY_H_
#define _UTILITY_H_

#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iostream> 
#include <tuple>
#include <queue>
#include <set>
#include <cassert>

#define NDEBUG 



#define _RR1_ // control whether use RR1
#define _RR2_ // control whether use RR2		
#define _RR3_ // control whether use RR3		
#define _FAST_VCL_ // control whether use the proposed fast VCL ordering construction method

using ui = unsigned int;  
using ept = unsigned int; 

#define pb push_back
#define mp make_pair

#define mmax(a, b) ((a) > (b) ? (a) : (b))
#define mmin(a, b) ((a) < (b) ? (a) : (b))

class Utility
{
public:
	static FILE *open_file(const char *file_name, const char *mode)
	{
		FILE *f = fopen(file_name, mode);
		if (f == nullptr)
		{
			printf("Can not open file: %s\n", file_name);
			exit(1);
		}
		return f;
	}

	static std::string integer_to_string(long long number)
	{
		std::vector<ui> sequence;
		if (number == 0)
			sequence.pb(0);
		while (number > 0)
		{
			sequence.pb(number % 1000);
			number /= 1000;
		}
		char buf[5];
		std::string res;
		for (ui i = sequence.size(); i > 0; i--)
		{
			if (i == sequence.size())
				sprintf(buf, "%u", sequence[i - 1]);
			else
				sprintf(buf, ",%03u", sequence[i - 1]);
			res += std::string(buf);
		}
		return res;
	}
};

#endif
