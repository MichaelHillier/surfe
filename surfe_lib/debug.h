#ifndef debug_h
#define debug_h

#include <mpirxx.h>
#include <modelling_input.h>
#include <modeling_methods.h>

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <Windows.h>

typedef mpf_class mpfc;

void open_console_window();

template <class T> void outc(const std::vector< T > &v, const std::string &description);
template <class T> void outc(const T &s, const std::string &description);
template <class T> void outc(const std::vector < std::vector < T > > &m, const std::string &description);

template <class T>
void outc(const std::vector< T > &v)
{
	if (GetConsoleWindow() == NULL) open_console_window();
	std::cout << description << std::endl;
	typedef T type;
	for (int j = 0; j < (int)v.size(); j++){
		if (type == mpfc) std::cout << "v[" << j << "] = " << v[j].get_d() << std::endl;
		else std::cout << "v[" << j << "] = " << v[j] << std::endl;
	}
}

template <class T>
void outc(const T &s, const std::string &description)
{
	if (GetConsoleWindow() == NULL) open_console_window();
	std::cout << description;
	typedef T type;
	if (type == mpfc) std::cout << " = " << s.get_d() << std::endl;
	else std::cout << " = " << s << std::endl;
}

template <class T>
void outc(const std::vector < std::vector < T > > &m, const std::string &description)
{
	if (GetConsoleWindow() == NULL) open_console_window();
	std::cout << description << std::endl;
	typedef T type;
	for (int j = 0; j < (int)m.size(); j++){
		for (int k = 0; k < (int)m[j].size(); k++) {
			if (type == mpfc) std::cout << setw(15) << v[j].get_d();
			else std::cout << setw(15) << v[j];
		}
		std::cout << std::endl;
	}
}

#endif