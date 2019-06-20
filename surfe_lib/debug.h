// SURFace Estimator(SURFE) - Terms and Conditions of Use

// Unless otherwise noted, computer program source code of the SURFace
// Estimator(SURFE) is covered under Crown Copyright, Government of Canada, and
// is distributed under the MIT License.

// The Canada wordmark and related graphics associated with this distribution
// are protected under trademark law and copyright law.No permission is granted
// to use them outside the parameters of the Government of Canada's corporate
// identity program. For more information, see
// http://www.tbs-sct.gc.ca/fip-pcim/index-eng.asp

// Copyright title to all 3rd party software distributed with the SURFace
// Estimator(SURFE) is held by the respective copyright holders as noted in
// those files.Users are asked to read the 3rd Party Licenses referenced with
// those assets.

// MIT License

// Copyright(c) 2017 Government of Canada

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files(the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and / or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions :

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT
// NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
// PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef debug_h
#define debug_h

/*#include <gmpxx.h>*/
#include <modeling_methods.h>
#include <modelling_input.h>
#include <typeinfo>
// if on windows
#ifdef _WIN32
#include <Windows.h>
#endif
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

//typedef mpf_class mpfc;

void open_console_window();

// template <class T>
// void outc(const std::vector<T> &v, const std::string &description);
// template <class T>
// void outc(const T &s, const std::string &description);
// template <class T>
// void outc(const std::vector<std::vector<T> > &m,
// 	const std::string &description);
// 
// template <class T>
// void outc(const std::vector<T> &v, const std::string &description) {
// #ifdef _WIN32
// 	if (GetConsoleWindow() == NULL) open_console_window();
// #endif
// 	std::cout << description << std::endl;
// 	typedef T type;
// 	for (int j = 0; j < (int)v.size(); j++) {
// 		if (typeid(type) == typeid(mpfc))
// 			std::cout << "v[" << j << "] = " << v[j].get_d() << std::endl;
// 		else
// 			std::cout << "v[" << j << "] = " << v[j] << std::endl;
// 	}
// }
// 
// template <class T>
// void outc(const T &s, const std::string &description) {
// #ifdef _WIN32
// 	if (GetConsoleWindow() == NULL) open_console_window();
// #endif
// 	std::cout << description;
// 	typedef T type;
// 	if (typeid(type) == typeid(mpfc))
// 		std::cout << " = " << s.get_d() << std::endl;
// 	else
// 		std::cout << " = " << s << std::endl;
// }
// 
// template <class T>
// void outc(const std::vector<std::vector<T> > &m,
// 	const std::string &description) {
// #ifdef _WIN32
// 	if (GetConsoleWindow() == NULL) open_console_window();
// #endif
// 	std::cout << description << std::endl;
// 	typedef T type;
// 	for (int j = 0; j < (int)m.size(); j++) {
// 		for (int k = 0; k < (int)m[j].size(); k++) {
// 			if (typeid(type) == typeid(mpfc))
// 				std::cout << setw(15) << m[j].get_d();
// 			else
// 				std::cout << setw(15) << m[j];
// 		}
// 		std::cout << std::endl;
// 	}
// }

#endif
