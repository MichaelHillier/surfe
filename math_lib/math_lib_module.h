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

#ifndef MATH_LIB_EXPORT_H
#define MATH_LIB_EXPORT_H

#ifdef MATH_LIB_STATIC_DEFINE
#define MATH_LIB_EXPORT
#define MATH_LIB_NO_EXPORT
#else
#ifndef MATH_LIB_EXPORT
#ifdef math_lib_EXPORTS
/* We are building this library */
#define MATH_LIB_EXPORT __declspec(dllexport)
#else
/* We are using this library */
#define MATH_LIB_EXPORT __declspec(dllimport)
#endif
#endif

#ifndef MATH_LIB_NO_EXPORT
#define MATH_LIB_NO_EXPORT
#endif
#endif

#ifndef MATH_LIB_DEPRECATED
#define MATH_LIB_DEPRECATED __declspec(deprecated)
#define MATH_LIB_DEPRECATED_EXPORT MATH_LIB_EXPORT __declspec(deprecated)
#define MATH_LIB_DEPRECATED_NO_EXPORT MATH_LIB_NO_EXPORT __declspec(deprecated)
#endif

#define DEFINE_NO_DEPRECATED 0
#if DEFINE_NO_DEPRECATED
#define MATH_LIB_NO_DEPRECATED
#endif

#endif