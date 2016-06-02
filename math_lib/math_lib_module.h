#ifndef MATH_LIB_EXPORT_H
#define MATH_LIB_EXPORT_H

#ifdef MATH_LIB_STATIC_DEFINE
#  define MATH_LIB_EXPORT
#  define MATH_LIB_NO_EXPORT
#else
#  ifndef MATH_LIB_EXPORT
#    ifdef math_lib_EXPORTS
/* We are building this library */
#      define MATH_LIB_EXPORT __declspec(dllexport)
#    else
/* We are using this library */
#      define MATH_LIB_EXPORT __declspec(dllimport)
#    endif
#  endif

#  ifndef MATH_LIB_NO_EXPORT
#    define MATH_LIB_NO_EXPORT 
#  endif
#endif

#ifndef MATH_LIB_DEPRECATED
#  define MATH_LIB_DEPRECATED __declspec(deprecated)
#  define MATH_LIB_DEPRECATED_EXPORT MATH_LIB_EXPORT __declspec(deprecated)
#  define MATH_LIB_DEPRECATED_NO_EXPORT MATH_LIB_NO_EXPORT __declspec(deprecated)
#endif

#define DEFINE_NO_DEPRECATED 0
#if DEFINE_NO_DEPRECATED
# define MATH_LIB_NO_DEPRECATED
#endif



#endif