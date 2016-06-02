#ifndef SURFE_LIB_EXPORT_H
#define SURFE_LIB_EXPORT_H

#ifdef SURFE_LIB_STATIC_DEFINE
#  define SURFE_LIB_EXPORT
#  define SURFE_LIB_NO_EXPORT
#else
#  ifndef SURFE_LIB_EXPORT
#    ifdef surfe_lib_EXPORTS
/* We are building this library */
#      define SURFE_LIB_EXPORT __declspec(dllexport)
#    else
/* We are using this library */
#      define SURFE_LIB_EXPORT __declspec(dllimport)
#    endif
#  endif

#  ifndef SURFE_LIB_NO_EXPORT
#    define SURFE_LIB_NO_EXPORT 
#  endif
#endif

#ifndef SURFE_LIB_DEPRECATED
#  define SURFE_LIB_DEPRECATED __declspec(deprecated)
#  define SURFE_LIB_DEPRECATED_EXPORT SURFE_LIB_EXPORT __declspec(deprecated)
#  define SURFE_LIB_DEPRECATED_NO_EXPORT SURFE_LIB_NO_EXPORT __declspec(deprecated)
#endif

#define DEFINE_NO_DEPRECATED 0
#if DEFINE_NO_DEPRECATED
# define SURFE_LIB_NO_DEPRECATED
#endif



#endif