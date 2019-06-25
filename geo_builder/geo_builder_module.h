#ifndef GEO_BUILDER_EXPORT_H
#define GEO_BUILDER_EXPORT_H
#ifdef _WIN32
#ifdef GEO_BUILDER_STATIC_DEFINE
#define GEO_BUILDER_EXPORT
#define GEO_BUILDER_NO_EXPORT
#else
#ifndef GEO_BUILDER_EXPORT
#ifdef math_lib_EXPORTS
/* We are building this library */
#define GEO_BUILDER_EXPORT __declspec(dllexport)
#else
/* We are using this library */
#define GEO_BUILDER_EXPORT __declspec(dllimport)
#endif
#endif

#ifndef GEO_BUILDER_NO_EXPORT
#define GEO_BUILDER_NO_EXPORT
#endif
#endif

#ifndef GEO_BUILDER_DEPRECATED
#define GEO_BUILDER_DEPRECATED __declspec(deprecated)
#define GEO_BUILDER_DEPRECATED_EXPORT GEO_BUILDER_EXPORT __declspec(deprecated)
#define GEO_BUILDER_DEPRECATED_NO_EXPORT GEO_BUILDER_NO_EXPORT __declspec(deprecated)
#endif

#define DEFINE_NO_DEPRECATED 0
#if DEFINE_NO_DEPRECATED
#define GEO_BUILDER_NO_DEPRECATED
#endif
#endif
#ifndef _WIN32
#define GEO_BUILDER_EXPORT
#endif
#endif
