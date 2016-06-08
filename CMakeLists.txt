cmake_minimum_required(VERSION 2.8)
 
PROJECT(surfe)

# Get MPIR dependency
find_path(MPIR_INCLUDE_DIR mpir.h mpirxx.h)
find_library(MPIR_LIBRARY_C NAMES mpir)
find_library(MPIR_LIBRARY_CPP NAMES mpirxx)

if(MPIR_INCLUDE_DIR AND MPIR_LIBRARY_C AND MPIR_LIBRARY_CPP)
    get_filename_component(MPIR_LIBRARY_C_DIR ${MPIR_LIBRARY_C} PATH)
	get_filename_component(MPIR_LIBRARY_CPP_DIR ${MPIR_LIBRARY_CPP} PATH)
    set(MPIR_FOUND TRUE)
endif()

if(MPIR_FOUND)
   MESSAGE(STATUS "Found MPIR_C: ${MPIR_LIBRARY_C}")
   MESSAGE(STATUS "Found MPIR_CPP: ${MPIR_LIBRARY_CPP}")
else()
   MESSAGE(FATAL_ERROR "Could not find MPIR")
endif()

find_package(OpenMP)
if (OpenMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()	

# Set link directories
link_directories(${MPIR_LIBRARY_C_DIR} ${MPIR_LIBRARY_CPP_DIR})

# Setup math_lib
FILE(GLOB MATH_LIB_HEADERS "math_lib/*.h")
FILE(GLOB MATH_LIB_SOURCES "math_lib/*.cpp")

add_library(math_lib SHARED ${MATH_LIB_HEADERS} ${MATH_LIB_SOURCES})
target_include_directories(math_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/math_lib)

#Setup surfe_lib
FILE(GLOB SURFE_LIB_HEADERS "surfe_lib/*.h")
FILE(GLOB SURFE_LIB_SOURCES "surfe_lib/*.cpp")

add_library(surfe_lib SHARED ${SURFE_LIB_HEADERS} ${SURFE_LIB_SOURCES})
target_include_directories(surfe_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/surfe_lib
                                            ${CMAKE_CURRENT_SOURCE_DIR}/math_lib)

include_directories(${MPIR_INCLUDE_DIR})

target_link_libraries(math_lib ${MPIR_LIBRARY_C} ${MPIR_LIBRARY_CPP})
target_link_libraries(surfe_lib ${MPIR_LIBRARY_C} ${MPIR_LIBRARY_CPP})