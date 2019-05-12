cmake_minimum_required(VERSION 2.8)
 
PROJECT(surfe)
#SET(BUILD_SHARED_LIBS_ON)
SET(CMAKE_CXX_STANDARD 11 )
# Get GMP dependency
find_path(GMP_INCLUDE_DIR gmp.h gmpxx.h)
find_library(GMP_LIBRARY_C NAMES gmp)
find_library(GMP_LIBRARY_CPP NAMES gmpxx)

if(GMP_INCLUDE_DIR AND GMP_LIBRARY_C AND GMP_LIBRARY_CPP)
    get_filename_component(GMP_LIBRARY_C_DIR ${GMP_LIBRARY_C} PATH)
	get_filename_component(GMP_LIBRARY_CPP_DIR ${GMP_LIBRARY_CPP} PATH)
    set(GMP_FOUND TRUE)
endif()

if(GMP_FOUND)
   MESSAGE(STATUS "Found GMP_C: ${GMP_LIBRARY_C}")
   MESSAGE(STATUS "Found GMP_CPP: ${GMP_LIBRARY_CPP}")
else()
   MESSAGE(FATAL_ERROR "Could not find GMP")
endif()
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Get Eigen dependency
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# Set link directories
link_directories(${GMP_LIBRARY_C_DIR} ${GMP_LIBRARY_CPP_DIR})

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
#Setup test
FILE(GLOB TEST_HEADERS "test/*.h")
FILE(GLOB TEST_SOURCES "test/*.cpp")

#add_library(test SHARED ${TEST_HEADERS} ${TEST_SOURCES})

include_directories(${GMP_INCLUDE_DIR})

target_link_libraries(math_lib ${GMP_LIBRARY_C} ${GMP_LIBRARY_CPP})
target_link_libraries(surfe_lib math_lib ${GMP_LIBRARY_C} ${GMP_LIBRARY_CPP})
ADD_EXECUTABLE(test_func test/test.cpp)
target_link_libraries(test_func surfe_lib math_lib ${GMP_LIBRARY_C} ${GMP_LIBRARY_CPP})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/test ${CMAKE_CURRENT_SOURCE_DIR}/surfe_lib
                                            ${CMAKE_CURRENT_SOURCE_DIR}/math_lib)

add_subdirectory(pybind11)
#setup python
pybind11_add_module(surfepy surfe_pybindings/pybindings.cpp)
target_link_libraries(surfepy surfe_lib math_lib ${GMP_LIBRARY_C} ${GMP_LIBRARY_CPP})

