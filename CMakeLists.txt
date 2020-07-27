cmake_minimum_required(VERSION 2.8)
 
PROJECT(surfe)
#SET(BUILD_SHARED_LIBS_ON)
SET(CMAKE_CXX_STANDARD 11 )

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
#include_directories(${PROJECT_SOURCE_DIR}/eigen-git-mirror/Eigen)#add_subdirectory(eigen-git-mirror)

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
target_link_libraries(surfe_lib math_lib)

add_subdirectory(pybind11)
#setup python
pybind11_add_module(surfepy surfe_pybindings/pybindings.cpp)
target_link_libraries(surfepy surfe_lib math_lib)
install(TARGETS surfepy surfe_lib math_lib DESTINATION surfepy)
option(GEO_BUILDER FALSE)
if (GEO_BUILDER)
	find_package(VTK REQUIRED)
	if (NOT VTK_FOUND)
		MESSAGE(FATAL_ERROR "Could not find VTK")
	endif()
	MESSAGE(STATUS "VTK VERSION is: " ${VTK_VERSION})
	if (${VTK_VERSION} VERSION_LESS_EQUAL 8.2)
		include(${VTK_USE_FILE})
	endif()
	# QT library
	find_package(Qt5 COMPONENTS REQUIRED
		Core
		Gui
		Widgets
	)
	# To find the qt ui files , implementation files, and generated qt files (moc, h)
	set(Qt_FILES ${CMAKE_CURRENT_SOURCE_DIR}/geo_builder/ui)
	include_directories(${Qt_FILES} ${CMAKE_CURRENT_BINARY_DIR})
	qt5_wrap_ui(ui_wrap ${Qt_FILES}/input.ui)
	qt5_wrap_cpp(moc_sources ${Qt_FILES}/inputImpl.h)
	qt5_add_resources(rcc_resources ${Qt_FILES}/resources.qrc)

	FILE(GLOB SRC_UI 
		"${Qt_FILES}/*.h"
		"${Qt_FILES}/*.cpp"
		"${Qt_FILES}/*.qrc"
	)
	source_group("UI Files"	FILES ${SRC_UI})
	
	FILE(GLOB GEO_BUILDER_HEADERS "geo_builder/*.h")
	FILE(GLOB GEO_BUILDER_SOURCES "geo_builder/*.cpp")

	add_library(geo_builder SHARED ${GEO_BUILDER_HEADERS} ${GEO_BUILDER_SOURCES} ${SRC_UI} ${ui_wrap} ${moc_sources} ${rcc_resources})
	target_include_directories(geo_builder PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/geo_builder
											${CMAKE_CURRENT_SOURCE_DIR}/surfe_lib
                                            ${CMAKE_CURRENT_SOURCE_DIR}/math_lib)
	target_link_libraries(geo_builder surfe_lib math_lib ${VTK_LIBRARIES} Qt5::Core Qt5::Gui Qt5::Widgets)
	
	#Setup test
	FILE(GLOB TEST_HEADERS "test/*.h")
	FILE(GLOB TEST_SOURCES "test/*.cpp")

	ADD_EXECUTABLE(test test/main.cpp)
	target_link_libraries(test geo_builder surfe_lib math_lib)
	include_directories(${CMAKE_CURRENT_SOURCE_DIR}/test
		${CMAKE_CURRENT_SOURCE_DIR}/geo_builder
		${CMAKE_CURRENT_SOURCE_DIR}/surfe_lib
		${CMAKE_CURRENT_SOURCE_DIR}/math_lib)
	
endif()

