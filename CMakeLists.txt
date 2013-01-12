# CMake version control ------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8)

# Core Project declaration ---------------------------------------------------------------
project(LandySzalayEstimator)
set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})

# TBB requirements -----------------------------------------------------------------------
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules/")
find_package(TBB REQUIRED)
include_directories(${TBB_INCLUDE_DIRS})
link_directories(${TBB_LIBRARY_DIRS})

# BOOST REQUIREMENTS ---------------------------------------------------------------------
find_package(Boost REQUIRED)
if (Boost_FOUND)
  message("-- Found Boost")
else ()
  message(FATAL_ERROR "-- Unable to Find Boost")
endif()
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS})

# EIGEN REQUIREMENTS ---------------------------------------------------------------------
find_package(Eigen3 REQUIRED)
if (EIGEN3_FOUND)
  message("-- Found Eigen3")
else ()
  message(FATAL_ERROR "-- Unable to find Eigen3")
endif()
include_directories(${EIGEN3_INCLUDE_DIR})

# C++11 ----------------------------------------------------------------------------------
add_definitions(-std=c++11)

# Set core files list --------------------------------------------------------------------
file(
	GLOB
	core_source_files
	sources/*
	GLOB
	core_include_files
	includes/*
)

# add executable for the core BLINK system
add_executable(
	LandySzalayEstimator
	${core_include_files}
	${core_source_files}
)

target_link_libraries(LandySzalayEstimator ${TBB_LIBRARIES})
