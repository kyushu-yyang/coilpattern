#####################################################################
## CMakeLists for compile the field quality code
##
## Author: Y.Yang
## Date  : 2020.03.07 (  create)
##         2020.07.17 (modified)
## Log   : modification 20170717
##         > example files are included
##         > delete the compile for main.cxx
##
#####################################################################

# define the cmake version
cmake_minimum_required( VERSION 2.8 )

# project name
project( CoilPattern )

# set the c++11
set(CMAKE_CXX_STANDARD, 14)

#
set(Boost_USE_STATIC_LIBS ON)
#
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
#
include_directories("package")

# setup
set( EXE_FILE CoilPattern.exe )
set( MAIN_FILE main.cxx )

# find the source code
file( GLOB_RECURSE EXAMPLES "examples/*cxx")
foreach( EACH_EXAMPLE IN LISTS EXAMPLES )
  message(STATUS "check example files: ${EACH_EXAMPLE}")
endforeach()

# include the header files
#add_subdirectory(src)

# for BOOST
find_package( Boost 1.50.0 REQUIRED )
message(STATUS "check Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}" )
message(STATUS "check Boost_LIBRARY_DIRS: ${Boost_LIBRARY_DIRS}" )

# for ROOT
#find_package( ROOT REQUIRED )
#message(STATUS "check ROOT_INCLUDE_DIRS: ${ROOT_INCLUDE_DIRS}" )
#message(STATUS "check ROOT_LIBRARIES   : ${ROOT_LIBRARIES}" )

# for Eigen
find_package( Eigen3 REQUIRED )
include_directories( ${PROJECT_SOURCE_DIR} ${EIGEN3_INCLUDE_DIR} )
message(STATUS "check EIGEN3_INCLUDE_DIR: ${EIGEN3_INCLUDE_DIR}")
#include_directories(/usr/include/eigen3)

# for OPENMP
find_package( OpenMP )
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )

# include the libraries head file
if (Boost_FOUND)
  include_directories( ${PROJECT_SOURCE_DIR}  ${Boost_INCLUDE_DIRS} )
else()
  message(FATAL_ERROR "boost is not found." )
endif()

#if (ROOT_FOUND)
#  include_directories( ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS} )
#else()
#  message(FATAL_ERROR "ROOT is not found.")
#endif()

# include the header files
add_subdirectory(src)

## create the executable file for examples
foreach( EACH_EXAMPLE IN LISTS EXAMPLES )
  ### extract the file name
  get_filename_component( EXEFILE_NAME ${EACH_EXAMPLE} NAME_WE )
  add_executable( ${EXEFILE_NAME}.exe ${EACH_EXAMPLE} )

  ### include the header files 
  target_link_libraries( ${EXEFILE_NAME}.exe CoilPattern )

  ### link the libraries (boost)
  if (Boost_FOUND)
    target_link_libraries( ${EXEFILE_NAME}.exe ${Boost_LIBRARIES} )
  endif()

  ### link the libraries (root)
  #if (ROOT_FOUND)
  #  target_link_libraries( ${EXEFILE_NAME}.exe ${ROOT_LIBRARIES} )
  #endif()

endforeach()

