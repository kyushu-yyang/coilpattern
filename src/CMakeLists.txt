# -------------------------------
# cmakelists for the source file
# -------------------------------

# option
option( COILPATTERN_BUILD_SHARED_LIBS "build supermadxx as a shared library" ON )

# set the name of this library
set( LIBRARY_NAME CoilPattern )

# find the source files
file( GLOB_RECURSE SOURCE_FILES "*cxx" )

## EIGEN
include_directories(/usr/include/eigen3)

## ROOT
#find_package( ROOT REQUIRED )
#include_directories( ${ROOT_INCLUDE_DIRS} )

# create a library
if (COILPATTERN_BUILD_SHARED_LIBS)
  add_library( ${LIBRARY_NAME} SHARED ${SOURCE_FILES} )
else()
  add_library( ${LIBRARY_NAME} STATIC ${SOURCE_FILES} )
endif()

# include the head files
target_include_directories( ${LIBRARY_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/include )

