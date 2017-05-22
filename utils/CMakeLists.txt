cmake_minimum_required(VERSION 2.8.9)
project(trajectoryanalysis)
#set(CMAKE_BUILD_TYPE Release)


# Set necessary flags
set(CMAKE_CXX_FLAGS		"-Wall -std=c++14")
set(CMAKE_CXX_FLAGS_RELEASE	"-O2 -DNDEBUG -funroll-loops -Wno-unknown-pragmas")
set(CMAKE_CXX_FLAGS_DEBUG	"-g")

#Set output directories
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

# set find paths
LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/Modules")


#FIND_PACKAGES for GSL
FIND_PACKAGE(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS} ${GSLCBLAS_INCLUDE_DIRS})
list(APPEND LIBS ${GSL_LIB} ${GSL_CBLAS_LIB})
#GSL include and header files
#include_directories(/usr/include)
#link_directories(/usr/lib64)

#Add defintions
add_definitions(-DIO_XDR)

#XDR stuff
#to read xtc files
link_directories(/home/flatinwo/software/local/lib)
include_directories(/home/flatinwo/software/local/include/xdrfile)

#include the headers
include_directories(src/trajectoryanalysis/src)
include_directories(src/trajectoryanalysis/chirality)
include_directories(src/trajectoryanalysis/water)

#add all source files
file(GLOB SOURCES "src/trajectoryanalysis/src/*.cpp"
		  "src/trajectoryanalysis/chirality/*.cpp"
		  "src/trajectoryanalysis/water/*.cpp")


#add executable
file(GLOB BCHRLTY "src/trajectoryanalysis/drivers/*/*.cpp")

#Generate static library
add_library(trajanalysis STATIC ${SOURCES})


#Generate executables

message(STATUS "THESE ARE THE EXTERNAL LIBRARIES INCLUDED: " ${LIBS})

add_executable(chiralityBOP ${BCHRLTY} ${SOURCES})
target_link_libraries(chiralityBOP ${LIBS} xdrfile)
