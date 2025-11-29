cmake_minimum_required(VERSION 3.22.1)


# ######## MEL #########
set(PROJECT_NAME abstractSquidgetTool)
project(${PROJECT_NAME})

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(MSVC)
elseif(APPLE)
  add_definitions(-Wno-shorten-64-to-32)
  add_definitions(-Wsign-compare)
endif()

set(MAYA_WANT_CPP_17 YES)

include(FetchContent)

# ######## Maya #########
# Requires MAYA devkit
include($ENV{DEVKIT_LOCATION}/cmake/pluginEntry.cmake)
message("DEVKIT_LOCATION: $ENV{DEVKIT_LOCATION}")

# ######## Eigen3 #########
message(STATUS "EIGEN_INCLUDE: $ENV{EIGEN_INCLUDE}")
include_directories($ENV{EIGEN_INCLUDE})


# ######## OpenCV #########
# Used for registration function
find_package(OpenCV REQUIRED COMPONENTS core flann imgproc surface_matching)
message(STATUS "OpenCV_INCLUDE_DIRS: ${OpenCV_INCLUDE_DIRS}")
set(PACKAGE_LIBS ${PACKAGE_LIBS} ${OpenCV_LIBS})
message(STATUS "PACKAGE_LIBS ${PACKAGE_LIBS} ${OpenCV_LIBS}")
set(INCLUDE_DIRS ${INCLUDE_DIRS} ${OpenCV_INCLUDE_DIRS})

# ######## Libigl #########
# Used for deformation
include(libigl)
include_directories(${libigl_SOURCE_DIR}/include)

# ######## Workspace #########
set(WORKSPACE_INCLUDE "${CMAKE_SOURCE_DIR}/src")
include_directories(${WORKSPACE_INCLUDE})
message(STATUS "Include src ${WORKSPACE_INCLUDE}")

# ######## Source #########
file(GLOB_RECURSE SRC_CPP "${CMAKE_SOURCE_DIR}/src/*.cpp")

message(STATUS ${SRC_CPP})
set(SOURCE_FILES pluginMain.cpp ${SRC_CPP})

# ${MEL_FILES})
# set(MEL_FILES
#     "./mel/abstractToolProperties.mel"
#     "./mel/abstractToolValues.mel")

set(LIBRARIES OpenMaya OpenMayaAnim OpenMayaUI OpenMayaRender OpenMayaFX
  Foundation)

find_opengl()
build_plugin()
