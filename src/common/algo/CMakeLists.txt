cmake_minimum_required(VERSION 3.20)
project(algo)

# Add library after setting up PCL
add_library(${PROJECT_NAME} STATIC
    registration.hpp
    registration.cpp
    # visualization.hpp
    # visualization.cpp
)

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# # Json
# target_link_libraries(${PROJECT_NAME} PRIVATE nlohmann_json::nlohmann_json)

# OpenCV
find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PRIVATE igl::core ${OpenCV_LIBS})
