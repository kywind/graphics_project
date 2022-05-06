CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(FinalProject CXX)

IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release)
ENDIF()

FIND_PACKAGE(OpenCV REQUIRED)
INCLUDE_DIRECTORIES(${OpenCV_INCLUDE_DIRS})

ADD_SUBDIRECTORY(deps/vecmath)

SET(FinalProject_SOURCES
        src/image.cpp
        src/main.cpp
        src/mesh.cpp
        src/scene_parser.cpp)

SET(FinalProject_INCLUDES
        include/camera.hpp
        include/group.hpp
        include/hit.hpp
        include/image.hpp
        include/light.hpp
        include/material.hpp
        include/mesh.hpp
        include/object3d.hpp
        include/plane.hpp
        include/ray.hpp
        include/scene_parser.hpp
        include/sphere.hpp
        include/transform.hpp
        include/triangle.hpp
        include/surfacelight.hpp
        include/meshtriangle.hpp
        include/texture.hpp
        )

SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

ADD_EXECUTABLE(${PROJECT_NAME} ${FinalProject_SOURCES} ${FinalProject_INCLUDES})

TARGET_LINK_LIBRARIES(${PROJECT_NAME} vecmath ${OpenCV_LIBS})

TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PRIVATE include)
