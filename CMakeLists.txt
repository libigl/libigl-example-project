cmake_minimum_required(VERSION 2.8.6)

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# We need C++11. Put this directive after CGAL's include.
if(CMAKE_COMPILER_IS_GNUCXX)
    message(STATUS "GCC detected, adding compile flags")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -Wall -Wextra -Wno-reorder -Wno-unknown-pragmas -Wno-deprecated-declarations")
endif()
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif()

find_package(EIGEN REQUIRED)
include_directories( "${EIGEN_INCLUDE_DIR}" )
# Must be found before trying to find libraries in libigl/external
find_package(LIBIGL REQUIRED)
include_directories( "${LIBIGL_INCLUDE_DIR}" )
find_package(OpenGL REQUIRED)
find_package(GLFWH REQUIRED)
include_directories( ${GLFW_INCLUDE_DIR})
find_package(GLEW REQUIRED)
include_directories( ${GLEW_INCLUDE_DIR})

set(NANOGUI_BUILD_PYTHON OFF CACHE BOOL " " FORCE)
set(NANOGUI_BUILD_EXAMPLE OFF CACHE BOOL " " FORCE)
set(NANOGUI_BUILD_SHARED OFF CACHE BOOL " " FORCE)
find_package(NANOGUI QUIET)

add_subdirectory(   "${NANOGUI_INCLUDE_DIR}/../" "nanogui")
include_directories("${NANOGUI_INCLUDE_DIR}/../include")
include_directories("${NANOGUI_INCLUDE_DIR}/../ext/nanovg/src")
list(APPEND LIBS "nanogui")

if(APPLE)
	set(CMAKE_CXX_LINK_FLAGS "-framework OpenGL -framework Cocoa -framework IOKit -framework CoreVideo")
         set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.2")
endif (APPLE) #APPLE

link_directories(
  /usr/local/lib
  /opt/local/lib
  ${EIGEN_DIRS}
)

set(LIBS 
  ${LIBS}
  ${LIBIGL_LIBRARIES} 
  ${OPENGL_LIBRARIES} 
)

FILE(GLOB SRCFILES *.cpp)
add_executable(example_bin ${SRCFILES})
target_link_libraries(example_bin ${LIBS})
