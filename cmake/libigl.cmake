if(TARGET igl::core)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG v2.3.0
)

# Note: In libigl v3.0.0, the following will become a one-liner:
# FetchContent_MakeAvailable(libigl)

FetchContent_GetProperties(libigl)
if(NOT libigl_POPULATED)
    FetchContent_Populate(libigl)
endif()
list(PREPEND CMAKE_MODULE_PATH "${libigl_SOURCE_DIR}/cmake")
include(${libigl_SOURCE_DIR}/cmake/libigl.cmake)
