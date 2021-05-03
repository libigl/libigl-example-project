if(TARGET igl::core)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG b0fd49d59840b5a48e79096462d3198a5c13eee9
)
FetchContent_MakeAvailable(libigl)
