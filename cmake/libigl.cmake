if(TARGET igl::core)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG 10b3cd442e69e0a4e0287387453f7f0bf60ae4b2
)
FetchContent_MakeAvailable(libigl)
