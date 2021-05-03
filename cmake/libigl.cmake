if(TARGET igl::core)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG ac42b6de0d13aabc61a7ce3a2a95b914aaca383d
)
FetchContent_MakeAvailable(libigl)
