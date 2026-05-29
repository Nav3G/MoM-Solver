cmake_minimum_required(VERSION 3.14)
project(mom_solver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)
FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG        3.4.0
)
FetchContent_MakeAvailable(eigen)

add_executable(mom_solver src/main.cpp src/io.cpp src/geometry.cpp src/quadrature.cpp src/matrix_fill_v2.cpp src/excitation_v2.cpp)
target_include_directories(mom_solver PRIVATE include)
target_link_libraries(mom_solver PRIVATE Eigen3::Eigen)