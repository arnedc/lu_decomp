project(lu_decomp)
cmake_minimum_required(VERSION 2.8)
add_executable(lu_decomp solve_mat.c import_mat.c matrix_mult.c decomp_mat.c print_mat.c create_mat.c main.c)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

set(os Linux)
set(mach 64)
find_package(MKL REQUIRED)
include_directories(${MKL_INCLUDE_DIRS})
link_directories(${MKL_LIBRARIES})

set(LIBS ${LIBS} ${MKL_LIBRARIES})

target_link_libraries(lu_decomp ${LIBS})
include (${CMAKE_ROOT}/Modules/FindBLAS.cmake)

