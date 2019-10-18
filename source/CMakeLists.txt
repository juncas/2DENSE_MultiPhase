# CMake 最低版本号要求
cmake_minimum_required (VERSION 2.8)

# 项目信息
project (2DENSE)

# 启动FORTRAN语言
enable_language(Fortran)
set(CMAKE_Fortran_COMPILER "ifort")

find_package(OpenMP REQUIRED)

set(HDF5_ROOT "/home/pengjun/hdf5")
find_package(HDF5 COMPONENTS Fortran REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})

set(CMAKE_Fortran_FLAGS "-O2 ${CMAKE_Fortran_FLAGS}")

# 将所有FORTRAN文件插入
file(GLOB_RECURSE src_files ./*.f90)

# 生成目标文件
add_executable("2DENSE.out" ${src_files})
target_link_libraries("2DENSE.out" ${HDF5_Fortran_LIBRARIES} ${OpenMP_Fortran_LIBRARIES})
