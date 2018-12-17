cmake_minimum_required(VERSION 3.13)
project(torcs)
add_executable(torcs main.cpp carrier.cpp carrier.h)
include(FindArmadillo)
find_package(Armadillo REQUIRED)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /usr/bin/clang++ -std=gnu++11 -Wall")  # -c -o  -Xpreprocessor -fopenmp -lomp -o
message(${CMAKE_CXX_FLAGS})
link_directories(/usr/local/lib)
if (Armadillo_FOUND)
    include_directories(${ARMADILLO_INCLUDE_DIR})
    target_link_libraries(torcs ${ARMADILLO_LIBRARY})
endif ()








