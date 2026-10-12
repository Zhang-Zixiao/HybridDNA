cmake_minimum_required(VERSION 3.14)
project(hybriddna_poa LANGUAGES CXX)
include(FetchContent)
set(spoa_build_exe OFF CACHE BOOL "" FORCE)
set(spoa_build_tests OFF CACHE BOOL "" FORCE)
set(spoa_install OFF CACHE BOOL "" FORCE)
set(spoa_optimize_for_native OFF CACHE BOOL "" FORCE)
set(spoa_optimize_for_portability OFF CACHE BOOL "" FORCE)
FetchContent_Declare(spoa
    GIT_REPOSITORY https://github.com/rvaser/spoa.git
    GIT_TAG d9d887b4cf0115a48ff1d12516db397df6d69064)
FetchContent_MakeAvailable(spoa)
add_executable(poa_consensus poa_consensus.cpp)
set_target_properties(poa_consensus PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)
target_link_libraries(poa_consensus PRIVATE spoa)
set_target_properties(poa_consensus PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_CURRENT_SOURCE_DIR}")
if(MINGW)
    target_link_options(poa_consensus PRIVATE -static)
endif()
