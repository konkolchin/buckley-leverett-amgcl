cmake_minimum_required(VERSION 3.16)
project(buckley_leverett LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# AMGCL is header-only. Its own CMakeLists always searches for Boost, which this
# program does not need (AMGCL_NO_BOOST). Point the include path at a local checkout.
set(AMGCL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/third_party/amgcl" CACHE PATH "AMGCL source root")
if(NOT EXISTS "${AMGCL_DIR}/amgcl/make_solver.hpp")
    message(FATAL_ERROR
        "AMGCL headers not found in ${AMGCL_DIR}. Clone them with:\n"
        "  git clone --depth 1 --branch 1.4.3 https://github.com/ddemidov/amgcl.git third_party/amgcl")
endif()

add_executable(buckley src/buckley_leverett.cpp)
target_include_directories(buckley PRIVATE ${AMGCL_DIR})
target_compile_definitions(buckley PRIVATE AMGCL_NO_BOOST NOMINMAX)

if(MSVC)
    target_compile_options(buckley PRIVATE /W3 /bigobj)
else()
    target_compile_options(buckley PRIVATE -Wall -Wextra -Wpedantic)
endif()
