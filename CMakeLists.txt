cmake_minimum_required(VERSION 3.10)
project(libntl_julia CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(JlCxx REQUIRED)

# Find NTL
find_path(NTL_INCLUDE_DIR NAMES NTL/ZZ.h)
find_library(NTL_LIBRARY NAMES ntl)

if(NOT NTL_INCLUDE_DIR)
    message(FATAL_ERROR "NTL include directory not found")
endif()

if(NOT NTL_LIBRARY)
    message(FATAL_ERROR "NTL library not found")
endif()

add_library(libntl_julia SHARED libntl_julia.cpp)

target_link_libraries(libntl_julia
    PRIVATE JlCxx::cxxwrap_julia
    PRIVATE ${NTL_LIBRARY}
)

target_include_directories(libntl_julia
    PRIVATE ${NTL_INCLUDE_DIR}
)

target_compile_definitions(libntl_julia PRIVATE NTL_EXCEPTIONS=on)

set_target_properties(libntl_julia PROPERTIES
    PREFIX ""
    OUTPUT_NAME "libntl_julia"
)

install(TARGETS libntl_julia
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)
