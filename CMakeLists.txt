cmake_minimum_required(VERSION 3.27)

project(free_fermions LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)

set(CMAKE_CXX_FLAGS "-pthread")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g") 
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# -------------------------------- #
# ------ Fetching externals ------ #
# -------------------------------- #

include(FetchContent)

set(EXTERNALS "")

if (DEFINED DATAFRAME_SOURCE)
    set(dataframe_SOURCE_DIR ${DATAFRAME_SOURCE})
    add_subdirectory(${DATAFRAME_SOURCE} "dataframe_src")
else()
    FetchContent_Declare(
        dataframe
        GIT_REPOSITORY https://github.com/eliotheinrich/dataframe
        GIT_TAG master
        GIT_SHALLOW TRUE
    )
    list(APPEND EXTERNALS dataframe)
endif()

if (DEFINED UTILS_SOURCE)
    set(utils_SOURCE_DIR ${UTILS_SOURCE})
    add_subdirectory(${UTILS_SOURCE} "utils_src")
else()
    FetchContent_Declare(
        utils 
        GIT_REPOSITORY https://github.com/eliotheinrich/utils
        GIT_TAG main 
        GIT_SHALLOW TRUE
    )
    list(APPEND EXTERNALS utils)
endif()

FetchContent_Declare(
  glaze
  GIT_REPOSITORY https://github.com/stephenberry/glaze.git
  GIT_TAG main
)
list(APPEND EXTERNALS glaze)

FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen
    GIT_TAG master
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE
)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
set(EIGEN_BUILD_DOC OFF)
set(BUILD_TESTING OFF)
set(EIGEN_BUILD_PKGCONFIG OFF)
list(APPEND EXTERNALS eigen)

FetchContent_MakeAvailable(${EXTERNALS})

include_directories(
    ${dataframe_SOURCE_DIR}/src
    ${eigen_SOURCE_DIR}
    ${utils_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/src
    ${glaze_SOURCE_DIR}/include
)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(fmt)
if (BUILDING_PYFERMION)
    find_package(Python 3.10 
        REQUIRED COMPONENTS Interpreter Development.Module
        OPTIONAL_COMPONENTS Development.SABIModule
    )

    # Finding nanobind
    execute_process(
        COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
        OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE NB_DIR
    )
    list(APPEND CMAKE_PREFIX_PATH "${NB_DIR}")
    find_package(nanobind CONFIG REQUIRED)

    nanobind_add_module(
        pyfermion_bindings 
        STABLE_ABI
        NB_STATIC
        src/PyFermion.cpp
    )

    target_link_libraries(pyfermion_bindings PUBLIC dataframe)

    install(
        TARGETS pyfermion_bindings
        LIBRARY DESTINATION pyfermion 
    )
else()
    add_executable(test src/test.cpp)
    target_link_libraries(test PUBLIC fmt::fmt simulator_display)
endif()
