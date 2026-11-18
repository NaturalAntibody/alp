cmake_minimum_required(VERSION 3.15)
project(alp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# Link the MSVC C++ runtime statically: Python ships vcruntime140.dll but not
# msvcp140.dll, so a dynamically linked libalp.dll could fail to load.
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(ALP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib/alp)
set(ALP_SRCS
    ${ALP_DIR}/njn_dynprogprob.cpp
    ${ALP_DIR}/njn_dynprogproblim.cpp
    ${ALP_DIR}/njn_dynprogprobproto.cpp
    ${ALP_DIR}/njn_ioutil.cpp
    ${ALP_DIR}/njn_localmaxstat.cpp
    ${ALP_DIR}/njn_localmaxstatmatrix.cpp
    ${ALP_DIR}/njn_localmaxstatutil.cpp
    ${ALP_DIR}/njn_random.cpp
    ${ALP_DIR}/sls_alignment_evaluer.cpp
    ${ALP_DIR}/sls_alp.cpp
    ${ALP_DIR}/sls_alp_data.cpp
    ${ALP_DIR}/sls_alp_regression.cpp
    ${ALP_DIR}/sls_alp_sim.cpp
    ${ALP_DIR}/sls_basic.cpp
    ${ALP_DIR}/sls_pvalues.cpp
)
# Third-party public-domain NCBI code: silence its warnings.
if(MSVC)
    set_source_files_properties(${ALP_SRCS} PROPERTIES COMPILE_OPTIONS "/w")
else()
    set_source_files_properties(${ALP_SRCS} PROPERTIES COMPILE_OPTIONS "-w")
endif()

add_library(alp SHARED src/alp_api.cpp ${ALP_SRCS})
target_include_directories(alp PRIVATE ${ALP_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/src)
# File name is libalp.so / libalp.dylib / libalp.dll. The library is loaded with
# ctypes, so all symbols are exported on Windows, and the output goes straight
# into the build dir (the $<0:> genex suppresses per-config subdirectories).
set_target_properties(alp PROPERTIES
    PREFIX "lib"
    WINDOWS_EXPORT_ALL_SYMBOLS ON
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}$<0:>"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}$<0:>"
)
if(MSVC)
    target_compile_options(alp PRIVATE /O2)
else()
    target_compile_options(alp PRIVATE -O2)
endif()
if(NOT WIN32)
    target_link_libraries(alp PRIVATE m)
endif()

# Used by scikit-build-core to place the library inside the `alp` package.
install(TARGETS alp
    LIBRARY DESTINATION alp
    RUNTIME DESTINATION alp
)
