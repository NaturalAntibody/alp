cmake_minimum_required(VERSION 3.14)
project(alp CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE)
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
set_source_files_properties(${ALP_SRCS} PROPERTIES COMPILE_OPTIONS "-w")

add_library(alp SHARED src/alp_api.cpp ${ALP_SRCS})
target_include_directories(alp PRIVATE ${ALP_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_compile_options(alp PRIVATE -O2)
target_link_libraries(alp PRIVATE m)
