cmake_minimum_required(VERSION 3.10)

project(myproject)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_STANDARD 17)

include(tools)

add_external_package(Eigen3 VERSION 3.4.0 CONFIG)

set(MYPROJECT_HDRS
        # Add header files here
        # ${CMAKE_CURRENT_SOURCE_DIR}/my_file.h
        )

set(MYPROJECT_SRCS
        # Add source files here
        # ${CMAKE_CURRENT_SOURCE_DIR}/my_file.cpp
        )

set(MYPROJECT_executable_SRCS
        ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
        )

add_subdirectory(tests)

# Simulation code
add_executable(myproject ${MYPROJECT_SRCS} ${MYPROJECT_HDRS} ${MYPROJECT_executable_SRCS})
target_link_libraries(myproject PUBLIC Eigen3::Eigen)

# Documentation
option(BUILD_DOC, "Build documentation" OFF)
find_package(Doxygen)
if (DOXYGEN_FOUND)
    # Find include directories required to compile against the headers of the target
    #get_target_property(MYPROJECT_PUBLIC_HEADER_DIR myproject INTERFACE_INCLUDE_DIRECTORIES)
    # Find header files in include directories
    #file(GLOB_RECURSE MYPROJECT_PUBLIC_HEADERS ${MYPROJECT_PUBLIC_HEADER_DIR}/*.h)
    set(DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR})
    # Doxygen documentation will be saved in a subdirectory of the build directory
    set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen)
    file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR})
    set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/html/index.html)
    # Main Doxygen controls are in Doxyfile.in
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    add_custom_target(doc_doxygen ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM)
else (DOXYGEN_FOUND)
    message("Need Doxygen to build the documentation")
endif (DOXYGEN_FOUND)
