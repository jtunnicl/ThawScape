#
# CMake script to run a single test case
#
message(STATUS "Running GrateRipCLI test")
message(STATUS "  Test run directory: ${TEST_RUN_DIR}")
message(STATUS "  Test src directory: ${TEST_SRC_DIR}")
message(STATUS "  Test ref directory: ${TEST_REF_DIR}")
message(STATUS "  Test binary: ${TEST_BINARY}")

#
# make the test directory
#
execute_process(COMMAND ${CMAKE_COMMAND} -E remove_directory ${TEST_RUN_DIR})
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${TEST_RUN_DIR})

#
# copy input files
#
file(COPY ${TEST_SRC_DIR}/ThawScape.ini DESTINATION ${TEST_RUN_DIR})
file(COPY ${TEST_SRC_DIR}/topo.asc DESTINATION ${TEST_RUN_DIR})
file(COPY ${TEST_SRC_DIR}/FA.asc DESTINATION ${TEST_RUN_DIR})

#
# run the code
#
execute_process(
    COMMAND ${CMAKE_COMMAND} -E chdir ${TEST_RUN_DIR} ${TEST_BINARY}
    RESULT_VARIABLE status
)
if (status)
    message(FATAL_ERROR "Error running ThawScape: '${status}'")
endif (status)

#
# check the results
#
# first get the list of reference files and then loop over them
message(STATUS "Comparing result with reference...")
file(GLOB REF_FILES RELATIVE ${TEST_REF_DIR} ${TEST_REF_DIR}/erosion_*.asc)
foreach (FILENAME ${REF_FILES})
    message(STATUS "Comparing file: ${FILENAME}")

    # first check the corresponding file exists in the run directory
    if (EXISTS "${TEST_RUN_DIR}/${FILENAME}")
        # then compare the two file
        execute_process(
            COMMAND ${CMAKE_COMMAND} -E compare_files --ignore-eol
                "${TEST_REF_DIR}/${FILENAME}" "${TEST_RUN_DIR}/${FILENAME}"
            RESULT_VARIABLE compare_result
        )
        if (compare_result)
            message(FATAL_ERROR "Output file does not match reference: ${FILENAME}")
        endif()
    else()
        message(FATAL_ERROR "Output file does not exist in run directory: ${FILENAME}")
    endif()
endforeach()
