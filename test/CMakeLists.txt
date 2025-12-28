file(GLOB_RECURSE SRCS *.cpp)

add_executable(dicer_test ${SRCS})
target_link_libraries(dicer_test dicer)
