set(Tests
	sim_creation
	poisson
	derivative_du2dx
	derivative_dv2dy
	derivative_duvdx
	derivative_duvdy
	auxiliary
)

include_directories(${CMAKE_SOURCE_DIR}/src)

foreach (test_file IN LISTS Tests)
	set(this_test ${test_file}_target)

	add_executable(${this_test} ${test_file}.c)
	add_dependencies(${this_test} simlib)
	target_link_libraries(${this_test} simlib m)
	target_compile_options(${this_test} PRIVATE -Wall -Wextra -pedantic)

	add_test(NAME ${test_file}_test COMMAND ${test_file}_target)
endforeach()
