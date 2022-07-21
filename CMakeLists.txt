cmake_minimum_required(VERSION 3.0)
project(NVS2D
	LANGUAGES C
	DESCRIPTION "A two-dimensional Navier-Stokes solver on a staggered grid"
)

# just list the sources that we need to build our lib
set(simlib_source
	src/boundary.c
	src/poisson.c
	src/setup.c
	src/velocities.c
	src/IO.c
)

# adds a library target called fieldlib to be built from the source files
add_library(simlib ${simlib_source})
target_include_directories(simlib PUBLIC src)
target_compile_options(simlib PRIVATE -Wall -Wextra -pedantic)

# Adds an executable target called simulator to be built from src/main.c
add_executable(simulator src/main.c)
target_link_libraries(simulator PUBLIC simlib)

enable_testing()
add_subdirectory(tests)
