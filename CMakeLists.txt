cmake_minimum_required(VERSION 3.0)
project(Simulator)

# just lists the sources that we need to build our lib
# note that "main" is exluded here
set(FIELDLIB_SOURCES
	src/fields.c
	src/fields.h
	src/poisson.c
	src/poisson.h
	src/boundary.c
	src/boundary.h
	src/IO.c
	src/IO.h
	src/particle.c
	src/particle.h
	src/types.h
        src/parallel.c
        src/parallel.h
)

set(OPTIONS
	-Wall
	-Wextra
	-pedantic
	-O3
)

# adds a library target called fieldlib to be built from the source files
add_library(fieldlib ${FIELDLIB_SOURCES})
target_include_directories(fieldlib PUBLIC src)
target_compile_options(fieldlib PRIVATE ${OPTIONS})

# Adds an executable target called simulator to be built from the source/main.c
add_executable(simulator src/main.c)
target_compile_options(simulator PRIVATE ${OPTIONS})
target_link_libraries(simulator PUBLIC fieldlib png)

find_package(MPI REQUIRED)
target_link_libraries(simulator PUBLIC MPI::MPI_C)


#add_subdirectory(tests)

