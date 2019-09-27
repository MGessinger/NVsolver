cmake_minimum_required(VERSION 3.0)
project(Simulation)

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
)

# adds a library target called fieldlib to be built from the source files
add_library(fieldlib ${FIELDLIB_SOURCES})
# Specifies the include directory to use when compiling our target fieldlib.
target_include_directories(fieldlib PUBLIC src)

# Catch Bugs Before Testing
# Compiler warnings are often indicators of future bugs
# that you would otherwise see only at runtime.
if(MSVC)
  target_compile_options(fieldlib PRIVATE /W4 /WX)
else(MSVC)
  target_compile_options(fieldlib PRIVATE -Wall -Wextra -pedantic -g2)
endif(MSVC)

# Adds an executable target called blatt1 to be built from the source/main.c
add_executable(simulator src/main.c)
# links our executable with the source library fieldlib
target_link_libraries(simulator fieldlib)
# Add a subdirectory to the build (where another CMakeLists.txt is located)
add_subdirectory(tests)

