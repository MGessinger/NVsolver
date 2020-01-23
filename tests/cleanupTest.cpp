#include <iostream>

#include "gtest/gtest.h"
#include <mpi/mpi.h>

extern "C" {
	#include "types.h"
}

TEST(Cleanup, FInalize)
{
	EXPECT_EQ(MPI_Finalize(),MPI_SUCCESS);
}
