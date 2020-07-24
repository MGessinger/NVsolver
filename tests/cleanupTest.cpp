#include <iostream>

#include "gtest/gtest.h"
#include <mpi.h>

extern "C" {
	#include "types.h"
}

TEST(Cleanup, Finalize)
{
	EXPECT_EQ(MPI_Finalize(),MPI_SUCCESS);
}
