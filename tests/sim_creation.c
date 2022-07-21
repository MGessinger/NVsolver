#include <stdlib.h>
#include <time.h>
#include "setup.h"

int main () {
	srand(time(NULL));

	for (int it = 0; it < 20; it++) {
		int imax = rand() % 256;
		int jmax = rand() % 256;
		simulation * S = newSimulation(imax, jmax, 0.01, 0.01, 1000);

		clearSimulation(S);
	}

	return EXIT_SUCCESS;
}
