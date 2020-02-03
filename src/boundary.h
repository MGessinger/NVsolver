#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include "types.h"

#define C_F (0x0F)
#define C_B (0x00)
#define B_N (0x01)
#define B_O (0x02)
#define B_S (0x04)
#define B_W (0x08)
#define B_NS (0x05)
#define B_OW (0x0a)

#define LO_NIBBLE(x) ((x) & 0x0F)
#define HI_NIBBLE(x) ((x) & 0xF0)

#define NOSLIP   (0u)
#define FREESLIP (1u)
#define OUTFLOW  (2u)

bndCond createBoundCond (int wl, int wr, int wt, int wb);
/* Create and destroy a bndCond structure */

void    applyPbndCond (REAL **P, lattice *grid, char **FLAG);
void    setBCond (REAL **U, REAL **V, lattice *grid, bndCond *bCond);
void    setSpecBCond (REAL **U, REAL **V, lattice *grid, const char *problem);
/* Set the boundary conditions for the given problem */

void    initFlags (const char *problem, char **FLAG, lattice *grid, MPI_Comm Region);
/* Set flags for fluid and obstacle cells */

#endif /* BOUNDARY_H_ */
