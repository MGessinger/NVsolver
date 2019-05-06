#ifndef REAL_H_
#define REAL_H_

#define REAL double

typedef struct lattice {
    REAL xlength;
    REAL ylength;
    int imax;
    int jmax;
} lattice;

typedef struct fluidSimulation {
    REAL Re;
    REAL t_end;
    REAL GX;
    REAL GY;
    REAL eps;
    REAL omega;
    REAL alpha;
    int itmax;
} fluidSim;

#define OUTPUT (0x10)
#define PRINT  (0x01)
#define SILENT (0x0)

#endif /* REAL_H_ */
