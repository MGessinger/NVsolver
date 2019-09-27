#ifndef REAL_H_
#define REAL_H_

#define REAL double

typedef struct lattice {
    int imax, jmax;
    int il, ir;
    int jb, jt;
    REAL delx;
    REAL dely;
} lattice;

typedef struct fluidSimulation {
    REAL Re;
    REAL tau;
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
