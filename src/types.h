#ifndef TYPES_H_
#define TYPES_H_

typedef enum {
	STICK,
	SLIP,
	OUTFLOW,
	SPECIAL
} boundary;

#define B_NORTH	(0x10)
#define B_EAST	(0x20)
#define B_SOUTH	(0x40)
#define B_WEST	(0x80)

typedef struct lattice {
	double dx, dy;		/* Step Width */
	int imax, jmax;		/* Grid size */
	boundary bc_top, bc_bottom;
	boundary bc_right, bc_left;
} lattice;

typedef struct fluid {
	double **U, **V, **P;	/* Velocity and pressure */
	double Reynolds;
} fluid;

typedef struct simulation {
	lattice * G;
	fluid * F;
	double **auxF, **auxG;	/* Auxiliary fields */
	double **RHS;		/* Poisson Eq. RHS */
	double GX, GY;		/* External forces */
	double tau, dt;		/* Time parameters */
	double alpha, omega;	/* Solver Parameters */
	double solver_tol;	/* Poisson tolerance */
	int max_iterations;
} simulation;

static inline double sqr (double x) {
	return x * x;
}

static inline double maxd (double a, double b) {
	if (a < b)
		return b;
	else
		return a;
}

static inline double mind (double a, double b) {
	if (a < b)
		return a;
	else
		return b;
}

static inline double absd (double a) {
	return maxd(a, -a);
}

#endif /* TYPES_H_ */
