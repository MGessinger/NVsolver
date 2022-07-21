#ifndef TYPES_H_
#define TYPES_H_

typedef struct lattice {
	double dx, dy;	/* Step Width */
	int imax, jmax;
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
	if (a >= b)
		return a;
	else
		return b;
}

static inline double mind (double a, double b) {
	if (a >= b)
		return b;
	else
		return a;
}

static inline double absd (double a) {
	if (a >= 0)
		return a;
	else
		return -a;
}

#endif /* TYPES_H_ */
