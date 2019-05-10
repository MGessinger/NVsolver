#include "boundary.h"

void applyHomogeneousNeumannBC(REAL **p, int imax, int jmax)
{
    if (p == NULL)
        return;
    /* Set the corners to zero because they are irrelevant */
    p[0][0] = p[0][jmax+1] = p[imax+1][0] = p[imax+1][jmax+1] = 0;
    /* Ghost cells just copy the value next to them */
    for (int i = 1; i <= imax; i++)
    {
        p[i][0] = p[i][1];
        p[i][jmax+1] = p[i][jmax];
    }
    for (int j = 1; j <= jmax; j++)
    {
        p[0][j] = p[1][j];
        p[imax+1][j] = p[imax][j];
    }
    return;
}

void applyHomogeneousDirichletBC(REAL **p, int imax, int jmax)
{
    if (p == NULL)
        return;
    /* Set the corners to zero because they are irrelevant */
    p[0][0] = p[0][jmax+1] = p[imax+1][0] = p[imax+1][jmax+1] = 0;
    /* Ghost cells copy the negative value of the cell next to them */
    for (int i = 1; i <= imax; i++)
    {
        p[i][0] = -p[i][1];
        p[i][jmax+1] = -p[i][jmax];
    }
    for (int j = 1; j <= jmax; j++)
    {
        p[0][j] = -p[1][j];
        p[imax+1][j] = -p[imax][j];
    }
    return;
}

void setBCond(double **U, double **V, int imax, int jmax, boundaryCond *bCond)
{
    if (bCond == NULL) /* OUTFLOW on all four edges is the default. */
    {
        applyHomogeneousNeumannBC(U,imax-1,jmax);
        applyHomogeneousNeumannBC(V,imax,jmax-1);
        return;
    }
    for (int i = 1; i <= imax; i++)
    {
        U[i][0] = (bCond->wb == NOSLIP) ? -U[i][1] : U[i][1];
        U[i][jmax+1] = (bCond->wt == NOSLIP) ? -U[i][jmax] : U[i][jmax];
        V[i][0] = (bCond->wb == OUTFLOW) ? V[i][1] : 0;
        V[i][jmax] = (bCond->wt == OUTFLOW) ? V[i][jmax-1] : 0;
    }
    for (int j = 1; j <= jmax; j++)
    {
        U[0][j] = (bCond->wl == OUTFLOW) ? U[1][j] : 0;
        U[imax][j] = (bCond->wr == OUTFLOW) ? U[imax-1][j] : 0;
        V[0][j] = (bCond->wl == NOSLIP) ? -V[1][j] : V[1][j];
        V[imax+1][j] = (bCond->wr == NOSLIP) ? -V[imax][j] : V[imax][j];
    }
    return;
}

void setSpecBCond(REAL **U, REAL **V, int imax, int jmax, char *problem)
{
    if (problem == NULL)
        return;
    if (strcmp(problem,"Driven Cavity") == 0)
    {
        for (int i = 1; i <= imax; i++)
            U[i][jmax+1] = 2-U[i][jmax];
    }
    return;
    V = V;
}
