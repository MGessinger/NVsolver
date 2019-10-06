#include "boundary.h"

boundaryCond* createBoundCond(int wl, int wr, int wt, int wb)
{
    boundaryCond *bCond = malloc(sizeof(boundaryCond));
    if (bCond == NULL)
        return NULL;
    bCond->wl = wl;
    bCond->wr = wr;
    bCond->wt = wt;
    bCond->wb = wb;
    bCond->FLAG = NULL;
    return bCond;
}

void destroyBoundCond(boundaryCond *bCond, int imax)
{
    if (bCond == NULL)
        return;
    destroy2DIntegerField(bCond->FLAG,imax);
    free(bCond);
    return;
}

void applyHomogeneousNeumannBC(REAL **p, int imax, int jmax)
{
    if (p == NULL)
        return;
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

void setBCond(REAL **U, REAL **V, lattice *grid, boundaryCond *bCond)
{
    if (bCond == NULL) /* OUTFLOW on all four edges and trivial geometry is the default. */
    {
        applyHomogeneousNeumannBC(U,grid->imax-1,grid->jmax);
        applyHomogeneousNeumannBC(V,grid->imax,grid->jmax-1);
        return;
    }
    /* Boundary conditions on the actual boundary of the region */
    if (grid->edges & BOTTOM)
        for (int i = 1; i <= grid->deli; i++)
        {
            U[i][0] = (bCond->wb == NOSLIP) ? -U[i][1] : U[i][1];
            V[i][0] = (bCond->wb == OUTFLOW) ? V[i][1] : 0;
        }
    if (grid->edges & TOP)
        for (int i = 1; i <= grid->deli; i++)
        {
            U[i][grid->delj+1] = (bCond->wt == NOSLIP) ? -U[i][grid->delj] : U[i][grid->delj];
            V[i][grid->delj] = (bCond->wt == OUTFLOW) ? V[i][grid->delj-1] : 0;
        }
    if (grid->edges & LEFT)
        for (int j = 1; j <= grid->delj; j++)
        {
            U[0][j] = (bCond->wl == OUTFLOW) ? U[1][j] : 0;
            V[0][j] = (bCond->wl == NOSLIP) ? -V[1][j] : V[1][j];
        }
    if (grid->edges & RIGHT)
        for (int j = 1; j <= grid->delj; j++)
        {
            U[grid->deli][j] = (bCond->wr == OUTFLOW) ? U[grid->deli-1][j] : 0;
            V[grid->deli+1][j] = (bCond->wr == NOSLIP) ? -V[grid->deli][j] : V[grid->deli][j];
        }
    short flag;
    /* Boundary condtions on obstacles */
    for (int i = 1; i <= grid->deli; i++)
        for (int j = 1; j <= grid->delj; j++)
        {
            flag = bCond->FLAG[i-1][j-1];
            if (flag == C_F)
                continue;
            /* Deep inside an obstacle */
            if (flag == C_B)
            {
                U[i][j] = V[i][j] = 0;
                continue;
            }
            /* Set boundary conditions */
            switch (flag^C_B)
            {
            case B_N:
                V[i][j] = 0;
                U[i][j] = -U[i][j+1];
                U[i-1][j] = -U[i-1][j+1];
                break;
            case B_O:
                U[i][j] = 0;
                V[i][j] = V[i+1][j];
                V[i][j-1] = V[i+1][j-1];
                break;
            case B_S:
                V[i][j-1] = 0;
                U[i][j] = -U[i][j-1];
                U[i-1][j] = -U[i-1][j-1];
                break;
            case B_W:
                U[i-1][j] = 0;
                V[i][j] = V[i-1][j];
                V[i][j-1] = V[i-1][j-1];
                break;
            case B_N | B_O:
                U[i][j] = V[i][j] = 0;
                U[i-1][j] = -U[i-1][j+1];
                V[i][j-1] = -V[i+1][j-1];
                break;
            case B_N | B_W:
                U[i-1][j] = V[i][j] = 0;
                U[i-1][j] = -U[i-1][j+1];
                V[i][j-1] = V[i-1][j-1];
                break;
            case B_S | B_O:
                U[i][j] = V[i][j-1] = 0;
                U[i-1][j] = -U[i-1][j-1];
                V[i][j-1] = V[i+1][j-1];
                break;
            case B_S | B_W:
                U[i-1][j] = V[i][j-1] = 0;
                U[i-1][j] = -U[i-1][j-1];
                V[i][j-1] = V[i-1][j-1];
                break;
            default:
                U[i][j] = V[i][j] = 0;
                break;
            }
        }
    return;
}

void setSpecBCond(REAL **U, REAL **V, lattice *grid, const char *problem)
{
    /* Set special (e.g. inflow) conditions */
    if (problem == NULL)
        return;
    if (strcmp(problem,"Driven Cavity") == 0)
    {
        if (!(grid->edges & TOP))
            return;
        for (int i = 1; i <= grid->deli; i++)
            U[i][grid->delj+1] = 2-U[i][grid->delj];
        return;
    }
    if (strcmp(problem,"Step") == 0)
    {
        if (!(grid->edges & LEFT))
            return;
        for (int j = 0; j <= grid->delj; j++)
        {
            if (j < grid->jmax/2)
                continue;
            U[0][j] = 1;
            V[0][j] = V[1][j];
        }
        return;
    }
    if (strcmp(problem,"Tunnel") == 0 || strcmp(problem,"Von Karman") == 0)
    {
        if (!(grid->edges & LEFT))
            return;
        for (int j = 0; j <= grid->delj; j++)
        {
            U[0][j] = 1;
        }
        return;
    }
    return;
}

void initFlags(const char *problem, short **FLAG, lattice *grid)
{
    /* Manually sets the flag field for arbitrary generalised geometries.
     * If the flags are read from a file, set problem to "Image"! */
    if (!FLAG || !grid)
        return;
    if (strcmp(problem,"Step") == 0)
    {
        for (int i = 0; i < grid->deli; i++)
        {
            if (i + grid->il >= grid->jmax/2)
                break;
            for (int j = 0; j < grid->delj; j++)
            {
                if (j + grid->jb >= grid->jmax/2)
                    break;
                FLAG[i][j] = C_B;
            }
        }
    }
    else if (strcmp(problem,"Von Karman") == 0)
    {
        for (int i = 0; i < grid->delj; i++)
        {
            if (i + grid->il < grid->jmax/3)
                continue;
            if (i +grid->il >= 2*grid->jmax/3)
                break;
            for (int j = -(grid->jmax/4)/2; j <= (grid->jmax/4)/2; j++)
            {
                int x = i+j+grid->il+grid->jb;
                if (x < grid->jmax/3)
                    continue;
                if (x >= 2*grid->jmax/3 || x >= grid->deli)
                    break;
                FLAG[i+j][i] = C_B;
            }
        }
    }
    for (int i = 0; i < grid->deli; i++)
    {
        for (int j = 0; j < grid->delj; j++)
        {
            if (FLAG[i][j] == C_F)
                continue;
            if ((j+1) != grid->delj && FLAG[i][j+1] == C_F)
                FLAG[i][j] |= B_N;
            else if (j != 0 && FLAG[i][j-1] == C_F)
                FLAG[i][j] |= B_S;
            if ((i+1) != grid->delj && FLAG[i+1][j] == C_F)
                FLAG[i][j] |= B_O;
            else if (i != 0 && FLAG[i-1][j] == C_F)
                FLAG[i][j] |= B_W;
        }
    }
    return;
}
