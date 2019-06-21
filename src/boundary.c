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

void setBCond(REAL **U, REAL **V, int imax, int jmax, boundaryCond *bCond)
{
    if (bCond == NULL) /* OUTFLOW on all four edges and trivial geometry is the default. */
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
    short flag;
    for (int i = 1; i <= imax; i++)
        for (int j = 1; j <= jmax; j++)
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

void setSpecBCond(REAL **U, REAL **V, int imax, int jmax, char *problem)
{
    if (problem == NULL)
        return;
    if (strcmp(problem,"Driven Cavity") == 0)
    {
        for (int i = 1; i <= imax; i++)
            U[i][jmax+1] = 2-U[i][jmax];
        return;
    }
    if (strcmp(problem,"Step") == 0)
    {
        for (int j = jmax/2; j <= jmax; j++)
        {
            U[0][j] = 1;
            V[0][j] = V[1][j];
        }
        return;
    }
    if (strcmp(problem,"Ramp") == 0)
    {
        for (int i = 2; i < imax/3; i++)
            V[i][0] = -1;
    }
    if (strcmp(problem,"Tunnel") == 0 || strcmp(problem,"Von Karman") == 0)
    {
        for (int j = 0; j <= jmax+1; j++)
        {
            U[0][j] = 1;
        }
        return;
    }
    return;
}

void initFlags(const char *problem, short **FLAG, int imax, int jmax)
{
    if (FLAG == NULL)
        return;
    /* Manually sets the flag field for arbitrary generalised geometries.
     * If the flags are read from a file, set problem to "Image" */
    if (strcmp(problem,"Step") == 0)
    {
        for (int i = 0; i < jmax/2; i++)
            for (int j = 0; j < jmax/2; j++)
                FLAG[i][j] = C_B;
    }
    else if (strcmp(problem,"Von Karman") == 0){
        for (int i = jmax/3; i < 2*jmax/3; i++)
            for (int j = -(jmax/4)/2; j <= (jmax/4)/2; j++)
            {
                if (i+j < jmax/3 || i+j >= 2*jmax/3)
                    continue;
                FLAG[i+j][i] = C_B;
            }
    }
    else if (strcmp(problem,"Ramp") == 0)
    {
        for (int j = jmax/2+5; j < jmax; j++)
        {
            FLAG[0][j] = FLAG[1][j] = C_B;
            FLAG[imax/3][j] = FLAG[imax/3+1][j] = C_B;
        }
        FLAG[0][jmax/2+1] = FLAG[0][jmax/2+2] = C_B;
        FLAG[0][jmax/2+3] = FLAG[0][jmax/2+4] = C_B;
        FLAG[1][jmax/2+1] = FLAG[1][jmax/2+2] = C_B;
        FLAG[1][jmax/2+3] = FLAG[1][jmax/2+4] = C_B;
        for (int i = 0; i < imax/3+5; i++)
            FLAG[i][jmax/2-1] = FLAG[i][jmax/2] = C_B;

    }
    for (int i = 0; i < imax; i++)
        for (int j = 0; j < jmax; j++)
        {
            if (FLAG[i][j] == C_F)
                continue;
            if (FLAG[i][j+1] == C_F)
                FLAG[i][j] |= B_N;
            else if (j != 0 && FLAG[i][j-1] == C_F)
                FLAG[i][j] |= B_S;
            if (FLAG[i+1][j] == C_F)
                FLAG[i][j] |= B_O;
            else if (i != 0 && FLAG[i-1][j] == C_F)
                FLAG[i][j] |= B_W;
        }
    return;
}
