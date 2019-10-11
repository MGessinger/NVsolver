#include "particle.h"

particle* createParticleArray(int partcount)
{
    particle *parts = malloc(partcount*sizeof(particle));
    if (parts == NULL)
        return NULL;
    for (int i = 0; i < partcount; i++)
    {
        parts[i].onScreen = 0;
    }
    return parts;
}

void destroyParticleArray(particle *parts)
{
    free(parts);
    return;
}

void ParticleSeed (particle *parts, REAL posx1, REAL posx2, REAL posy1, REAL posy2, int partcount, int num)
{
    /* Place particles on the line going from (posx1,posy1) to (posx2,posy2) */
    if (num > partcount)
        num = partcount;
    REAL delx = (posx2 - posx1)/num;
    REAL dely = (posy2 - posy1)/num;
    for (int i = 0; i < partcount; i++)
    {
        if (i == num)
            break;
        if (parts[i].onScreen)
            continue;
        parts[i].onScreen = 1;
        parts[i].x = posx1;
        parts[i].y = posy1;
        posx1 += delx;
        posx2 += dely;
    }
    return;
}

void ParticleVelocity (REAL **U, REAL **V, particle *parts, lattice *grid, short **FLAG, int partcount)
{
    int i,j;
    REAL x2, y2;
    for (int p = 0; p < partcount; p++)
    {
        if (parts[p].onScreen == 0)
            continue;
        i = parts[p].x/grid->delx + 1;
        j = parts[p].y/grid->dely+1.5;
        if((i < 1 || i >= grid->imax) || (j < 1 || j >= grid->jmax))
        {
            parts[p].onScreen = 0;
            continue;
        }
        if (FLAG[i][j] != C_F)
        {
            parts[p].onScreen = 0;
            continue;
        }
        x2 = i*grid->delx;
        y2 = (j-0.5)*grid->dely;
        parts[p].u = (x2 - parts[p].x)*(y2 - parts[p].y)*U[i-1][j-1];
        parts[p].u += (parts[p].x - x2 - grid->delx)*(y2 - parts[p].y)*U[i][j-1];
        parts[p].u += (x2 - parts[p].x)*(parts[p].y - y2-grid->dely)*U[i-1][j];
        parts[p].u += (parts[p].x - x2- grid->delx)*(parts[p].y - y2 - grid->dely)*U[i][j];
        parts[p].u /= grid->delx*grid->dely;

        i = parts[p].x/grid->delx + 1.5;
        j = parts[p].y/grid->dely + 1;
        if((i < 1 || i >= grid->imax) || (j < 1 || j >= grid->jmax))
        {
            parts[p].onScreen = 0;
            continue;
        }
        if (FLAG[i][j] != C_F)
        {
            parts[p].onScreen = 0;
            continue;
        }
        x2 = (i - 0.5)*grid->delx;
        y2 = j*grid->dely;
        parts[p].v = (x2 - parts[p].x)*(y2 - parts[p].y)*V[i-1][j-1];
        parts[p].v += (parts[p].x - x2 - grid->delx)*(y2 - parts[p].y)*V[i][j-1];
        parts[p].v += (x2 - parts[p].x)*(parts[p].y - y2 - grid->dely)*V[i-1][j];
        parts[p].v += (parts[p].x - x2 - grid->delx)*(parts[p].y - y2 - grid->dely)*V[i][j];
        parts[p].v /= grid->delx*grid->dely;
    }
    return;
}

void ParticleTransport (particle *parts, int partcount, REAL delt)
{
    for (int p = 0; p < partcount; p++)
    {
        if (parts[p].onScreen == 0)
            continue;
        parts[p].x += parts[p].u*delt;
        parts[p].y += parts[p].v*delt;
    }
    return;
}
