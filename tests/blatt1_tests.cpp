#include <iostream>

#include "gtest/gtest.h"

extern "C" {
// Fieldlib is a C library and requires C linkage
#include "blas.h"
#include "fields.h"
#include "real.h"
#include "IO.h"
}

TEST(Fieldlib,Blas)
{
    REAL **A = createMatrix(10,10);
    REAL *x = createVector(10);
    REAL *y = createVector(10);
    for (int i = 0; i < 10; i++)
    {
        for (int j = 0; j < 10; j++)
            A[i][j] = (i+1)*(j+1);
    }
    fill1Dfield(1,x,10);

    printf("A:\n");
    print2Dfield(A,10,10);
    gemv(1,A,x,0,y,10,10);
    printf("nrm2(A*x): %lg\n\n",nrm2(y,10));
    EXPECT_TRUE(isEqualScalar(nrm2(y,10),1079.177928,1e-06));

    copy(x,y,10);
    gemv(1,A,x,5,y,10,10);
    scal(1.5,x,10);
    printf("(5*x + A*x) * (1.5*x): %lg\n\n",dot(x,y,10));
    EXPECT_TRUE(isEqualScalar(dot(x,y,10),4612.500000,1e-06));
    destroyMatrix(A,10);
    destroyVector(x);
    destroyVector(y);
    return;
}
/*
TEST(Fieldlib, VTK){
    int sizeX = 60, sizeY = 40;
    REAL **A = create2Dfield(sizeX,sizeY);
    for (int i = 0; i < sizeX; i++)
    {
        for (int j = 0; j < sizeY; j++)
            A[i][j] = (REAL)i/60*M_PI;
    }
    writeVTKfileFor2DscalarField("gradientfield.vtk",
                                 "Gradientfield",
                                 A,sizeX,sizeY,
                                 0.1,0.1);
    applyFunctionTo2Dfield(sin,A,sizeX,sizeY);
    writeVTKfileFor2DscalarField("sinusfield.vtk",
                                 "Sinusfield",
                                 A,sizeX,sizeY,
                                 0.1,0.1);
    destroy2Dfield(A,sizeX);
    return;
}

TEST(Aufgabe_6,Teil_c){
    int sizeXU, sizeYU;
    int sizeXV, sizeYV;
    REAL **U = read2Dfield("fieldU.dat",&sizeXU,&sizeYU);
    REAL **V = read2Dfield("fieldV.dat",&sizeXV,&sizeYV);

    EXPECT_EQ(sizeXU,60);
    EXPECT_EQ(sizeXV,60);
    EXPECT_EQ(sizeYU,60);
    EXPECT_EQ(sizeYV,60);

    writeVTKfileFor2DvectorField("vectorfield.vtk",
                                 "Vectorfield",
                                 U, V, sizeXU, sizeYU,
                                 0.1, 0.1);
    destroy2Dfield(U,sizeXU);
    destroy2Dfield(V,sizeXV);
}*/

TEST(Fieldlib,IO)
{
    int sizeXA = 60, sizeYA = 60;
    int sizeXB = 0, sizeYB = 0;
    REAL **A = createMatrix(sizeXA,sizeYA);
    for (int i = 0; i < sizeXA; i++)
        for (int j = 0; j < sizeYA; j++)
            A[i][j] = M_2_PI*(i-j)/(i+j+1);
    write2Dfield("matrix.dat",A,sizeXA,sizeYA);
    REAL **B = read2Dfield("matrix.dat",&sizeXB,&sizeYB);
    EXPECT_EQ(sizeXA,sizeXB);
    EXPECT_EQ(sizeYA,sizeYB);
    EXPECT_TRUE(isEqual2Dfield(A,B,sizeXA,sizeYA,1e-20));

    destroyMatrix(A,sizeXA);
    destroyMatrix(B,sizeXB);
    return;
}
