/*
FindSepAttPoints.h
*/

#include "lib/icMatrix.h"

double CalDeterminant(icMatrix2x2 mat);
void GetEigenVectors(icMatrix2x2 mat, double evalues[2], icVector2 ev[2]);

void CalEigenVecsForSCC();
void CalAllSpecialPoints();

void CalEigenVecsForMesh();

void CalJacobianForWholeMesh();
void CalJacobianForVertex(int vertID);
void CalAllSpecialPointsForASCC(int scc_index);
void CalEigenForMatrix2x2(icMatrix2x2, icVector2 []);

void CalAllSpecialPoints_alltraingles();
