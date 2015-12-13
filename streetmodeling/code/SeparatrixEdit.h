
////SeparatrixEdit.h

#include "lib/icvector.h"


void SeparatrixModify(int saddleID, int trajID);

void GetInitSeparatrixStrip(int saddleID, int sepindex);

void GetTheRegionForSeparatrixModification(int inorout);

void GetInnerVertices();

void SetVectorsOnNewBoundary(int inorout);

void AllocSeparatrixEdit();

void InitSeparatrixEdit();

void Separatrix_Growing(int type);

void DiscretizeSeparatrix(int trajID);  ////Discretize the separatrix to get bunch of controlling points

void GetBoundaryofTriangleStrip();

void GetTheNewTriangleStrip();


void SetVectorOnVertices(int inorout);

icVector2 GetCurveDirection(int triangleID);

icVector2 GetAVector(icVector2 curveorient, icVector2 normal, double ang, int inorout);
