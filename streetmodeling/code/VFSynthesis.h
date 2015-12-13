////VFSynthesis.h

/*This file provide routines fulfilling the design of vector field through user interface
The functionalities it provides include add/edit elements (singular or regular), edit limit cycle
performing region smoothing

The inputs of this module are underneath mesh, user input (generally are some mouse positions)
The output is the vector field on each vertex of the underneath mesh (under both 2D local and global frame)
*/

#include "lib/icVector.h"
#include "lib/icMatrix.h"

////
//The followings are some possible useful routines from previous program
void AddElement(int type, double x, double y);
void InitEditBox(int index, double x, double y);

void InitialJacobian(icMatrix3x3 &Jacobian, int type);

void UpdateJacobian();
void UpdateTransform(int, int);
void UpdateRegElemTransform(int TransformType, int elem_id);
void UpdateSinElemTransform(int TransformType, int elem_id);
void UpdateEditBox(int);


void create_Jac_Field_Single(double a, double b, double c, 
							 double d, double x0, double y0); /*Jacobian based design */

/*------------------------------------------------------*/
////New routines for create the field through region smoothing
void SetVectorsAtTriangleforNewElem(int type, int triangle);
void SetVectorsForNewRegularElem(int type, int triangle, icVector2 global_vec);
void AddElement(int type, int triangle);
void BuildtheField();
void BuildRegionforAddElement();

/*------------------------------------------------------*/

//void AddEditBox(CElementList *celem, CPoint center);
void ScreenToWorld(double position[3], CPoint center);

/*------------------------------------------------------*/
//For regular element adding
void SetBasis(CPoint point, int type);  ////set the bases of the element
void SetBasis(int triangle, CPoint point, int type);
void SetDirect(CPoint point);  ////set the direction of the element
void SetTransform(
		double base[2],
		icVector2 Direct,
		icMatrix3x3 &transform_matrix, 
		icMatrix3x3 &transposeRot
	 );   ////set transformation matrix for convergent or divergent elements

/*------------------------------------------------------*/
//For special element adding
//void SetSpecBasis(CPoint point);
//void SetSpecDirect(CPoint point);

/*------------------------------------------------------*/
////Calculate the vector on the vertices of the triangular mesh
void getVector(double x, double y, icVector2 &vec, double &mag);
void CalVectors();
void GetLocalVector();
void GlobalToLocal(double v[2], icVector2 &iv, int face_id);

/*------------------------------------------------------*/
////Clear all the field with white noise
void ClearAllField();

////Rotate or reflect the whole field
void RotateTheWholeField(double);
void ReflectTheWholeField(double);

/*------------------------------------------------------*/
void SaveFieldPerVertex(const char *);
void ReloadFieldPerVertex(const char *);              ////here we suppose the file name is constant "savedfield.txt"

void SaveFieldElements(const char *filename);
void ReloadFieldElements(const char *);
void WriteSingularElem(FILE*, int);
void ReadSingularElem(FILE*, int);
void WriteRegularElem(FILE*, int);
void ReadRegularElem(FILE*, int);
void SaveBeforeCal();

void InitUnderneathMesh();
void NormalizeField();

void SaveFieldAfterSmoothing();

void CalVectors_formula(double t);

/////Testing routine
void TestRoutine();

