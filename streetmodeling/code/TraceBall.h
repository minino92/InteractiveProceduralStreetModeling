#pragma once
#include <math.h>
#include <malloc.h>

#define RENORMCOUNT 97
#define TRACKBALLSIZE  (0.8)
#define X 0
#define Y 1
#define Z 2
#define W 3

typedef float Quaternion[4];
typedef float Matrix[4][4];


class CTraceBall
{
public:
	CTraceBall(void);
	~CTraceBall(void);

	void quat_to_mat(Quaternion q, Matrix mat);
    void mat_to_quat(Matrix mat, Quaternion q);
    float *this_vnew();
    void this_vcopy(const float *v1, float *v2);
    void this_vset(float *v, float x, float y, float z);
    void this_vzero(float *v);
    float this_vlength(const float *v);
    void this_vscale(float *v, float div);
    void this_vmult(const float *src1, const float *src2, float *dst);
    void this_vnormal(float *v);
   void this_vadd(const float *src1, const float *src2, float *dst);
    void this_vsub(const float *src1, const float *src2, float *dst);
    float this_vdot(const float *v1, const float *v2);
    void this_vcross(const float *v1, const float *v2, float *cross);
    void vhalf(const float *v1, const float *v2, float *half);
    void vdirection(const float *v1, float *dir);
    void vreflect(const float *in, const float *mirror, float *out);
    void vmultmatrix(const Matrix m1, const Matrix m2, Matrix prod);
    void vtransform(const float *v, const Matrix mat, float *vt);
    void vtransform4(const float *v, const Matrix mat, float *vt);
    void mcopy(const Matrix m1, Matrix m2);
    void minvert(const Matrix mat, Matrix result);

    float tb_project_to_sphere(float, float, float);
    void normalize_quat(float [4]);

	void trackball(float q[4], float p1x, float p1y, float p2x, float p2y);
    void axis_to_quat(float a[3], float phi, float q[4]);

	void add_quats(float q1[4], float q2[4], float dest[4]);

	void build_rotmatrix(float m[4][4], float q[4]);


};
