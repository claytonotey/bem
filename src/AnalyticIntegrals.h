#ifndef ANALYTICINTEGRALS_H
#define ANALYTICINTEGRALS_H

#include "Vector3.h"
#include "ComplexVector3.h"
#include "MathUtils.h"

void z0_S0_S1_G0_G1_G3(Real *S0, Real *S1, Real *G0x, Real *G0y, Real *G1x, Real *G1y, Real *G3x, Real *G3y, Real v2x, Real v2y, Real x0, Real y0, Real v1x, Real v1y, Real yx0, Real yy0);

void z0_S0_S1_G0_G1(Real *S0, Real *S1, Real *G0x, Real *G0y, Real *G1x, Real *G1y, Real v2x, Real v2y, Real x0, Real y0, Real v1x, Real v1y, Real yx0, Real yy0);

Real z0_gdotG3(Real v2x, Real v2y);

void G3overlog_edgelimit(Real *G3x, Real *G3y, Real v2x, Real v2y, Real x0, Real y0, Real v1x, Real v1y);

void C3overlog_edgelimit(Real *C3z, Real v2x, Real v2y, Real x0, Real y0);

void G3(Real *G3x, Real *G3y, Real v2x, Real v2y, Real x0, Real y0, Real z0, Real v1x, Real v1y, Real yx0, Real yy0, Real yz0);

void C3(Real *C3z, Real v2x, Real v2y, Real x0, Real y0, Real z0);

void S3_G3(Real *S3, Real *G3x, Real *G3y, Real v2x, Real v2y, Real x0, Real y0, Real z0, Real v1x, Real v1y, Real yx0, Real yy0, Real yz0);

void S3_G3_C3(Real *S3, Real *G3x, Real *G3y, Real *C3z, Real v2x, Real v2y, Real x0, Real y0, Real z0, Real v1x, Real v1y, Real yx0, Real yy0, Real yz0);

void S0_S1_G0_G1_C0_C1(Real *S0, Real *S1, Real *G0x, Real *G0y, Real *G1x, Real *G1y, Real *C0z, Real *C1z, Real v2x, Real v2y, Real x0, Real y0, Real z0, Real v1x, Real v1y, Real yx0, Real yy0, Real yz0);

void S0_S1_G0_G1(Real *S0, Real *S1, Real *G0x, Real *G0y, Real *G1x, Real *G1y, Real v2x, Real v2y, Real x0, Real y0, Real z0, Real v1x, Real v1y, Real yx0, Real yy0, Real yz0);



#endif
