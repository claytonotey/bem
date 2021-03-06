#ifndef TETRAHEDRONQUADRATURE_H
#define TETRAHEDRONQUADRATURE_H

#include "Vector3.h"
#include "FunctionArgs.h"

template<class T, class F>
T tetquad1(F f, FunctionArgs *args, Vector3 *v0in, Vector3 &v1, Vector3 &v2, Vector3 &v3)
{
	Vector3 v0;
	Vector3 v(0.0,0.0,0.0);
	if(v0in == NULL) v0 = v; else v0 = *v0in;
	Real v0x = v0.x;
	Real v0y = v0.y;
	Real v0z = v0.z;
	Real v1x = v1.x;
	Real v1y = v1.y;
	Real v1z = v1.z;
	Real v2x = v2.x;
	Real v2y = v2.y;
	Real v2z = v2.z;
	Real v3x = v3.x;
	Real v3y = v3.y;
	Real v3z = v3.z;
	T q;
	v.x = v0x + 0.25 * v1x + 0.25 * v2x + 0.25 * v3x;
	v.y = v0y + 0.25 * v1y + 0.25 * v2y + 0.25 * v3y;
	v.z = v0z + 0.25 * v1z + 0.25 * v2z + 0.25 * v3z;
	f(&v,&q,args);
	return q;
}

template<class T, class F>
T tetquad4(F f, FunctionArgs *args, Vector3 *v0in, Vector3 &v1, Vector3 &v2, Vector3 &v3)
{
	Vector3 v0;
	Vector3 v(0.0,0.0,0.0);
	if(v0in == NULL) v0 = v; else v0 = *v0in;
	Real v0x = v0.x;
	Real v0y = v0.y;
	Real v0z = v0.z;
	Real v1x = v1.x;
	Real v1y = v1.y;
	Real v1z = v1.z;
	Real v2x = v2.x;
	Real v2y = v2.y;
	Real v2z = v2.z;
	Real v3x = v3.x;
	Real v3y = v3.y;
	Real v3z = v3.z;
	T q;
	T r;
	v.x = v0x + 0.585410196624969 * v1x + 0.138196601125011 * v2x + 0.138196601125011 * v3x;
	v.y = v0y + 0.585410196624969 * v1y + 0.138196601125011 * v2y + 0.138196601125011 * v3y;
	v.z = v0z + 0.585410196624969 * v1z + 0.138196601125011 * v2z + 0.138196601125011 * v3z;
	f(&v,&q,args);
	q *= 0.25;
	v.x = v0x + 0.138196601125011 * v1x + 0.138196601125011 * v2x + 0.138196601125011 * v3x;
	v.y = v0y + 0.138196601125011 * v1y + 0.138196601125011 * v2y + 0.138196601125011 * v3y;
	v.z = v0z + 0.138196601125011 * v1z + 0.138196601125011 * v2z + 0.138196601125011 * v3z;
	f(&v,&r,args);
	r *= 0.25;
	q += r;
	v.x = v0x + 0.138196601125011 * v1x + 0.138196601125011 * v2x + 0.585410196624969 * v3x;
	v.y = v0y + 0.138196601125011 * v1y + 0.138196601125011 * v2y + 0.585410196624969 * v3y;
	v.z = v0z + 0.138196601125011 * v1z + 0.138196601125011 * v2z + 0.585410196624969 * v3z;
	f(&v,&r,args);
	r *= 0.25;
	q += r;
	v.x = v0x + 0.138196601125011 * v1x + 0.585410196624969 * v2x + 0.138196601125011 * v3x;
	v.y = v0y + 0.138196601125011 * v1y + 0.585410196624969 * v2y + 0.138196601125011 * v3y;
	v.z = v0z + 0.138196601125011 * v1z + 0.585410196624969 * v2z + 0.138196601125011 * v3z;
	f(&v,&r,args);
	r *= 0.25;
	q += r;
	return q;
}

template<class T, class F>
T tetquad11(F f, FunctionArgs *args, Vector3 *v0in, Vector3 &v1, Vector3 &v2, Vector3 &v3)
{
	Vector3 v0;
	Vector3 v(0.0,0.0,0.0);
	if(v0in == NULL) v0 = v; else v0 = *v0in;
	Real v0x = v0.x;
	Real v0y = v0.y;
	Real v0z = v0.z;
	Real v1x = v1.x;
	Real v1y = v1.y;
	Real v1z = v1.z;
	Real v2x = v2.x;
	Real v2y = v2.y;
	Real v2z = v2.z;
	Real v3x = v3.x;
	Real v3y = v3.y;
	Real v3z = v3.z;
	T q;
	T r;
	v.x = v0x + 0.25 * v1x + 0.25 * v2x + 0.25 * v3x;
	v.y = v0y + 0.25 * v1y + 0.25 * v2y + 0.25 * v3y;
	v.z = v0z + 0.25 * v1z + 0.25 * v2z + 0.25 * v3z;
	f(&v,&q,args);
	q *= -0.0789333333333333;
	v.x = v0x + 0.785714285714286 * v1x + 0.0714285714285714 * v2x + 0.0714285714285714 * v3x;
	v.y = v0y + 0.785714285714286 * v1y + 0.0714285714285714 * v2y + 0.0714285714285714 * v3y;
	v.z = v0z + 0.785714285714286 * v1z + 0.0714285714285714 * v2z + 0.0714285714285714 * v3z;
	f(&v,&r,args);
	r *= 0.0457333333333333;
	q += r;
	v.x = v0x + 0.0714285714285714 * v1x + 0.0714285714285714 * v2x + 0.0714285714285714 * v3x;
	v.y = v0y + 0.0714285714285714 * v1y + 0.0714285714285714 * v2y + 0.0714285714285714 * v3y;
	v.z = v0z + 0.0714285714285714 * v1z + 0.0714285714285714 * v2z + 0.0714285714285714 * v3z;
	f(&v,&r,args);
	r *= 0.0457333333333333;
	q += r;
	v.x = v0x + 0.0714285714285714 * v1x + 0.0714285714285714 * v2x + 0.785714285714286 * v3x;
	v.y = v0y + 0.0714285714285714 * v1y + 0.0714285714285714 * v2y + 0.785714285714286 * v3y;
	v.z = v0z + 0.0714285714285714 * v1z + 0.0714285714285714 * v2z + 0.785714285714286 * v3z;
	f(&v,&r,args);
	r *= 0.0457333333333333;
	q += r;
	v.x = v0x + 0.0714285714285714 * v1x + 0.785714285714286 * v2x + 0.0714285714285714 * v3x;
	v.y = v0y + 0.0714285714285714 * v1y + 0.785714285714286 * v2y + 0.0714285714285714 * v3y;
	v.z = v0z + 0.0714285714285714 * v1z + 0.785714285714286 * v2z + 0.0714285714285714 * v3z;
	f(&v,&r,args);
	r *= 0.0457333333333333;
	q += r;
	v.x = v0x + 0.100596423833201 * v1x + 0.399403576166799 * v2x + 0.399403576166799 * v3x;
	v.y = v0y + 0.100596423833201 * v1y + 0.399403576166799 * v2y + 0.399403576166799 * v3y;
	v.z = v0z + 0.100596423833201 * v1z + 0.399403576166799 * v2z + 0.399403576166799 * v3z;
	f(&v,&r,args);
	r *= 0.149333333333333;
	q += r;
	v.x = v0x + 0.399403576166799 * v1x + 0.100596423833201 * v2x + 0.399403576166799 * v3x;
	v.y = v0y + 0.399403576166799 * v1y + 0.100596423833201 * v2y + 0.399403576166799 * v3y;
	v.z = v0z + 0.399403576166799 * v1z + 0.100596423833201 * v2z + 0.399403576166799 * v3z;
	f(&v,&r,args);
	r *= 0.149333333333333;
	q += r;
	v.x = v0x + 0.399403576166799 * v1x + 0.399403576166799 * v2x + 0.100596423833201 * v3x;
	v.y = v0y + 0.399403576166799 * v1y + 0.399403576166799 * v2y + 0.100596423833201 * v3y;
	v.z = v0z + 0.399403576166799 * v1z + 0.399403576166799 * v2z + 0.100596423833201 * v3z;
	f(&v,&r,args);
	r *= 0.149333333333333;
	q += r;
	v.x = v0x + 0.399403576166799 * v1x + 0.100596423833201 * v2x + 0.100596423833201 * v3x;
	v.y = v0y + 0.399403576166799 * v1y + 0.100596423833201 * v2y + 0.100596423833201 * v3y;
	v.z = v0z + 0.399403576166799 * v1z + 0.100596423833201 * v2z + 0.100596423833201 * v3z;
	f(&v,&r,args);
	r *= 0.149333333333333;
	q += r;
	v.x = v0x + 0.100596423833201 * v1x + 0.399403576166799 * v2x + 0.100596423833201 * v3x;
	v.y = v0y + 0.100596423833201 * v1y + 0.399403576166799 * v2y + 0.100596423833201 * v3y;
	v.z = v0z + 0.100596423833201 * v1z + 0.399403576166799 * v2z + 0.100596423833201 * v3z;
	f(&v,&r,args);
	r *= 0.149333333333333;
	q += r;
	v.x = v0x + 0.100596423833201 * v1x + 0.100596423833201 * v2x + 0.399403576166799 * v3x;
	v.y = v0y + 0.100596423833201 * v1y + 0.100596423833201 * v2y + 0.399403576166799 * v3y;
	v.z = v0z + 0.100596423833201 * v1z + 0.100596423833201 * v2z + 0.399403576166799 * v3z;
	f(&v,&r,args);
	r *= 0.149333333333333;
	q += r;
	return q;
}

#endif
