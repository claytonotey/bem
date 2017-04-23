#include "Interaction.h"

// FIXME: limits the number of triangles to 2^16=65536
Key getKey(Key t, Key s)
{
  return (t<<16)+s;
}

// Triangle-Point
TrianglePointInteraction :: TrianglePointInteraction(Vector3& v0_, Vector3& v1_, Vector3& v2_, int n, int useTotal) 
{
  v0 = v0_;
  v1 = v1_;
  v2 = v2_;
  this->n = n;
  this->useTotal = useTotal;
  scalar.resize(n);
  gradient.resize(n);
  bSet = false;
  useCount = 0;
}

TrianglePointInteraction :: ~TrianglePointInteraction()
{
}

void TrianglePointInteraction :: start()
{
  i = 0;
}

void TrianglePointInteraction :: end()
{
  useCount++;
  bSet = true;
}

void TrianglePointInteraction :: next()
{
  i++;
}

void TrianglePointInteraction :: getScalar(Complex& z)
{
  z = scalar[i];
}

void TrianglePointInteraction :: getGradient(ComplexVector3& v)
{
  v = gradient[i];
}

void TrianglePointInteraction :: setScalar(Complex &z)
{
  scalar[i] = z;
}

void TrianglePointInteraction :: setGradient(ComplexVector3 &v)
{
  gradient[i] = v;
}


// Triangle-Triangle
TriangleTriangleInteraction :: TriangleTriangleInteraction() {}
TriangleTriangleInteraction :: TriangleTriangleInteraction(Vector3& v0t_, Vector3& v1t_, Vector3& v2t_, Vector3& v0s_, Vector3& v1s_, Vector3& v2s_, int nt, int ns, int useTotal)
{
  v0t = v0t_;
  v1t = v1t_;
  v2t = v2t_;
  v0s = v0s_;
  v1s = v1s_;
  v2s = v2s_;
  this->nt = nt;
  this->ns = ns;
  this->useTotal = useTotal;
  scalar.resize(nt*ns);
  gradient.resize(nt*ns);
  bSet = false;
  useCount = 0;
}

TriangleTriangleInteraction :: ~TriangleTriangleInteraction()
{
}

void TriangleTriangleInteraction :: start()
{
  i = 0;
}

void TriangleTriangleInteraction :: end()
{
  useCount++;
  bSet = true;
}

void TriangleTriangleInteraction :: next()
{
  i++;
}

void TriangleTriangleInteraction :: getScalar(Complex &z)
{
  z = scalar[i];
}

void TriangleTriangleInteraction :: getGradient(ComplexVector3 &v)
{
  v = gradient[i];
}

void TriangleTriangleInteraction :: setScalar(Complex &z)
{
  scalar[i] = z;
}

void TriangleTriangleInteraction :: setGradient(ComplexVector3 &v)
{
  gradient[i] = v;
}
