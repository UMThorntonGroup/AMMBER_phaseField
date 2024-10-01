#ifndef INCLUDE_VARTYPEENUMS_H_
#define INCLUDE_VARTYPEENUMS_H_

enum fieldType
{
  SCALAR,
  VECTOR
};

enum PDEType
{
  EXPLICIT_TIME_DEPENDENT,
  IMPLICIT_TIME_DEPENDENT,
  TIME_INDEPENDENT,
  AUXILIARY
};

#endif
