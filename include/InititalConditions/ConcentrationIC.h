#ifndef CONCENTRATIONIC_H
#define	CONCENTRATIONIC_H

#include "InitialCondition.h"

class ConcentrationIC;

template<> InputParameters validParams<ConcentrationIC>();

class ConcentrationIC : public InitialCondition
{
public:
  ConcentrationIC(const InputParameters & parameters);
  virtual Real value(const Point & p);
  
private:
  Real _y_IC;
  Real _PT_IC;
  Real _T_IC;
};

#endif //CONCENTRATIONIC_H
