#ifndef BEDMASSACCUMULATION_H
#define BEDMASSACCUMULATION_H

#include "TimeDerivative.h"

//Forward Declarations
class BedMassAccumulation;

template<>
InputParameters validParams<BedMassAccumulation>();

class BedMassAccumulation : public TimeDerivative
{
public:
  
  BedMassAccumulation(const InputParameters & parameters);
  
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  
private:
  int _index;
  const MaterialProperty<std::vector<Real> > & _retardation;
};

#endif //BEDMASSACCUMULATION_H
