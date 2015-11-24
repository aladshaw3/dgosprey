#include "BedMassAccumulation.h"

template<>
InputParameters validParams<BedMassAccumulation>()
{
  InputParameters params = validParams<TimeDerivative>();
  params.addParam<int>("index", 0, "The index of the coupling variable. Must be given in same order of appearance as in the FlowProperties Material block. Indexing starts from 0. 0 is default value.");
  return params;
}


BedMassAccumulation::BedMassAccumulation(const InputParameters & parameters)
:TimeDerivative(parameters),
_index(getParam<int>("index")),
_retardation(getMaterialProperty<std::vector<Real> >("retardation"))
{

}

Real BedMassAccumulation::computeQpResidual()
{
  return _retardation[_qp][_index] * TimeDerivative::computeQpResidual();
}

Real BedMassAccumulation::computeQpJacobian()
{
  return _retardation[_qp][_index] * TimeDerivative::computeQpJacobian();
}
