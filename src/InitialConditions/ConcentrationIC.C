#include "ConcentrationIC.h"

template<> InputParameters validParams<ConcentrationIC>()
{
  InputParameters params = validParams<InitialCondition>();
  params.addRequiredParam<Real>("initial_mole_frac","The initial mole fraction of the gas species in the column. Note: All mole fractions must sum to one, otherwise errors will occur that will not be caught!!!");
  params.addRequiredParam<Real>("initial_press","Initial condition for total pressure");
  params.addRequiredParam<Real>("initial_temp","Initial condition for temperature");
  return params;
}

ConcentrationIC::ConcentrationIC(const InputParameters & parameters) :
InitialCondition(parameters),
_y_IC(getParam<Real>("initial_mole_frac")),
_PT_IC(getParam<Real>("initial_press")),
_T_IC(getParam<Real>("initial_temp"))
{
  
}

Real ConcentrationIC::value(const Point & p)
{
  //Note: we will not check to ensure the mole fraction is valid. This may cause error.
  return (_PT_IC * _y_IC) / (8.3144621 * _T_IC);
}