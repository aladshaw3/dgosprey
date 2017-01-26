/*!
 *  \file AdsorptionMassTransfer.h
 *	\brief Standard kernel for the transfer of mass via adsorption
 *  \author Austin Ladshaw
 *	\date 01/29/2015
 *	\copyright This kernel was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science and was developed for use
 *			   by Idaho National Laboratory and Oak Ridge National Laboratory
 *			   engineers and scientists. Portions Copyright (c) 2016, all
 *             rights reserved.
 *
 *			   Austin Ladshaw does not claim any ownership or copyright to the
 *			   MOOSE framework in which these kernels are constructed, only
 *			   the kernels themselves. The MOOSE framework copyright is held
 *			   by the Battelle Energy Alliance, LLC (c) 2010, all rights reserved.
 */

/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/


#include "AdsorptionMassTransfer.h"

template<>
InputParameters validParams<AdsorptionMassTransfer>()
{
	InputParameters params = validParams<Kernel>();
	params.addCoupledVar("solid_conc","Coupled variables for the solid adsorption concentrations");
	params.addCoupledVar("solid_pert","Coupled variables for the solid adsorption concentrations");
	return params;
}

AdsorptionMassTransfer::AdsorptionMassTransfer(const InputParameters & parameters) :
Kernel(parameters),
_porosity(getMaterialProperty<Real>("porosity")),
_pellet_density(getMaterialProperty<Real>("pellet_density")),
_solid(coupledValue("solid_conc")),
_solid_old(coupledValueOld("solid_conc"))
//_solid_old(coupledValueOlder("solid_conc"))
//_solid_old(coupledValue("solid_pert"))
{

}

Real AdsorptionMassTransfer::computeQpResidual()
{
	//std::cout << _u[_qp] << "\t" << _solid_old[_qp] << "\t" << _solid[_qp] << "\t" << _dt_old << std::endl;
	//return 0.0;
	//return -(1.0-_porosity[_qp])*_pellet_density[_qp]*(_solid[_qp]-_solid_old[_qp])*0.125*_dt_old*_test[_i][_qp];
	return -(1.0-_porosity[_qp])*_pellet_density[_qp]*((_solid[_qp]-_solid_old[_qp])/(1.0*_dt_old))*_test[_i][_qp];
	//return -(1.0-_porosity[_qp])*_pellet_density[_qp]*( (_solid_old[_qp] - _solid[_qp]) / sqrt(DBL_EPSILON) ) * _u_dot[_qp]*_test[_i][_qp];
}

Real AdsorptionMassTransfer::computeQpJacobian()
{
	return 0.0;
	//return -(1.0-_porosity[_qp])*_pellet_density[_qp]*fabs( (_solid_old[_qp] - _solid[_qp]) / sqrt(DBL_EPSILON) ) * _test[_i][_qp]*_phi[_j][_qp]*_du_dot_du[_qp];
}

