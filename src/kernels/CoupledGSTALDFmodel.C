/*!
 *  \file CoupledGSTLDFAmodel.h
 *	\brief Standard kernel for coupling non-linear variables via the GSTA model with LDF kinetics
 *	\details This file creates a standard MOOSE kernel for the coupling of non-linear variables
 *			together via the GSTA model, and applies linear driving force kinetics for the rate
 *			of adsorption.
 *
 *			This kernel extends the CoupledGSTAmodel kernel by calculating the model parameters from
 *			information in the ThermodynamicProperties material. In addition, it extends the LinearDrivingForce
 *			kernel by estimating the overall LDF rate coefficient from material properties.
 *			The Kno parameters (described below) are to be estimated from the site enthalpies (dHno)
 *			and entropies (dSno) using the van't Hoff expression (shown below). Thus, this model
 *			is inherently a function of temperature and will require a different form of coupling
 *			with the temperature parameter.
 *
 *			In addition, the linear driving force parameter (k) is estimated using the Resistance-in-Series
 *			model, which couples film mass transfer (kf), pore diffusion (Dp), and surface diffusion (Dc)
 *			into a single lumped rate parameter (k). Those parameters all come from material properties
 *			files in the DGOSPREY framework.
 *
 *			Resistance-in-series: (1/k) = (rhop*q*rp/(3*kf*C)) + (rhop*q*rp*rp/(15*ep*Dp)) + (rc*rc/(15*Dc))
 *			where rhop is the particle density, q is the adsorption, rp is the particle radius, kf is the film
 *			mass transfer parameter, C is the concentration in the gas phase, ep is the particle porosity,
 *			Dp is the pore diffusion parameter, rc is the adsorbent crystal radius, and Dc is the surface
 *			diffusion parameter.
 *
 *			van't Hoff: ln(Kno) = -dHno/(R*T) + dSno/R
 *			where R is the gas law constant and T is the column temperature.
 *
 *			GSTA isotherm: q = (q_max / m) * SUM(n*Kno*(p/Po)^n)/(1+SUM(Kno*(p/Po)^n))
 *			where q is amount adsorbed, q_max is the maximum capacity, m is the number of adsorption sites
 *			and Kno are the dimensionless equilibrium parameters. Also, p is partial pressure of gas and Po
 *			is taken as a reference state pressure (100 kPa).
 *
 *	\note	For the use of this kernel, our coupled variable with be a gas concentration in mol/L (C), therefore,
 *			we need to use ideal gas law to rewrite the GSTA model in terms of C as opposed to p. Thus, we are also
 *			forced to couple with kernel with column temperature.
 *
 *			Ideal Gas Law: p = C*R*T
 *
 *  \author Austin Ladshaw, Alexander Wiechert
 *	\date 08/28/2017
 *	\copyright This kernel was designed and built at the Georgia Institute
 *             of Technology by Alexander Wiechert for PhD research in the area
 *             of adsorption and surface science and was developed for use
 *			   by Idaho National Laboratory and Oak Ridge National Laboratory
 *			   engineers and scientists. Portions Copyright (c) 2017, all
 *             rights reserved.
 *
 *			   Alexander Wiechert does not claim any ownership or copyright to the
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

#include "CoupledGSTALDFmodel.h"

template<>
InputParameters validParams<CoupledGSTALDFmodel>()
{
	InputParameters params = validParams<CoupledGSTAmodel>();
	return params;
}

CoupledGSTALDFmodel::CoupledGSTALDFmodel(const InputParameters & parameters)
: CoupledGSTAmodel(parameters),
_coupled_u_old(coupledValueOld("coupled_gas")),
_coupled_temp_old(coupledValueOld("coupled_temp")),
_pellet_density(getMaterialProperty<Real>("pellet_density")),
_pellet_diameter(getMaterialProperty<Real>("pellet_diameter")),
_crystal_radius(getMaterialProperty<Real>("crystal_radius")),
_binder_porosity(getMaterialProperty<Real>("binder_porosity")),
_binder_fraction(getMaterialProperty<Real>("binder_ratio")),
_film_transfer(getMaterialProperty<std::vector<Real> >("film_transfer")),
_pore_diff(getMaterialProperty<std::vector<Real> >("pore_diffusion")),
_surf_diff(getMaterialProperty<std::vector<Real> >("surface_diffusion"))
{
	
}

void CoupledGSTALDFmodel::computeLDFcoeff()
{
	CoupledGSTAmodel::computeGSTAparams();
	computeGSTAequilibriumOld();
	
	double part_coef = 1.0, Co = 100.0 / (8.3144621 * _coupled_temp_old[_qp]);
	double Henry = (_maxcap * _gstaparam[0])/(_numsites*Co);
	if (std::isnan(_ads_equil/_coupled_u_old[_qp]) || std::isinf(_ads_equil/_coupled_u_old[_qp]))
		part_coef = _pellet_density[_qp]*Henry;
	else
		part_coef = _pellet_density[_qp]*_ads_equil/_coupled_u_old[_qp];
	
	double filmres, poreres, surfres;
	filmres = part_coef * _pellet_diameter[_qp] / (6.0*_film_transfer[_qp][_index]);
	if (_binder_fraction[_qp] == 0.0)
		poreres = part_coef * _pellet_diameter[_qp] * _pellet_diameter[_qp] / (60.0*_pore_diff[_qp][_index]*_binder_porosity[_qp]);
	else
		poreres = part_coef * _pellet_diameter[_qp] * _pellet_diameter[_qp] / (60.0*_pore_diff[_qp][_index]*_binder_porosity[_qp]*_binder_fraction[_qp]);
	if (_surf_diff[_qp][_index] == 0.0)
		surfres = 0.0;
	else
		surfres = _crystal_radius[_qp] * _crystal_radius[_qp] / (15.0 * _surf_diff[_qp][_index]);
	
	double k = filmres + poreres + surfres;
	//std::cout << 1.0/k << std::endl;
	_ldf_coeff = (10.0/k)*(_ads_equil/_maxcap) + (1000000.0/k)*(1.0 - (_ads_equil/_maxcap));
	//std::cout << _ldf_coeff << std::endl;
	//_ldf_coeff = 0.1;
	
}

void CoupledGSTALDFmodel::computeGSTAequilibriumOld()
{
	double top = 0.0, bot = 1.0, Co = 100.0 / (8.3144621 * _coupled_temp_old[_qp]);

	for (int n = 0; n<(int)_numsites; n++)
	{
		top = top + ( (double)(n+1) * _gstaparam[n] * std::pow((_coupled_u_old[_qp]/Co),(double)(n+1)) );
		bot = bot + ( _gstaparam[n] * std::pow((_coupled_u_old[_qp]/Co),(double)(n+1)) );
	}
	
	_ads_equil = (_maxcap/_numsites)*(top/bot);
}

Real CoupledGSTALDFmodel::computeQpResidual()
{
	computeLDFcoeff();
	return _ldf_coeff*CoupledGSTAisotherm::computeQpResidual();
}

Real CoupledGSTALDFmodel::computeQpJacobian()
{
	computeLDFcoeff();
	return _test[_i][_qp]*_ldf_coeff*_phi[_j][_qp];
}

Real CoupledGSTALDFmodel::computeQpOffDiagJacobian(unsigned int jvar)
{
	computeLDFcoeff();
	
	// Off-diagonal element for coupled gas
	if (jvar == _coupled_var_u)
	{
		return -_test[_i][_qp]*_ldf_coeff*CoupledGSTAisotherm::computeGSTAconcDerivative();
	}
	
	// Off-diagonal element for coupled temperature
	if (jvar == _coupled_var_temp)
	{
		return -_test[_i][_qp]*_ldf_coeff*CoupledGSTAmodel::computeGSTAtempDerivative();
	}
	
	
	return 0.0;
}


