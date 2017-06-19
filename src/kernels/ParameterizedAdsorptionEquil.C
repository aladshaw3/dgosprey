/*!
 *  \file ParameterizedAdsorptionEquil.C
 *	\brief Standard kernel for  simulating multispecies adsorption equilibria
 *	\details This file creates a standard MOOSE kernel for the call the MAGPIE adsorption equilibria function
 *
 *
 *  \author Austin Ladshaw, Alexander Wiechert
 *	\date 06/12/2017
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


#include "ParameterizedAdsorptionEquil.h"

template<>
InputParameters validParams<ParameterizedAdsorptionEquil>()
{
    InputParameters params = validParams<CoupledLangmuirForcingFunction>();
    params.addParam<unsigned int>("index",0,"Index of the species that we are interested in.");
	params.addRequiredCoupledVar("ads_est","Name of the AuxVariable estimating adsorption");
    return params;
}

ParameterizedAdsorptionEquil::ParameterizedAdsorptionEquil(const InputParameters & parameters)
: CoupledLangmuirForcingFunction(parameters),
_index(getParam<unsigned int>("index")),
_magpie_dat(getMaterialProperty< MAGPIE_DATA >("magpie_data")),
_q_est(coupledValue("ads_est"))
{
    
}

Real ParameterizedAdsorptionEquil::computeQpResidual()
{
	if (_magpie_dat[_qp].gsta_dat[_index].qmax == 0.0)
		return 0.0;
	
    _maxcap = _magpie_dat[_qp].gsta_dat[_index].qmax;
	/*
    if (_magpie_dat[_qp].gpast_dat[_index].y == 0.0 || _q_est[_qp] == 0.0)
    {
        double K = dq_dp (0.0, (void *)&_magpie_dat[_qp], _index)*8.3144621*_magpie_dat[_qp].sys_dat.T;
        _langmuircoef = K/_maxcap;
    }
    else
    {
        double pi = _magpie_dat[_qp].sys_dat.PT*_magpie_dat[_qp].gpast_dat[_index].y;
        double ci = pi/8.3144621/_magpie_dat[_qp].sys_dat.T;
        _langmuircoef = fabs(_q_est[_qp]/(ci*(_maxcap+1e-6-_q_est[_qp])));
    }
	 */
	/*
	double m,T,p,c;
	T = _magpie_dat[_qp].sys_dat.T;
	p = _magpie_dat[_qp].sys_dat.PT*_magpie_dat[_qp].gpast_dat[_index].y;
	c = p/R/T;
	m = dq_dp(p, (void *)&_magpie_dat[_qp], _index);
	double Ak,B,C;
	Ak = m*R*T*c*c;
	B = (2.0*m*R*T*c) - _maxcap;
	C = m*R*T;
	if (m == 0.0)
		_langmuircoef = 0.0;
	if (c == 0.0)
		_langmuircoef = (m*R*T)/_maxcap;
	else
	{
		double root = std::pow(((B*B) - (4.0*Ak*C)),0.5);
		if ((-B+root) >= 0.0)
			_langmuircoef = (-B+root)/(2.0*Ak);
		else
			_langmuircoef = (-B-root)/(2.0*Ak);
	}
	_langmuircoef = (m*R*T)/_maxcap;
	 */
	double m,T,p,theta;
	T = _magpie_dat[_qp].sys_dat.T;
	p = _magpie_dat[_qp].sys_dat.PT*_magpie_dat[_qp].gpast_dat[_index].y;
	m = dq_dp(p, (void *)&_magpie_dat[_qp], _index);
	theta = _q_est[_qp]/_maxcap;
	_langmuircoef = (((R*T)/_maxcap)*m)/((1.0-theta)*(1.0-theta));

	//std::cout << _index << "\t" << _langmuircoef << std::endl;
    return CoupledLangmuirForcingFunction::computeQpResidual();
}

Real ParameterizedAdsorptionEquil::computeQpJacobian()
{
	if (_magpie_dat[_qp].gsta_dat[_index].qmax == 0.0)
		return 0.0;
    return CoupledLangmuirForcingFunction::computeQpJacobian();
}

Real ParameterizedAdsorptionEquil::computeQpOffDiagJacobian(unsigned int jvar)
{
	if (_magpie_dat[_qp].gsta_dat[_index].qmax == 0.0)
		return 0.0;
	
    if (jvar == _coupled_var)
    {
        _maxcap = _magpie_dat[_qp].gsta_dat[_index].qmax;
		/*
        if (_magpie_dat[_qp].gpast_dat[_index].y == 0.0 || _q_est[_qp] == 0.0)
        {
            double K = dq_dp (0.0, (void *)&_magpie_dat[_qp], _index)*8.3144621*_magpie_dat[_qp].sys_dat.T;
            _langmuircoef = K/_maxcap;
        }
        else
        {
            double pi = _magpie_dat[_qp].sys_dat.PT*_magpie_dat[_qp].gpast_dat[_index].y;
            double ci = pi/8.3144621/_magpie_dat[_qp].sys_dat.T;
            _langmuircoef = fabs(_q_est[_qp]/(ci*(_maxcap+1e-6-_q_est[_qp])));
        }
        */
		/*
		double m,T,p,c;
		T = _magpie_dat[_qp].sys_dat.T;
		p = _magpie_dat[_qp].sys_dat.PT*_magpie_dat[_qp].gpast_dat[_index].y;
		c = p/R/T;
		m = dq_dp(p, (void *)&_magpie_dat[_qp], _index);
		double Ak,B,C;
		Ak = m*R*T*c*c;
		B = (2.0*m*R*T*c) - _maxcap;
		C = m*R*T;
		if (m == 0.0)
			_langmuircoef = 0.0;
		if (c == 0.0)
			_langmuircoef = (m*R*T)/_maxcap;
		else
		{
			double root = std::pow(((B*B) - (4.0*Ak*C)),0.5);
			if ((-B+root) >= 0.0)
				_langmuircoef = (-B+root)/(2.0*Ak);
			else
				_langmuircoef = (-B-root)/(2.0*Ak);
		}
		_langmuircoef = (m*R*T)/_maxcap;
		 */
		double m,T,p,theta;
		T = _magpie_dat[_qp].sys_dat.T;
		p = _magpie_dat[_qp].sys_dat.PT*_magpie_dat[_qp].gpast_dat[_index].y;
		m = dq_dp(p, (void *)&_magpie_dat[_qp], _index);
		theta = _q_est[_qp]/_maxcap;
		_langmuircoef = (((R*T)/_maxcap)*m)/((1.0-theta)*(1.0-theta));
		
        return CoupledLangmuirForcingFunction::computeQpOffDiagJacobian(jvar);
    }
	
    return 0.0;
}


