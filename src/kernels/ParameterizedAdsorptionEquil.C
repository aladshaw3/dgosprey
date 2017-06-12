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
    return params;
}

ParameterizedAdsorptionEquil::ParameterizedAdsorptionEquil(const InputParameters & parameters)
: CoupledLangmuirForcingFunction(parameters),
_index(getParam<unsigned int>("index")),
_magpie_dat(getMaterialProperty< MAGPIE_DATA >("magpie_data"))
{
    
}

Real ParameterizedAdsorptionEquil::computeQpResidual()
{
    double q_est = 0.0;
    MAGPIE_DATA magpie_copy;
    magpie_copy = _magpie_dat[_qp];
    
    //Call MAGPIE Simulation
    if (_magpie_dat[_qp].gsta_dat[_index].qmax > 0.0)
    {
        int success = 0;
        success = MAGPIE( (void *)&magpie_copy );
        if (success < 0 || success > 3)
        {
            for (int i=0; i<magpie_copy.sys_dat.N; i++)
            {
                if (i != _index)
                {
                    magpie_copy.gpast_dat[i].y = 0.0;
                    magpie_copy.sys_dat.Carrier = true;
                }
            }
            success = MAGPIE( (void *)&magpie_copy );
            if (success < 0 || success > 3)
            {
                mError(simulation_fail);
            }
        }
        else success = 0;
        
        q_est = magpie_copy.gpast_dat[_index].q;
        
    }
    else
    {
        return 0.0;
    }
    
    _maxcap = _magpie_dat[_qp].gsta_dat[_index].qmax;
    if (magpie_copy.gpast_dat[_index].y == 0.0)
    {
        double K = dq_dp (0.0, (void *)&magpie_copy, _index)*8.3144621*_magpie_dat[_qp].sys_dat.T;
        _langmuircoef = K/_maxcap;
    }
    else
    {
        double pi = _magpie_dat[_qp].sys_dat.PT*magpie_copy.gpast_dat[_index].y;
        double ci = pi/8.3144621/_magpie_dat[_qp].sys_dat.T;
        _langmuircoef = q_est/(ci*(_maxcap+1e-6-q_est));
    }
    
    return CoupledLangmuirForcingFunction::computeQpResidual();
}

Real ParameterizedAdsorptionEquil::computeQpJacobian()
{
    return CoupledLangmuirForcingFunction::computeQpJacobian();
}

Real ParameterizedAdsorptionEquil::computeQpOffDiagJacobian(unsigned int jvar)
{
   
    if (jvar == _coupled_var)
    {
        double q_est = 0.0;
        MAGPIE_DATA magpie_copy;
        magpie_copy = _magpie_dat[_qp];
        
        //Call MAGPIE Simulation
        if (_magpie_dat[_qp].gsta_dat[_index].qmax > 0.0)
        {
            int success = 0;
            success = MAGPIE( (void *)&magpie_copy );
            if (success < 0 || success > 3)
            {
                for (int i=0; i<magpie_copy.sys_dat.N; i++)
                {
                    if (i != _index)
                    {
                        magpie_copy.gpast_dat[i].y = 0.0;
                        magpie_copy.sys_dat.Carrier = true;
                    }
                }
                success = MAGPIE( (void *)&magpie_copy );
                if (success < 0 || success > 3)
                {
                    mError(simulation_fail);
                }
            }
            else success = 0;
            
            q_est = magpie_copy.gpast_dat[_index].q;
            
        }
        else
        {
            return 0.0;
        }
        
        _maxcap = _magpie_dat[_qp].gsta_dat[_index].qmax;
        if (magpie_copy.gpast_dat[_index].y == 0.0)
        {
            double K = dq_dp (0.0, (void *)&magpie_copy, _index)*8.3144621*_magpie_dat[_qp].sys_dat.T;
            _langmuircoef = K/_maxcap;
        }
        else
        {
            double pi = _magpie_dat[_qp].sys_dat.PT*magpie_copy.gpast_dat[_index].y;
            double ci = pi/8.3144621/_magpie_dat[_qp].sys_dat.T;
            _langmuircoef = q_est/(ci*(_maxcap+1e-6-q_est));
        }
        
        return CoupledLangmuirForcingFunction::computeQpOffDiagJacobian(jvar);
    }
    
    return 0.0;
}


