//----------------------------------------
//  Created by Austin Ladshaw on 09/14/2015
//  Copyright (c) 2015
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

#ifndef AdsorptionHeatAccumulation_H
#define AdsorptionHeatAccumulation_H

#include "Kernel.h"

//Forward Declarations
class AdsorptionHeatAccumulation;

template<>
InputParameters validParams<AdsorptionHeatAccumulation>();

class AdsorptionHeatAccumulation : public Kernel
{
public:
	
	AdsorptionHeatAccumulation(const InputParameters & parameters);
	
protected:
	virtual Real computeQpResidual();
	virtual Real computeQpJacobian();
	
private:
	const MaterialProperty<Real> & _porosity;
	const MaterialProperty<Real> & _pellet_density;
	std::vector<VariableValue *> _solid_heat;
	std::vector<VariableValue *> _solid_heat_old;
};

#endif