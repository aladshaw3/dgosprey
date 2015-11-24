#include "GAnisotropicDiffusion.h"

template<>
InputParameters validParams<GAnisotropicDiffusion>()
{
	InputParameters params = validParams<Kernel>();
	params.addParam<Real>("Dxx",0,"xx-component of diffusion tensor");
	params.addParam<Real>("Dxy",0,"xy-component of diffusion tensor");
	params.addParam<Real>("Dxz",0,"xz-component of diffusion tensor");
	params.addParam<Real>("Dyx",0,"yx-component of diffusion tensor");
	params.addParam<Real>("Dyy",0,"yy-component of diffusion tensor");
	params.addParam<Real>("Dyz",0,"yz-component of diffusion tensor");
	params.addParam<Real>("Dzx",0,"zx-component of diffusion tensor");
	params.addParam<Real>("Dzy",0,"zy-component of diffusion tensor");
	params.addParam<Real>("Dzz",0,"zz-component of diffusion tensor");
	return params;
}

GAnisotropicDiffusion::GAnisotropicDiffusion(const InputParameters & parameters) :
Kernel(parameters),
_Dxx(getParam<Real>("Dxx")),
_Dxy(getParam<Real>("Dxy")),
_Dxz(getParam<Real>("Dxz")),
_Dyx(getParam<Real>("Dyx")),
_Dyy(getParam<Real>("Dyy")),
_Dyz(getParam<Real>("Dyz")),
_Dzx(getParam<Real>("Dzx")),
_Dzy(getParam<Real>("Dzy")),
_Dzz(getParam<Real>("Dzz"))
{
	_Diffusion(0,0) = _Dxx;
	_Diffusion(0,1) = _Dxy;
	_Diffusion(0,2) = _Dxz;
	
	_Diffusion(1,0) = _Dyx;
	_Diffusion(1,1) = _Dyy;
	_Diffusion(1,2) = _Dyz;
	
	_Diffusion(2,0) = _Dzx;
	_Diffusion(2,1) = _Dzy;
	_Diffusion(2,2) = _Dzz;
}

Real GAnisotropicDiffusion::computeQpResidual()
{
	return _Diffusion*_grad_test[_i][_qp]*_grad_u[_qp];
}

Real GAnisotropicDiffusion::computeQpJacobian()
{
	return _Diffusion*_grad_test[_i][_qp]*_grad_phi[_j][_qp];
}
