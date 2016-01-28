 [GlobalParams]
vx = 2.0 #R-direction (positive leaves outward, negative pushes inward)
vy = 0.0 #Z-direction (positive moves bottom to top, negative moves top to bottom)
vz = 0.0 #Not used in RZ system

Dxx = 0.01 #Radial Diffusion in cylinder
Dxy = 0.0
Dxz = 0.0

Dyx = 0.0
Dyy = 0.0 #Axial Diffusion in cylinder
Dyz = 0.0

Dzx = 0.0
Dzy = 0.0
Dzz = 0.0
[]


[Mesh]
type = GeneratedMesh
dim = 1
nx = 100
#nx = 10
xmax = 1.0
[]

[Variables]

[./u]
order = CONSTANT
family = MONOMIAL
#order = FIRST
#family = L2_LAGRANGE
initial_condition = 0.0
[../]

[]

[AuxVariables]

[]

[ICs]

[]

[Kernels]

[./u_dot]
type = CoefTimeDerivative
variable = u
Coefficient = 1.0
[../]

[./u_gadv]
type = GAdvection
variable = u
[../]

[./u_gdiff]
type = GAnisotropicDiffusion
variable = u
[../]

[]

[]

[DGKernels]

[./u_dgadv]
type = DGAdvection
variable = u
[../]

[./u_dgdiff]
type = DGAnisotropicDiffusion
variable = u
[../]

[]


[AuxKernels]

[]

[BCs]

[./u_bc_bot]
type = DGFluxLimitedBC
variable = u
boundary = 0
u_input = 1.0
[../]

[./u_bc_top]
type = DGFluxLimitedBC
variable = u
boundary = 1
[../]

[]

[Materials]

[]

[Postprocessors]

[./u_out]
type = SideAverageValue
boundary = 1
variable = u
[../]

[]

[Executioner]

type = Transient
#scheme = implicit-euler
scheme = bdf2

nl_rel_tol = 1e-6
nl_abs_tol = 1e-4
nl_rel_step_tol = 1e-10
nl_abs_step_tol = 1e-10
l_tol = 1e-6
l_max_its = 100

solve_type = PJFNK
line_search = bt
start_time = 0.0
end_time = 1.0
petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
petsc_options_value = 'hypre boomeramg 100'

[./TimeStepper]
#type = ConstantDT
type = SolutionTimeAdaptiveDT
dt = 0.01
[../]

[]

[]

[Outputs]
exodus = true
csv = true
print_linear_residuals = true
[]
