 [GlobalParams]

[]


[Mesh]
type = GeneratedMesh
dim = 1
nx = 10
xmax = 1.0
[]

[Variables]

[./u]
order = FIRST
family = LAGRANGE
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

[./u_adv]
type = Convection
variable = u
x=2.0
y=0.0
z=0.0
[../]

[./u_diff]
type = CoefDiffusion
variable = u
coef=0.01
[../]

[]


[AuxKernels]

[]

[BCs]

[./u_bc_bot]
type = DirichletBC
variable = u
boundary = 0
value = 1.0
[../]

[./u_bc_top]
type = NeumannBC
variable = u
boundary = 1
value = 0.0
[../]

[]

[Materials]

[]

[Postprocessors]

[./u_out]
type = AverageNodalVariableValue
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
