# NTEC representative cell: 2-D liquid shear layer over a contact-electrified
# catalyst interface (Tang et al., ACS Nano 2022 mechanism hypothesis).
#
#   y = film_thickness_m      sheared, CH4-saturated top boundary moving at U
#   0 < y < film_thickness_m  liquid (block 1): incompressible Navier-Stokes + CH4
#   y = 0                     catalyst interface: contact charge sigma, CH4 sink
#   -d < y < 0                dielectric/cocatalyst layer (block 0): Poisson
#   y = -d                    grounded back contact
#   x                         periodic (fully developed layer, length 4 h)
#
# Governing equations (steady, isothermal):
#   rho (u.grad) u = -grad p + mu lap u,  div u = 0           (liquid)
#   div(D grad c) = 0  (u.grad c vanishes for the fully developed layer)
#   div(eps grad phi) = 0,  eps dphi/dn = sigma at the interface
#   D dc/dn = -k_s c at the interface (first-order surface consumption)
#
# Exchanged with the runner's Cantera step each coupling iteration:
#   interface_charge_C_m2       sigma(wall shear rate), contact-electrification closure
#   interface_rate_constant_m_s k_s = r_Cantera(c_s) exp(min(beta|E|, cap)/kT) / c_s
# Iteration 1 uses the case guesses. Every top-level variable is overridden
# from hydrogen_case.json; defaults only make the file independently checkable.

# --- geometry / operating / properties (hydrogen_case.json) ---
reactor_volume_m3 = 1e-3
interface_area_m2 = 0.1
temperature_K = 373.15
pressure_Pa = 101325
methane_mass_flow_kg_s = 1e-7
shear_rate_s_inv = 100
mechanical_power_W_kg = 20
liquid_viscosity_Pa_s = 2e-3
permittivity_F_m = 7e-11
ionic_conductivity_S_m = 1

# --- model_parameters (hydrogen_case.json) ---
liquid_density_kg_m3 = 6100
film_thickness_m = 1e-3
dielectric_thickness_m = 1e-5
methane_diffusivity_m2_s = 2e-9
methane_henry_mol_m3_Pa = 1e-5
interface_charge_guess_C_m2 = 0
interface_rate_constant_guess_m_s = 1e-8
mesh_elements_x = 4
mesh_elements_film = 40
mesh_elements_dielectric = 10
coupling_iteration = 1

# --- exchanged with Cantera (runner override from iteration 2) ---
interface_charge_C_m2 = ${interface_charge_guess_C_m2}
interface_rate_constant_m_s = ${interface_rate_constant_guess_m_s}

# --- derived ---
wall_speed = ${fparse shear_rate_s_inv * film_thickness_m}
layer_length = ${fparse 4 * film_thickness_m}
c_bulk = ${fparse methane_henry_mol_m3_Pa * pressure_Pa}

[Mesh]
  [layers]
    type = CartesianMeshGenerator
    dim = 2
    dx = '${layer_length}'
    ix = '${mesh_elements_x}'
    dy = '${dielectric_thickness_m} ${film_thickness_m}'
    iy = '${mesh_elements_dielectric} ${mesh_elements_film}'
    subdomain_id = '0 1'
  []
  [shift]
    type = TransformGenerator
    input = layers
    transform = TRANSLATE
    vector_value = '0 ${fparse -dielectric_thickness_m} 0'
  []
  [interface]
    type = SideSetsBetweenSubdomainsGenerator
    input = shift
    primary_block = 1
    paired_block = 0
    new_boundary = interface
  []
  [interface_dielectric]
    type = SideSetsBetweenSubdomainsGenerator
    input = interface
    primary_block = 0
    paired_block = 1
    new_boundary = interface_dielectric
  []
  [pin]
    type = ExtraNodesetGenerator
    input = interface_dielectric
    new_boundary = pressure_pin
    coord = '0 ${film_thickness_m} 0'
  []
[]

[Variables]
  [velocity]
    family = LAGRANGE_VEC
    block = 1
  []
  [p]
    block = 1
  []
  [c]
    block = 1
    initial_condition = ${c_bulk}
  []
  [phi]
    block = 0
  []
[]

[AuxVariables]
  [vel_x]
    block = 1
  []
  [shear_rate]
    family = MONOMIAL
    order = CONSTANT
    block = 1
  []
  [dissipation]
    family = MONOMIAL
    order = CONSTANT
    block = 1
  []
  [field_y]
    family = MONOMIAL
    order = CONSTANT
    block = 0
  []
[]

[AuxKernels]
  [vel_x]
    type = VectorVariableComponentAux
    variable = vel_x
    vector_variable = velocity
    component = x
  []
  [shear_rate]
    type = VariableGradientComponent
    variable = shear_rate
    gradient_variable = vel_x
    component = y
  []
  [dissipation]
    type = ParsedAux
    variable = dissipation
    coupled_variables = 'shear_rate'
    expression = '${liquid_viscosity_Pa_s} * shear_rate * shear_rate'
  []
  [field_y]
    # dphi/dy; the interfacial field magnitude |E| = sigma / eps
    type = VariableGradientComponent
    variable = field_y
    gradient_variable = phi
    component = y
  []
[]

[Kernels]
  [mass]
    type = INSADMass
    variable = p
  []
  [mass_pspg]
    type = INSADMassPSPG
    variable = p
  []
  [momentum_advection]
    type = INSADMomentumAdvection
    variable = velocity
  []
  [momentum_viscous]
    type = INSADMomentumViscous
    variable = velocity
  []
  [momentum_pressure]
    type = INSADMomentumPressure
    variable = velocity
    pressure = p
    integrate_p_by_parts = true
  []
  [momentum_supg]
    type = INSADMomentumSUPG
    variable = velocity
    velocity = velocity
  []
  [methane]
    type = ADMatDiffusion
    variable = c
    diffusivity = D_ch4
  []
  [electrostatics]
    type = ADMatDiffusion
    variable = phi
    diffusivity = eps
  []
[]

[Materials]
  [liquid]
    type = ADGenericConstantMaterial
    block = 1
    prop_names = 'rho mu D_ch4'
    prop_values = '${liquid_density_kg_m3} ${liquid_viscosity_Pa_s} ${methane_diffusivity_m2_s}'
  []
  [ins]
    type = INSADTauMaterial
    block = 1
    velocity = velocity
    pressure = p
  []
  [dielectric]
    type = ADGenericConstantMaterial
    block = 0
    prop_names = 'eps'
    prop_values = '${permittivity_F_m}'
  []
  [surface_consumption]
    type = ADParsedMaterial
    boundary = interface
    property_name = ch4_rate
    coupled_variables = 'c'
    constant_names = 'k'
    constant_expressions = '${interface_rate_constant_m_s}'
    expression = 'k * c'
  []
  [surface_sink]
    type = ADParsedMaterial
    boundary = interface
    property_name = ch4_flux_bc
    material_property_names = 'ch4_rate'
    expression = '-ch4_rate'
  []
[]

[BCs]
  [Periodic]
    [x]
      auto_direction = x
      variable = 'velocity p c phi'
    []
  []
  [no_slip]
    type = VectorFunctionDirichletBC
    variable = velocity
    boundary = interface
    function_x = 0
    function_y = 0
  []
  [sheared_top]
    type = VectorFunctionDirichletBC
    variable = velocity
    boundary = top
    function_x = ${wall_speed}
    function_y = 0
  []
  [pressure_pin]
    type = DirichletBC
    variable = p
    boundary = pressure_pin
    value = 0
  []
  [methane_saturated_top]
    type = DirichletBC
    variable = c
    boundary = top
    value = ${c_bulk}
  []
  [methane_sink]
    type = ADMatNeumannBC
    variable = c
    boundary = interface
    value = 1
    boundary_material = ch4_flux_bc
  []
  [contact_charge]
    # eps dphi/dn = sigma on the dielectric side of the interface
    type = NeumannBC
    variable = phi
    boundary = interface_dielectric
    value = ${interface_charge_C_m2}
  []
  [grounded_back_contact]
    type = DirichletBC
    variable = phi
    boundary = bottom
    value = 0
  []
[]

[Postprocessors]
  [interface_ch4_mol_m3]
    type = SideAverageValue
    variable = c
    boundary = interface
  []
  [ch4_consumption_total]
    type = ADSideIntegralMaterialProperty
    property = ch4_rate
    boundary = interface
  []
  [ch4_consumption_mol_m2_s]
    type = ParsedPostprocessor
    pp_names = 'ch4_consumption_total'
    expression = 'ch4_consumption_total / ${layer_length}'
  []
  [ch4_supply_total]
    type = SideDiffusiveFluxIntegral
    variable = c
    boundary = top
    diffusivity = ${methane_diffusivity_m2_s}
  []
  [ch4_supply_mol_m2_s]
    type = ParsedPostprocessor
    pp_names = 'ch4_supply_total'
    expression = 'ch4_supply_total / ${layer_length}'
  []
  [interface_field_V_m]
    type = SideAverageValue
    variable = field_y
    boundary = interface_dielectric
  []
  [interface_charge_C_m2]
    type = ConstantPostprocessor
    value = ${interface_charge_C_m2}
  []
  [gauss_charge_C_m2]
    type = ParsedPostprocessor
    pp_names = 'interface_field_V_m'
    expression = '${permittivity_F_m} * interface_field_V_m'
  []
  [wall_shear_rate_s_inv]
    type = SideAverageValue
    variable = shear_rate
    boundary = interface
  []
  [top_shear_rate_s_inv]
    type = SideAverageValue
    variable = shear_rate
    boundary = top
  []
  [wall_work_W_m2]
    type = ParsedPostprocessor
    pp_names = 'top_shear_rate_s_inv'
    expression = '${liquid_viscosity_Pa_s} * top_shear_rate_s_inv * ${wall_speed}'
  []
  [dissipation_total]
    type = ElementIntegralVariablePostprocessor
    variable = dissipation
    block = 1
  []
  [viscous_dissipation_W_m2]
    type = ParsedPostprocessor
    pp_names = 'dissipation_total'
    expression = 'dissipation_total / ${layer_length}'
  []
  [mean_velocity_m_s]
    type = ElementAverageValue
    variable = vel_x
    block = 1
  []
  [elements]
    type = NumElements
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu mumps'
  line_search = none
  automatic_scaling = true
  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-16
  nl_max_its = 50
[]

[Outputs]
  csv = true
  exodus = true
  print_linear_residuals = false
[]
