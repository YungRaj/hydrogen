# Electrochemical methane-to-carbon cell: 1-D tertiary current distribution.
#
#   x = 0                      x = electrolyte_thickness_m
#   anode | electrolyte gap (phi, c_CH4) | cathode + CH4-saturated bulk
#
# Anode:   CH4 -> C(s) + n (H+ or H2O-derived) + n e-   (n = electrons_per_ch4)
# Cathode: hydrogen evolution, n e- -> n/2 H2            (overall CH4 -> C + 2 H2)
#
# Governing equations (steady, isothermal, excess supporting electrolyte):
#   div(kappa grad phi) = 0          ionic charge conservation (Ohm's law)
#   div(D grad c)       = 0          dissolved CH4 (neutral: no migration)
# Butler-Volmer boundary currents, with the anode reaction first order in c:
#   i_a = i0_a [ (c/c_bulk) exp(aa F eta_a / RT) - exp(-ac F eta_a / RT) ]
#   eta_a = V_cell - phi - E_a          (cathode metal is the 0 V reference)
#   i_c = i0_c [ exp(-ac_c F eta_c / RT) - exp(aa_c F eta_c / RT) ],  eta_c = -phi - E_c
# The anode exchange current is exchanged with Cantera each coupling
# iteration (anode_exchange_current_A_m2); iteration 1 uses the case guess.
#
# Every top-level variable below is overridden from hydrogen_case.json by the
# runner. The defaults only make the file independently checkable.

# --- geometry / operating / properties (hydrogen_case.json) ---
electrode_area_m2 = 1e-2
electrolyte_thickness_m = 1e-3
temperature_K = 353.15
pressure_Pa = 101325
methane_mass_flow_kg_s = 1e-6
applied_potential_V = 1.4
ionic_conductivity_S_m = 10
electronic_conductivity_S_m = 1e5
methane_diffusivity_m2_s = 1e-9

# --- model_parameters (hydrogen_case.json) ---
methane_henry_mol_m3_Pa = 1.4e-5
anode_equilibrium_potential_V = 0.2
cathode_equilibrium_potential_V = 0.0
anode_transfer_coefficient = 0.5
cathode_transfer_coefficient = 0.5
cathode_exchange_current_A_m2 = 10
anode_exchange_current_guess_A_m2 = 1e-2
ch4_current_fraction = 1.0
cathode_h2_faradaic_efficiency = 1.0
electrons_per_ch4 = 4
mesh_elements = 200
coupling_iteration = 1

# --- exchanged with Cantera (runner override from iteration 2) ---
anode_exchange_current_A_m2 = ${anode_exchange_current_guess_A_m2}

# --- derived ---
faraday = 96485.33212
gas_constant = 8.314462618
f_rt = ${fparse faraday / (gas_constant * temperature_K)}
c_bulk = ${fparse methane_henry_mol_m3_Pa * pressure_Pa}

[Mesh]
  [gap]
    type = GeneratedMeshGenerator
    dim = 1
    nx = ${mesh_elements}
    xmax = ${electrolyte_thickness_m}
  []
[]

[Variables]
  [phi]
    initial_condition = ${fparse -cathode_equilibrium_potential_V}
  []
  [c]
    initial_condition = ${c_bulk}
  []
[]

[Kernels]
  [charge]
    type = ADMatDiffusion
    variable = phi
    diffusivity = kappa
  []
  [methane]
    type = ADMatDiffusion
    variable = c
    diffusivity = D_ch4
  []
[]

[Materials]
  [transport]
    type = ADGenericConstantMaterial
    prop_names = 'kappa D_ch4'
    prop_values = '${ionic_conductivity_S_m} ${methane_diffusivity_m2_s}'
  []
  [anode_kinetics]
    type = ADParsedMaterial
    boundary = left
    property_name = i_anode
    coupled_variables = 'phi c'
    constant_names = 'i0 aa ac frt V Ea cb'
    constant_expressions = '${anode_exchange_current_A_m2} ${anode_transfer_coefficient}
                            ${fparse 1 - anode_transfer_coefficient} ${f_rt}
                            ${applied_potential_V} ${anode_equilibrium_potential_V} ${c_bulk}'
    expression = 'i0 * ((c / cb) * exp(aa * frt * (V - phi - Ea)) - exp(-ac * frt * (V - phi - Ea)))'
  []
  [anode_ch4_sink]
    # D grad(c).n on the anode: CH4 consumed by the CH4 share of the anode current.
    type = ADParsedMaterial
    boundary = left
    property_name = ch4_flux_bc
    material_property_names = 'i_anode'
    constant_names = 'share n F'
    constant_expressions = '${ch4_current_fraction} ${electrons_per_ch4} ${faraday}'
    expression = '-share * i_anode / (n * F)'
  []
  [cathode_kinetics]
    type = ADParsedMaterial
    boundary = right
    property_name = i_cathode
    coupled_variables = 'phi'
    constant_names = 'i0 aa ac frt Ec'
    constant_expressions = '${cathode_exchange_current_A_m2} ${fparse 1 - cathode_transfer_coefficient}
                            ${cathode_transfer_coefficient} ${f_rt} ${cathode_equilibrium_potential_V}'
    expression = 'i0 * (exp(-ac * frt * (-phi - Ec)) - exp(aa * frt * (-phi - Ec)))'
  []
  [cathode_flux]
    type = ADParsedMaterial
    boundary = right
    property_name = cathode_charge_bc
    material_property_names = 'i_cathode'
    expression = '-i_cathode'
  []
[]

[BCs]
  # kappa grad(phi).n = +i_a at the anode (current enters the electrolyte)
  [anode_current]
    type = ADMatNeumannBC
    variable = phi
    boundary = left
    value = 1
    boundary_material = i_anode
  []
  # kappa grad(phi).n = -i_c at the cathode (current leaves the electrolyte)
  [cathode_current]
    type = ADMatNeumannBC
    variable = phi
    boundary = right
    value = 1
    boundary_material = cathode_charge_bc
  []
  [anode_methane]
    type = ADMatNeumannBC
    variable = c
    boundary = left
    value = 1
    boundary_material = ch4_flux_bc
  []
  [bulk_methane]
    type = DirichletBC
    variable = c
    boundary = right
    value = ${c_bulk}
  []
[]

[Postprocessors]
  [anode_current_density_A_m2]
    type = ADSideIntegralMaterialProperty
    property = i_anode
    boundary = left
  []
  [cathode_current_density_A_m2]
    type = ADSideIntegralMaterialProperty
    property = i_cathode
    boundary = right
  []
  [ionic_current_density_A_m2]
    # -kappa dphi/dx at the gap midplane, independent of both BCs
    type = SideDiffusiveFluxIntegral
    variable = phi
    boundary = right
    diffusivity = ${ionic_conductivity_S_m}
  []
  [ch4_supply_mol_m2_s]
    # diffusive CH4 flux entering from the bulk boundary
    type = SideDiffusiveFluxIntegral
    variable = c
    boundary = right
    diffusivity = ${methane_diffusivity_m2_s}
  []
  [ch4_consumption_mol_m2_s]
    type = ParsedPostprocessor
    pp_names = 'anode_current_density_A_m2'
    expression = '${ch4_current_fraction} * anode_current_density_A_m2 / (${electrons_per_ch4} * ${faraday})'
  []
  [anode_surface_ch4_mol_m3]
    type = PointValue
    variable = c
    point = '0 0 0'
  []
  [anode_electrolyte_potential_V]
    type = PointValue
    variable = phi
    point = '0 0 0'
  []
  [cathode_electrolyte_potential_V]
    type = PointValue
    variable = phi
    point = '${electrolyte_thickness_m} 0 0'
  []
  [anode_overpotential_V]
    type = ParsedPostprocessor
    pp_names = 'anode_electrolyte_potential_V'
    expression = '${applied_potential_V} - anode_electrolyte_potential_V - ${anode_equilibrium_potential_V}'
  []
  [cathode_overpotential_V]
    type = ParsedPostprocessor
    pp_names = 'cathode_electrolyte_potential_V'
    expression = '-cathode_electrolyte_potential_V - ${cathode_equilibrium_potential_V}'
  []
  [ohmic_drop_V]
    type = ParsedPostprocessor
    pp_names = 'anode_electrolyte_potential_V cathode_electrolyte_potential_V'
    expression = 'anode_electrolyte_potential_V - cathode_electrolyte_potential_V'
  []
  [elements]
    type = NumElements
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type -snes_force_iteration'
  petsc_options_value = 'lu mumps 1'
  line_search = bt
  automatic_scaling = true
  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-13
  nl_max_its = 100
[]

[Outputs]
  csv = true
  exodus = true
  print_linear_residuals = false
[]
