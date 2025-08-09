module Opacity
export opacity
export ProtonChannel, NeutronChannel
export StrongField, MediumField, WeakFieldOrLowTemperature, ZeroField, Unspecified

include("./utils/Constants.jl")
using Roots, QuadGK, BenchmarkTools

abstract type CaptureChannel end
struct ProtonChannel <: CaptureChannel end
struct NeutronChannel <: CaptureChannel end

abstract type Regime end
struct StrongField <: Regime end
struct MediumField <: Regime end
struct WeakFieldOrLowTemperature <: Regime end
struct ZeroField <: Regime end
struct Unspecified <: Regime end

electron_energy(electron_momentum_z, electron_level, magnetic_field) = sqrt(ELECTRON_MASS^2 + electron_momentum_z^2 + 2 * electron_level * magnetic_field)


"""
    electron_density_nonzero_field(electron_chemical_potential, magnetic_field, temperature)

Calculate the electron number density as a function of chemical potential.
"""
function electron_density_nonzero_field(electron_chemical_potential, magnetic_field, temperature)
    mue_gtr_me = electron_chemical_potential > ELECTRON_MASS
    energy_term = electron_chemical_potential*mue_gtr_me + TEMP_STEPS*temperature
    max_energy_level = floor(Int, (energy_term^2 - ELECTRON_MASS^2) / (2 * magnetic_field))
    
    result = 0.0
    for n_e in 0:max_energy_level
        k_z_max = sqrt(max(0, (mue_gtr_me*electron_chemical_potential + TEMP_STEPS*temperature)^2 - ELECTRON_MASS^2 - 2 * n_e * magnetic_field))
        level_contrib = magnetic_field / (4*pi^2) * quadgk(k_ze -> fermi_dirac(electron_energy(k_ze, n_e, magnetic_field),
                                                            electron_chemical_potential, temperature), -k_z_max, k_z_max; rtol=1e-4)[1]
        if n_e > 0
            level_contrib *= 2
        end
        result += level_contrib
    end
    return result
end


"""
    electron_density_zero_field(electron_chemical_potential, temperature)

Calculate the electron number density as a function of chemical potential, in the zero-field case.
"""
function electron_density_zero_field(electron_chemical_potential, temperature)
    upper_bound = max(20 * temperature, 20 * (temperature + electron_chemical_potential))
    return quadgk(k -> k^2 / (pi^2 * (exp((k - electron_chemical_potential) / temperature) + 1)), 0, upper_bound)[1]
end


"""
    find_electron_chemical_potential

Calculate the electron's chemical potential.
"""
find_electron_chemical_potential(::Union{StrongField, MediumField, WeakFieldOrLowTemperature}, proton_density, magnetic_field, temperature) =
    find_zero(mu -> proton_density - electron_density_nonzero_field(mu, magnetic_field, temperature), 5)

find_electron_chemical_potential(::ZeroField, proton_density, magnetic_field, temperature) = 
    find_zero(mu -> proton_density - electron_density_zero_field(mu, temperature), 20)

# function electron_density(electron_chemical_potential, magnetic_field, temperature)
#     mue_gtr_me = electron_chemical_potential > ELECTRON_MASS

#     n_max = round(Int, floor(((mue_gtr_me*electron_chemical_potential + TEMP_STEPS*temperature)^2 - ELECTRON_MASS^2) /
#                 (2 * magnetic_field)))
#     result = 0
#     for n_e in 0:n_max+1
#         k_z_max = sqrt(max(0, (mue_gtr_me*electron_chemical_potential + TEMP_STEPS*temperature)^2 - ELECTRON_MASS^2 - 2 * n_e * magnetic_field))
#         level_contrib = magnetic_field / (4*pi^2) * quadgk(k_ze -> fermi_dirac(electron_energy(k_ze, n_e, magnetic_field), electron_chemical_potential, temperature), -k_z_max, k_z_max; rtol=1e-4)[1]
#         if n_e > 0
#             level_contrib *= 2
#         end
#         result += level_contrib
#     end
#     return result
# end

"""
    fermi_dirac(E, mu, T)

Calculate the Fermi-Dirac distribution for a given ``E, mu, T``.
"""
function fermi_dirac(energy, chemical_potential, temperature)
    if (energy-chemical_potential)/temperature > 15 return 0 end;
    if (energy-chemical_potential)/temperature < -15 return 1 end;
    return 1 / (1 + exp((energy - chemical_potential) / temperature))
end


"""
    reduced_matrix_element(neutron_spin, proton_spin, neutrino_angle, proton_level, electron_level)

Calculate the reduced matrix element for the opacity.
"""
function reduced_matrix_element(neutron_spin, proton_spin, neutrino_angle; n_e_equals_0=false)
    cos_neutrino_angle = cos(neutrino_angle)

    # Both particles spin up
    if neutron_spin == proton_spin == 1
        plus_term = 2 * (G_V + G_A)^2 * (1 + cos_neutrino_angle) 
        minus_term = 2 * (G_V - G_A)^2 * (1 - cos_neutrino_angle) * (1-n_e_equals_0)
        return plus_term + minus_term
    end

    # Neutron spin up, proton spin down
    if neutron_spin == 1 && proton_spin == -1
        return 8 * G_A^2 * (1 + cos_neutrino_angle) * (1-n_e_equals_0)
    end

    # Neutron spin down, proton spin up
    if neutron_spin == -1 && proton_spin == 1
        return 8 * G_A^2 * (1 - cos_neutrino_angle)
    end

    # Both spin down
    if neutron_spin == proton_spin == -1
        plus_term = 2 * (G_V + G_A)^2 * (1 - cos_neutrino_angle) * (1-n_e_equals_0)
        minus_term = 2 * (G_V - G_A)^2 * (1 + cos_neutrino_angle)
        return plus_term + minus_term
    end

    throw(ArgumentError("neutron_spin and proton_spin should be +/- 1"))
end


"""
    function get_shifted_energy(channel, neutrino_energy, neutron_spin, proton_spin, magnetic_field; has_gminus2=false)

Calculate the shifted energy of the neutrino (corresponding to \$E_0^{\\pm}\$ in the paper).
"""
function get_shifted_energy(channel, neutrino_energy, neutron_spin, proton_spin, magnetic_field; has_gminus2=false)
    plus_or_minus = channel isa NeutronChannel ? 1 : -1
    energy_shift = plus_or_minus * (NUCLEON_MASS_SPLITTING -
            NEUTRON_GYROMAGNETIC_RATIO * neutron_spin * magnetic_field / (4 * NEUTRON_MASS * ELEM_CHARGE) +
            (PROTON_GYROMAGNETIC_RATIO - 2 * has_gminus2) * proton_spin * magnetic_field / (4 * PROTON_MASS * ELEM_CHARGE))
    neutrino_energy + energy_shift
end


# Blocking Terms
"""
    function blocking_term

Calculate the 'blocking term' (the term in parentheses within the sum over spins) for each capture channel and field/temperature regime.
"""
## Neutron, strong field
function blocking_term(::NeutronChannel, ::StrongField, _, proton_density, _, proton_spin, neutrino_energy, neutrino_angle, magnetic_field, M_times_T)
    neutrino_momentum_z = neutrino_energy * cos(neutrino_angle)
    neutrino_momentum_perp = neutrino_energy * sin(neutrino_angle)

    recurring_exponential = 1 - exp(-magnetic_field / M_times_T)
    term_1 = proton_density * exp((PROTON_GYROMAGNETIC_RATIO * proton_spin - 1) * magnetic_field / (2 * M_times_T)) /
                                        cosh((PROTON_GYROMAGNETIC_RATIO * magnetic_field) / (4 * M_times_T))
    term_2 = (pi / M_times_T)^(3/2) * sinh(magnetic_field / (2 * M_times_T))
    term_3 = 2 * M_times_T / (magnetic_field + M_times_T * recurring_exponential)
    term_4 = exp(-neutrino_momentum_perp^2 / (2 * magnetic_field + 2 * M_times_T) * 
                            (magnetic_field * (1 - recurring_exponential) / (magnetic_field + M_times_T * recurring_exponential)) - 
                            neutrino_momentum_z^2 / (4 * M_times_T))
    return 1 - term_1 * term_2 * term_3 * term_4
end

## Neutron, medium field
function blocking_term(::NeutronChannel, ::MediumField, _, proton_density, _, proton_spin, neutrino_energy, neutrino_angle, magnetic_field, M_times_T)
    term_1 = proton_density * exp((PROTON_GYROMAGNETIC_RATIO * proton_spin * magnetic_field) / (2 * M_times_T)) /
                                        cosh((PROTON_GYROMAGNETIC_RATIO * magnetic_field) / (4 * M_times_T))
    term_2 = (pi / M_times_T)^(3/2) * sinh(magnetic_field / (2 * M_times_T)) * M_times_T / magnetic_field
    term_3 = exp(-neutrino_energy^2 / (4 * M_times_T))
    return 1 - term_1 * term_2 * term_3
end

## Neutron, weak field/low temp
function blocking_term(::NeutronChannel, ::WeakFieldOrLowTemperature, _, proton_density, _, proton_spin, neutrino_energy, neutrino_angle, magnetic_field, M_times_T)
    term_1 = proton_density * exp((PROTON_GYROMAGNETIC_RATIO * proton_spin * magnetic_field) / (2 * M_times_T)) /
                                        (2 * cosh((PROTON_GYROMAGNETIC_RATIO * magnetic_field) / (4 * M_times_T)))
    term_2 = (pi / M_times_T)^(3/2) * exp(-neutrino_energy^2 / (4 * M_times_T))
    return 1 - term_1 * term_2
end


## Neutron, zero field
function blocking_term(::NeutronChannel, ::ZeroField, _, proton_density, _, proton_spin, neutrino_energy, neutrino_angle, magnetic_field, M_times_T)
    term_1 = (proton_density / 2) * (pi / M_times_T)^(3/2)
    term_2 = exp(-neutrino_energy^2 / (4 * M_times_T))
    return 1 - term_1 * term_2
end


## Proton, strong field
function blocking_term(::ProtonChannel, ::StrongField, neutron_density, _, neutron_spin, _, neutrino_energy, neutrino_angle, magnetic_field, M_times_T)
    neutrino_momentum_z = neutrino_energy * cos(neutrino_angle)
    neutrino_momentum_perp = neutrino_energy * sin(neutrino_angle)

    recurring_exponential = 1 - exp(-magnetic_field/M_times_T)

    term_1 = 1 / recurring_exponential
    term_2 = neutron_density * exp(-NEUTRON_GYROMAGNETIC_RATIO * neutron_spin * magnetic_field / (2 * M_times_T)) * (pi / M_times_T)^(3/2)
    term_3 = M_times_T / (magnetic_field + M_times_T * recurring_exponential)
    term_4 = exp(-neutrino_momentum_perp^2 / 2 * 
                    (recurring_exponential / (magnetic_field + M_times_T * recurring_exponential)) - 
                    neutrino_momentum_z^2 / (4 * M_times_T))
    return term_1 - term_2 * term_3 * term_4
end


## Proton, medium field/weak field/low temp
function blocking_term(::ProtonChannel, ::Union{MediumField, WeakFieldOrLowTemperature}, neutron_density, _, neutron_spin, _, neutrino_energy, neutrino_angle, magnetic_field, M_times_T)
    term_1 = neutron_density * exp((NEUTRON_GYROMAGNETIC_RATIO * neutron_spin * magnetic_field) / (2 * M_times_T)) /
                                        (2 * cosh((NEUTRON_GYROMAGNETIC_RATIO * magnetic_field) / (4 * M_times_T)))
    term_2 = (pi / M_times_T)^(3/2) * exp(-neutrino_energy^2 / (4 * M_times_T))
    return 1 - term_1 * term_2
end


## Proton, zero field
function blocking_term(::ProtonChannel, ::ZeroField, neutron_density, _, neutron_spin, _, neutrino_energy, neutrino_angle, magnetic_field, M_times_T)
    term_1 = (neutron_density / 2) * (pi / M_times_T)^(3/2)
    term_2 = exp(-neutrino_energy^2 / (4 * M_times_T))
    return 1 - term_1 * term_2
end



"""
    get_v_tilde_and_n_e_flag

Get the \$\\tilde{v}\$ term and whether all electrons are in the lowest energy level.
"""
get_v_tilde_and_n_e_flag(::StrongField, shifted_neutrino_energy, magnetic_field) = 1, true

function get_v_tilde_and_n_e_flag(::MediumField, shifted_neutrino_energy, magnetic_field)
    v_tilde = max(shifted_neutrino_energy^2 / magnetic_field, 1)
    n_e_equals_0 = v_tilde == 1
    return v_tilde, n_e_equals_0
end
get_v_tilde_and_n_e_flag(::WeakFieldOrLowTemperature, shifted_neutrino_energy, magnetic_field) = shifted_neutrino_energy^2 / magnetic_field, false
get_v_tilde_and_n_e_flag(::ZeroField, shifted_neutrino_energy, magnetic_field) = 1, false

# Exponential terms #
get_exponential_term(::NeutronChannel, neutron_spin, proton_spin, magnetic_field, M_times_T) = exp((NEUTRON_GYROMAGNETIC_RATIO * neutron_spin * magnetic_field) / (2 * M_times_T))
get_exponential_term(::ProtonChannel, neutron_spin, proton_spin, magnetic_field, M_times_T) = exp((PROTON_GYROMAGNETIC_RATIO * proton_spin * magnetic_field) / (2 * M_times_T))

# Prefactors #
get_prefactor(::NeutronChannel, G_tilde, neutron_density, _, _magnetic_field, _M_times_T) =
    G_tilde * neutron_density / 2
get_prefactor(::ProtonChannel, G_tilde, _, proton_density, magnetic_field, M_times_T) =
    G_tilde * proton_density * sinh(magnetic_field / (2 * M_times_T)) * M_times_T / magnetic_field


"""
    function opacity(channel, baryon_density, proton_fraction, neutrino_energy, neutrino_angle, magnetic_field, temperature; regime=Unspecified())

Calculate the opacity of the capture process in question.
"""
function opacity(channel, baryon_density, proton_fraction, neutrino_energy, neutrino_angle, magnetic_field, temperature; regime=Unspecified())
    # Fermi-Dirac regime not yet implemented
    if (temperature / 10)^(3/2) * (197.3^3 * 0.16) / baryon_density < 30
        throw(ArgumentError("Inputs are outside the Maxwell-Boltzmann range."))
    end

    # Quantities we define for convenience
    M_times_T = AVG_NUCLEON_MASS * temperature
    proton_density = proton_fraction * baryon_density
    neutron_density = baryon_density - proton_density

    # Calculate prefactor
    G_tilde = (G_F^2 * COS2_CABIBBO_ANGLE * magnetic_field) /
                (4*pi * cosh((NEUTRON_GYROMAGNETIC_RATIO * magnetic_field) / (4 * M_times_T)))
    prefactor = get_prefactor(channel, G_tilde, neutron_density, proton_density, magnetic_field, M_times_T)

    # Set regime, unless a specific regime was passed as a keyword argument
    regime = set_regime_if_unspecified(regime, magnetic_field, temperature, neutrino_energy)

    # Calculate the electron's chemical potential
    electron_chemical_potential = find_electron_chemical_potential(regime, proton_density, magnetic_field, temperature)
    
    # Sum over nucleon spins
    sum_over_spins = 0.0
    for neutron_spin in (-1, 1), proton_spin in (-1, 1)
        shifted_neutrino_energy = get_shifted_energy(channel, neutrino_energy, neutron_spin, proton_spin, magnetic_field; has_gminus2=(regime isa StrongField))
        v_tilde_term, n_e_flag = get_v_tilde_and_n_e_flag(regime, shifted_neutrino_energy, magnetic_field)
        
        reduced_matrix_element_term = reduced_matrix_element(neutron_spin, proton_spin, neutrino_angle; n_e_equals_0=n_e_flag)
        fermi_dirac_term = 1 - fermi_dirac(shifted_neutrino_energy, electron_chemical_potential, temperature)
        exponential_term = get_exponential_term(channel, neutron_spin, proton_spin, magnetic_field, M_times_T)

        pauli_blocking_term = blocking_term(channel, regime, neutron_density, proton_density, neutron_spin, proton_spin, neutrino_energy, neutrino_angle, magnetic_field, M_times_T)
        
        sum_over_spins += v_tilde_term * reduced_matrix_element_term * fermi_dirac_term * exponential_term * pauli_blocking_term
    end

    return prefactor * sum_over_spins / RECIP_CM_TO_MEV
end

# Regime picker #
gg(a, b) = a / b > 10
gtrsim(a, b) = a / b > 0.1

strong_field_condition(magnetic_field, temperature) = gtrsim(magnetic_field, AVG_NUCLEON_MASS*temperature)

medium_field_condition(magnetic_field, temperature) = gg(AVG_NUCLEON_MASS*temperature, magnetic_field) && gtrsim(magnetic_field, temperature^2)

weak_field_condition(magnetic_field, temperature, neutrino_energy) = gg(max(neutrino_energy^2, temperature^2), magnetic_field) &&
    gtrsim(magnetic_field, NEUTRON_MASS*NUCLEON_MASS_SPLITTING/NEUTRON_GYROMAGNETIC_RATIO)

low_temperature_condition(magnetic_field, temperature, neutrino_energy) = gg(neutrino_energy^2, magnetic_field) &&
    gtrsim(magnetic_field, AVG_NUCLEON_MASS * temperature)

zero_field_condition(magnetic_field) = gg(AVG_NUCLEON_MASS*NUCLEON_MASS_SPLITTING/NEUTRON_GYROMAGNETIC_RATIO, magnetic_field)

function set_regime_if_unspecified(regime::Regime, magnetic_field, temperature, neutrino_energy)
    if !(regime isa Unspecified)
        return regime
    elseif strong_field_condition(magnetic_field, temperature)
        return StrongField()
    elseif medium_field_condition(magnetic_field, temperature)
        return MediumField()
    elseif weak_field_condition(magnetic_field, temperature, neutrino_energy) || low_temperature_condition(magnetic_field, temperature, neutrino_energy)
        return WeakFieldOrLowTemperature()
    elseif zero_field_condition(magnetic_field)
        return ZeroField()
    else
        error("Could not determine parameter-space regime for eB=$magnetic_field MeV^2, T=$temperature MeV, E_nu=$neutrino_energy MeV")
    end
end

# println(opacity(NeutronChannel(), (197.3^3 * 0.16 * 0.001), 0.25, 10, pi/2, ELEM_CHARGE*GAUSS_TO_MEV2*1e17, 3; regime=StrongField()))
# println(opacity(NeutronChannel(), (197.3^3 * 0.16 * 0.001), 0.25, 10, pi/2, ELEM_CHARGE*GAUSS_TO_MEV2*1e17, 3; regime=MediumField()))
# println(opacity(NeutronChannel(), (197.3^3 * 0.16 * 0.001), 0.25, 10, pi/2, ELEM_CHARGE*GAUSS_TO_MEV2*1e17, 3; regime=WeakFieldOrLowTemperature()))
# println(opacity(NeutronChannel(), (197.3^3 * 0.16 * 0.001), 0.25, 10, pi/2, ELEM_CHARGE*GAUSS_TO_MEV2*1e17, 3; regime=ZeroField()))
# println(opacity(ProtonChannel(), (197.3^3 * 0.16 * 0.001), 0.25, 10, pi/2, ELEM_CHARGE*GAUSS_TO_MEV2*1e17, 3; regime=StrongField()))
# println(opacity(ProtonChannel(), (197.3^3 * 0.16 * 0.001), 0.25, 10, pi/2, ELEM_CHARGE*GAUSS_TO_MEV2*1e17, 3; regime=MediumField()))
# println(opacity(ProtonChannel(), (197.3^3 * 0.16 * 0.001), 0.25, 10, pi/2, ELEM_CHARGE*GAUSS_TO_MEV2*1e17, 3; regime=WeakFieldOrLowTemperature()))
# println(opacity(ProtonChannel(), (197.3^3 * 0.16 * 0.001), 0.25, 10, pi/2, ELEM_CHARGE*GAUSS_TO_MEV2*1e17, 3; regime=ZeroField()))

# @btime opacity(NeutronChannel(), (197.3^3 * 0.16 * 0.001), 0.25, 10, pi/2, ELEM_CHARGE*GAUSS_TO_MEV2*1e17, 3; regime=StrongField())

end # module Opacity