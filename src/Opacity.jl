module Opacity

export opacity

include("./utils/Constants.jl")
# include("utils/HelperFunctions.jl")
using Roots, QuadGK

@enum Regime strong_field=1 medium_field=2 weak_field_or_low_temperature=3 zero_field=4
@enum Channel neutron=1 proton=2


# this code is dogshit
# i'll make it better Later tm
electron_energy(k_ze, n_e, B) = sqrt(ELECTRON_MASS^2 + k_ze^2 + 2 * n_e * ELEM_CHARGE * B)

"""
    electron_density(electron_chemical_potential, B, T)

Calculate the electron number density as a function of chemical potential.
"""
function electron_density(electron_chemical_potential, bare_magnetic_field, temperature)
    mue_gtr_me = electron_chemical_potential > ELECTRON_MASS
    n_max = round(Int, floor(((mue_gtr_me*electron_chemical_potential + TEMP_STEPS*temperature)^2 - ELECTRON_MASS^2) /
                (2 * ELEM_CHARGE * bare_magnetic_field)))
    result = 0
    for n_e in 0:n_max+1
        k_z_max = sqrt(max(0, (mue_gtr_me*electron_chemical_potential + TEMP_STEPS*temperature)^2 - ELECTRON_MASS^2 - 2 * n_e * ELEM_CHARGE * bare_magnetic_field))
        level_contrib = ELEM_CHARGE * bare_magnetic_field / (4*pi^2) * quadgk(k_ze -> fermi_dirac(electron_energy(k_ze, n_e, bare_magnetic_field), electron_chemical_potential, temperature), -k_z_max, k_z_max)[1]
        if n_e > 0
            level_contrib *= 2
        end
        result += level_contrib
    end
    return result
end


"""
    fermi_dirac(E, mu, T)

Calculate the Fermi-Dirac distribution for a given ``E, mu, T``.
"""
function fermi_dirac(energy, chemical_potential, temperature)
    if (energy-chemical_potential)/temperature > 50 return 0 end;
    if (energy-chemical_potential)/temperature < -50 return 1 end;
    return 1 / (1 + exp((energy - chemical_potential) / temperature))
end


"""
    reduced_matrix_element(neutron_spin, proton_spin, neutrino_angle, proton_level, electron_level)

Calculate the reduced matrix element for the opacity.
"""
function reduced_matrix_element(neutron_spin, proton_spin, neutrino_angle; is_strong_field=false)
    cos_neutrino_angle = cos(neutrino_angle)

    # Both particles spin up
    if neutron_spin == proton_spin == 0.5
        plus_term = 2 * (G_V + G_A)^2 * (1 + cos_neutrino_angle) 
        minus_term = 2 * (G_V - G_A)^2 * (1 - cos_neutrino_angle) * (1-is_strong_field)
        return plus_term + minus_term
    end

    # Neutron spin up, proton spin down
    if neutron_spin == 0.5 && proton_spin == -0.5
        return 8 * G_A^2 * (1 + cos_neutrino_angle) * (1-is_strong_field)
    end

    # Neutron spin down, proton spin up
    if neutron_spin == -0.5 && proton_spin == 0.5
        return 8 * G_A^2 * (1 - cos_neutrino_angle)
    end

    # Both spin down
    if neutron_spin == proton_spin == -0.5
        plus_term = 2 * (G_V + G_A)^2 * (1 - cos_neutrino_angle) * (1-is_strong_field)
        minus_term = 2 * (G_V - G_A)^2 * (1 + cos_neutrino_angle)
        return plus_term + minus_term
    end

    throw(ArgumentError("neutron_spin and proton_spin should be +/- 0.5"))
end

# function sign_sum(channel, regime, neutrino_energy, )

#     energy_shift = NUCLEON_MASS_SPLITTING +
#                     NEUTRON_AMM * neutron_spin * NUCLEAR_MAGNETON * magnetic_field / ELEM_CHARGE -
#                     PROTON_AMM * proton_spin * NUCLEAR_MAGNETON * magnetic_field / ELEM_CHARGE

#     E_0_plus = neutrino_energy + energy_shift
#     E_0_minus = neutrino_energy - energy_shift
#     E_0_plus_prime = E_0_plus + 2 * proton_spin * mu_times_B
#     E_0_minus_prime = E_0_minus - 2 * proton_spin * mu_times_B
# end

energy_shift(neutron_spin, proton_spin, magnetic_field, has_gminus2) = NUCLEON_MASS_SPLITTING -
    NEUTRON_GYROMAGNETIC_RATIO * neutron_spin * NUCLEAR_MAGNETON * magnetic_field / ELEM_CHARGE +
    (PROTON_GYROMAGNETIC_RATIO - 2 * has_gminus2) * proton_spin * NUCLEAR_MAGNETON * magnetic_field / ELEM_CHARGE


function opacity(channel, baryon_density, proton_fraction, neutrino_energy, neutrino_angle, magnetic_field, temperature)
    if (temperature / 10)^(3/2) * (197.3^3 * 0.16) / baryon_density < 30
        throw(ArgumentError("Inputs are outside the Maxwell-Boltzmann range."))
    end

    # Quantities we define for convenience
    relevant_nucleon_mass = channel == "p" ? PROTON_MASS : NEUTRON_MASS
    M_times_T = relevant_nucleon_mass * temperature

    neutrino_momentum_z = neutrino_energy * cos(neutrino_angle)
    neutrino_momentum_perp = neutrino_energy * sin(neutrino_angle)

    proton_density = proton_fraction * baryon_density
    neutron_density = baryon_density - proton_density
    
    z = - (neutrino_momentum_perp^2 * magnetic_field) / (2 * M_times_T * (magnetic_field + M_times_T))
    G_tilde = (G_F^2 * COS2_CABIBBO_ANGLE * magnetic_field) /
                (8*pi * cosh((NEUTRON_GYROMAGNETIC_RATIO * magnetic_field) / (4 * M_times_T)) * cosh((PROTON_GYROMAGNETIC_RATIO * magnetic_field) / (4 * M_times_T)))
    G_tilde_prime = (G_F^2 * COS2_CABIBBO_ANGLE * magnetic_field) /
                        (8*pi * cosh(NEUTRON_GYROMAGNETIC_RATIO * magnetic_field / (4 * M_times_T)) * cosh((PROTON_GYROMAGNETIC_RATIO-2) * magnetic_field / (4 * M_times_T)))

    regime = which_regime(magnetic_field, temperature, neutrino_energy)

    # Proton channel
    if channel == 'p'
        if regime == strong_field
            prefactor = G_tilde_prime * proton_density * tanh(magnetic_field / (2 * M_times_T))
            
            # sum opacities over nucleon spin pairs
            sum_over_spins = 0
            for neutron_spin in (-0.5, 0.5), proton_spin in (-0.5, 0.5)
                recurring_exponential = 1 - exp(-magnetic_field/M_times_T)
                
                # M_red and the exponential
                reduced_matrix_element_term = reduced_matrix_element(neutron_spin, proton_spin, neutrino_angle; is_strong_field=true)
                exponential_term = exp(-(PROTON_GYROMAGNETIC_RATIO-2) * proton_spin * magnetic_field / (2 * M_times_T))
                
                # long term inside the parentheses
                bracketed_term_1 = 1 / recurring_exponential
                bracketed_term_2 = neutron_density * exp(-NEUTRON_GYROMAGNETIC_RATIO * neutron_spin * magnetic_field / (2 * M_times_T)) * (pi / M_times_T)^(3/2)
                bracketed_term_3 = M_times_T / (magnetic_field + M_times_T * recurring_exponential)
                bracketed_term_4 = exp(-neutrino_momentum_perp^2 / 2 * 
                                        (recurring_exponential / (magnetic_field + M_times_T * recurring_exponential)) - 
                                        neutrino_momentum_z^2 / (4 * M_times_T))

                sum_over_spins += reduced_matrix_element_term * exponential_term * (bracketed_term_1 - bracketed_term_2 * bracketed_term_3 * bracketed_term_4)
            end
            return prefactor * sum_over_spins / RECIP_CM_TO_MEV
        else
            println("lol sorry nothing but strong field yet")
            return -1
        end

    # neutron channel
    elseif channel == 'n'
        if regime == strong_field
            prefactor = G_tilde_prime * neutron_density / 2
            electron_chemical_potential = find_zero(mu -> proton_density - electron_density(mu, magnetic_field/ELEM_CHARGE, temperature), (0, 200))
            
            # sum opacities over nucleon spin pairs
            sum_over_spins = 0
            for neutron_spin in (-0.5, 0.5), proton_spin in (-0.5, 0.5)
                E_0_plus_prime = neutrino_energy + energy_shift(neutron_spin, proton_spin, magnetic_field, true)
                recurring_exponential = 1 - exp(-magnetic_field/M_times_T)
                
                # first few factors inside the sum
                reduced_matrix_element_term = reduced_matrix_element(neutron_spin, proton_spin, neutrino_angle; is_strong_field=true)
                fermi_dirac_term = 1 - fermi_dirac(E_0_plus_prime, electron_chemical_potential, temperature)
                exponential_term = exp((NEUTRON_GYROMAGNETIC_RATIO * neutron_spin * magnetic_field) / (2 * M_times_T)) # dropping a 2 in the denominator since spin = \pm 1/2
                
                # long term inside the parentheses
                bracketed_term_1 = proton_density * exp((PROTON_GYROMAGNETIC_RATIO * proton_spin - 1) * magnetic_field / (2 * M_times_T)) /
                                        cosh((PROTON_GYROMAGNETIC_RATIO * magnetic_field) / (4 * M_times_T))
                bracketed_term_2 = (pi / M_times_T)^(3/2) * sinh(magnetic_field / (2 * M_times_T))
                bracketed_term_3 = 2 * M_times_T / (magnetic_field + M_times_T * recurring_exponential)
                bracketed_term_4 = exp(-neutrino_momentum_perp^2 / (2 * magnetic_field + 2 * M_times_T) * 
                                        (magnetic_field * (1 - recurring_exponential) / (magnetic_field + M_times_T * recurring_exponential)) - 
                                        neutrino_momentum_z^2 / (4 * M_times_T))
                sum_over_spins += reduced_matrix_element_term * fermi_dirac_term * exponential_term *
                                (1 - bracketed_term_1 * bracketed_term_2 * bracketed_term_3 * bracketed_term_4)
                
                return prefactor * sum_over_spins / RECIP_CM_TO_MEV

            end
        end
    end
end

gg(a, b) = a / b > 10
gtrsim(a, b) = a / b > 0.2

function which_regime(field_strength, temperature, neutrino_energy)
    if gtrsim(field_strength, PROTON_MASS*temperature)
        return strong_field
    elseif gg(PROTON_MASS*temperature, field_strength) && gtrsim(field_strength, temperature^2)
        return medium_field
    elseif gg(max(neutrino_energy^2, temperature^2), field_strength) && gtrsim(field_strength, NEUTRON_MASS*NUCLEON_MASS_SPLITTING/NEUTRON_GYROMAGNETIC_RATIO) ||
            gg(neutrino_energy^2, field_strength) && gtrsim(field_strength, PROTON_MASS*temperature)
        return weak_field_or_low_temp
    elseif gg(NEUTRON_MASS*NUCLEON_MASS_SPLITTING/NEUTRON_GYROMAGNETIC_RATIO, field_strength)
        return zero_field
    end
end

# which_regime(ELEM_CHARGE*GAUSS_TO_MEV2*5e16, 3, 10)

println(opacity('n', (197.3^3 * 0.16 * 0.001), 0.25, 10, pi/2, ELEM_CHARGE*GAUSS_TO_MEV2*1e17, 3))
println(opacity('p', (197.3^3 * 0.16 * 0.001), 0.25, 10, pi/2, ELEM_CHARGE*GAUSS_TO_MEV2*1e17, 1))

end # module Opacity