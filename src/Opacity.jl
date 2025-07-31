module Opacity

include("utils/Constants.jl")

"""
    reduced_matrix_element(neutron_spin, proton_spin, neutrino_angle, proton_level, electron_level)

Calculate the reduced matrix element for the opacity.
"""
function reduced_matrix_element(neutron_spin, proton_spin, neutrino_angle; is_strong_field=false)
    cos_neutrino_angle = cos(neutrino_angle)
    # is_p0 = proton_level == 0
    # is_e0 = electron_level == 0

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
        minus_term = 2 * (G_V - G_A)^2 * (1 + cos_neutrino_angle) * is_strong_field
        return plus_term + minus_term
    end

    @error "reduced_matrix_element: neutron_spin and proton_spin should be +/- 0.5"
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

function opacity(channel, baryon_density, proton_fraction, neutrino_energy, neutrino_angle, magnetic_field, temperature)
    # Quantities we define for convenience
    mu_times_B = NUCLEAR_MAGNETON * magnetic_field / ELEM_CHARGE
    M_times_T = AVG_NUCLEON_MASS * temperature

    neutrino_momentum_z = neutrino_energy * cos(neutrino_angle)
    neutrino_momentum_perp = neutrino_energy * sin(neutrino_angle)

    proton_density = proton_fraction * baryon_density
    neutron_density = baryon_density - proton_density

    # energy_shift = NUCLEON_MASS_SPLITTING +
    #                 NEUTRON_GYROMAGNETIC_RATIO * neutron_spin * NUCLEAR_MAGNETON * magnetic_field / ELEM_CHARGE -
    #                 PROTON_GYROMAGNETIC_RATIO * proton_spin * NUCLEAR_MAGNETON * magnetic_field / ELEM_CHARGE
    
    # E_0_plus = neutrino_energy + energy_shift
    # E_0_minus = neutrino_energy - energy_shift
    # E_0_plus_prime = E_0_plus + 2 * proton_spin * mu_times_B
    # E_0_minus_prime = E_0_minus - 2 * proton_spin * mu_times_B
    
    z = -(neutrino_momentum_perp^2 * magnetic_field) / (2 * M_times_T * (magnetic_field + M_times_T))
    G_tilde = (G_F^2 * COS2_CABIBBO_ANGLE * magnetic_field) /
                (16*pi * cosh(NEUTRON_GYROMAGNETIC_RATIO * magnetic_field/2) * cosh(PROTON_GYROMAGNETIC_RATIO * magnetic_field/2))
    num = (G_F^2 * COS2_CABIBBO_ANGLE * magnetic_field)
    den1 = 16*pi * cosh(NEUTRON_GYROMAGNETIC_RATIO * magnetic_field/4)
    den2 = cosh((PROTON_GYROMAGNETIC_RATIO-2) * magnetic_field/4)
    G_tilde_prime = (G_F^2 * COS2_CABIBBO_ANGLE * magnetic_field) /
                        (16*pi * cosh(NEUTRON_GYROMAGNETIC_RATIO * magnetic_field / (4*M_times_T)) * cosh((PROTON_GYROMAGNETIC_RATIO-2) * magnetic_field / (4*M_times_T)))

    if channel == 'p'
        if gtrsim(magnetic_field, PROTON_MASS*temperature)
            prefactor = G_tilde_prime * proton_density * tanh(magnetic_field / (2*M_times_T))
            sign_sum = 0
            for proton_spin in (-0.5, +0.5)
                for neutron_spin in (-0.5, +0.5)
                    energy_shift = NUCLEON_MASS_SPLITTING +
                                        NEUTRON_GYROMAGNETIC_RATIO * neutron_spin * NUCLEAR_MAGNETON * magnetic_field / ELEM_CHARGE -
                                        PROTON_GYROMAGNETIC_RATIO * proton_spin * NUCLEAR_MAGNETON * magnetic_field / ELEM_CHARGE
    
                    E_0_plus = neutrino_energy + energy_shift
                    E_0_minus = neutrino_energy - energy_shift
                    E_0_plus_prime = E_0_plus + 2 * proton_spin * mu_times_B
                    E_0_minus_prime = E_0_minus - 2 * proton_spin * mu_times_B
                                    #TODO: which masses go where?
                    recurring_exponential = 1 - exp(-magnetic_field/M_times_T)
                    reduced_matrix_element_single_spin_pair = reduced_matrix_element(neutron_spin, proton_spin, neutrino_angle; is_strong_field=true)
                    first_exponential = exp(-(PROTON_GYROMAGNETIC_RATIO-2) * proton_spin * magnetic_field / (2 * M_times_T))
                    bracketed_term_1 = 1/recurring_exponential
                    bracketed_term_2 = proton_density * exp(-NEUTRON_GYROMAGNETIC_RATIO * neutron_spin * magnetic_field / (2 * M_times_T)) * (pi / M_times_T)^(3/2)
                    bracketed_term_3 = M_times_T / (magnetic_field + M_times_T * recurring_exponential)
                    bracketed_term_4 = exp(-neutrino_momentum_perp^2/2 * 
                                            (recurring_exponential / (magnetic_field + M_times_T * recurring_exponential)) - 
                                            neutrino_momentum_z^2 / (4 * M_times_T))
                    sign_sum += reduced_matrix_element_single_spin_pair * first_exponential * (bracketed_term_1 - bracketed_term_2 * bracketed_term_3 * bracketed_term_4)
                end
            end
            return prefactor * sign_sum / RECIP_CM_TO_MEV
        else
            println("lol sorry nothing but strong field yet")
            return -1
        end
    elseif channel == 'n'
        println("lol no neutron channel yet")
        return -1
    end
end

gg(a, b) = a / b > 10
gtrsim(a, b) = a / b > 0.2

function which_regime(field_strength, temperature, neutrino_energy)
    println(field_strength / (PROTON_MASS*temperature))
    println("-")
    println(PROTON_MASS*temperature/field_strength, " ", field_strength/temperature^2)
    println("-")
    println(neutrino_energy^2/field_strength, " ", field_strength/(PROTON_MASS*temperature))
    if gtrsim(field_strength, PROTON_MASS*temperature)
        println("strong field")
    elseif gg(PROTON_MASS*temperature, field_strength) && gtrsim(field_strength, temperature^2)
        println("medium field")
    elseif gg(max(neutrino_energy^2, temperature^2), field_strength) && gtrsim(field_strength, NEUTRON_MASS*NUCLEON_MASS_SPLITTING/NEUTRON_GYROMAGNETIC_RATIO)
        println("weak field")
    elseif gg(neutrino_energy^2, field_strength) && gtrsim(field_strength, PROTON_MASS*temperature)
        println("low T")
    elseif gg(NEUTRON_MASS*NUCLEON_MASS_SPLITTING/NEUTRON_GYROMAGNETIC_RATIO, field_strength)
        println("zero field")
    end
end

# which_regime(ELEM_CHARGE*GAUSS_TO_MEV2*5e16, 3, 10)

println(opacity('p', (197.3^3 * 0.16 * 0.001), 0.25, 10, pi/4, ELEM_CHARGE*GAUSS_TO_MEV2*1e17, 3))

end # module Opacity