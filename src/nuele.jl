
function σ_tot_nuele( E, nu::Particle, T_th=0MeV )

    if E > 10GeV
        return σ_tot_nuele_Gandhi1995( E, nu )
    else
        return σ_tot_nuele_GiuntiKim( E, nu, T_th )
    end

end


# ==================================

# parametrizations from Giunti + Kim ()

function calc_neutrino_prod_E_th( m_sec, m_targ )
    ( m_sec^2 / 2m_targ ) - m_targ/2
end

"""
from Giunti + Kim Eq. 5.40.
"""
function σ_tot_nuele_to_numu_GiuntiKim( E ) 
    s = 2m_electron * E
    m_Wboson = Particle("W+").mass.value * c^2

    E_th = calc_neutrino_prod_E_th( m_muon, m_electron)
    if E < E_th
        return 0 * cm^2
    end

    return unnatural( cm^2,
            G_fermi^2 * s/4π * 4 *  ( 1 - m_muon^2/s )
    )
end

"""
from Giunti + Kim Eq. 5.32.
"""
function σ_tot_nuele_to_nuele_GiuntiKim( E, T_th, g1, g2 )

    σ0 = unnatural( cm^2, 2 * G_fermi^2 * m_electron^2 / pi )
    # σ0 = 88.06e-46*cm^2

    T_max = 2E^2 / (m_electron + 2E)

    return uconvert( cm^2, 
        σ0 / m_electron * (
            (g1^2 + g2^2) * (T_max - T_th) - 
            (g2^2 + g1 * g2 * m_electron/2E) * ( T_max^2 - T_th^2 ) / E + 
            1/3 * g2^2 * ( T_max^3 - T_th^3 ) / E^2 
        )
    )
end

"""
    σ_tot_nuele_GiuntiKim( E, nu::Particle, T_th=0MeV )

Total neutrino-electron cross-section, for neutrino energy `E` and minimum electron kinetic energy `T_th`.
"""
function σ_tot_nuele_GiuntiKim( E, nu::Particle, T_th=0MeV )

    g = [1.,1.]
    if nu == Particle(12); 
        g[1] = 1/2 + sin_sq_θW
        g[2] = sin_sq_θW
        return σ_tot_nuele_to_nuele_GiuntiKim( E, T_th, g... )

    elseif nu == Particle(-12);
        g[2] = 1/2 + sin_sq_θW
        g[1] = sin_sq_θW
        return σ_tot_nuele_to_nuele_GiuntiKim( E, T_th, g... )

    elseif nu == Particle(14)
        g[1] = -1/2 + sin_sq_θW
        g[2] = sin_sq_θW 

        return (
            σ_tot_nuele_to_nuele_GiuntiKim( E, T_th, g... ) +
            # add inverse muon decay
            σ_tot_nuele_to_numu_GiuntiKim( E )
        )

    elseif nu == Particle(16)
        return 0 * cm^2

    elseif nu in Particle.([-14, -16])
        g[2] = -1/2 + sin_sq_θW
        g[1] = sin_sq_θW 
        return σ_tot_nuele_to_nuele_GiuntiKim( E, T_th, g... )

    end
end

# ==================================

# parametrizations from Gandhi et. al (1995)

function σ_tot_nuele_Gandhi1995( E, nu::Particle )
    
    if ( nu == Particle(12) )
        return σ_tot_nuele_Gandhi1995( E, nu, Particle(11) )

    elseif  ( nu == Particle(-12) )
        return σ_tot_nuele_Gandhi1995( E, nu, Particle(11) ) +
            # all other W-boson decay channels
            ( 1 - 0.108166 ) * σ_tot_nuele_Gandhi1995( E, nu, Particle(13) )

    elseif  ( nu == Particle(14) )
        return (
            σ_tot_nuele_Gandhi1995( E, nu, Particle(11) ) + 
            # inverse muon decay
            σ_tot_nuele_Gandhi1995( E, nu, Particle(13) )
        )


    elseif  ( nu == Particle(16) )
        return (
            σ_tot_nuele_Gandhi1995( E, nu, Particle(11) ) + 
            # inverse muon decay
            σ_tot_nuele_Gandhi1995( E, nu, Particle(15) )
        )

    elseif  ( nu in ( Particle(-14), Particle(-16) ) )
        return σ_tot_nuele_Gandhi1995( E, nu, Particle(11) )
    end

end



function σ_tot_nuele_Gandhi1995( E, nu::Particle, l::Particle )

    # integrals over inelasticity
    # f1(A) = ∫₀¹dy (1+Ay)⁻² × (1-y)²
    # f2(A) = ∫₀¹dy (1+Ay)⁻²
    # f3(A,B) = ∫₀¹dy (1+Ay)⁻¹ × (1+B(1-y))⁻¹
    # f4(A) = ∫₀¹dy (1+Ay)⁻¹ × (1-y)²

    f1(A) = 1/A + 2/A^2 - (2A+2)/A^3 * log1p(A)
    f2(A) = 1/(A+1)
    f3(A,B) = ( log1p(B) + log1p(A) )/( A + (1+A)*B )
    f4(A) = (1/A + 2/A^2 + 1/A^3) * log1p(A) - 3/2A - 1/A^2

    s = 2 * m_electron * E
    σ0 = unnatural( cm^2, G_fermi^2 * s / 4pi )

    m_Zboson = Particle("Z0").mass.value * c^2
    m_Wboson = Particle("W+").mass.value * c^2
    ΓW = Particle("W+").width.value * c^2

    Re = 2 * sin_sq_θW
    Le = 2 * sin_sq_θW - 1

    A = s/m_Zboson^2
    B = s/m_Wboson^2 
    C = (ΓW / m_Wboson)^2 

    if ( nu in ( Particle(14), Particle(16) ) ) && ( l == Particle(11) )
        return σ0 * ( Re^2 * f1(A) + Le^2 * f2(A) )

    elseif ( nu in ( Particle(-14), Particle(-16) ) ) && ( l == Particle(11) )
            return σ0 * ( Re^2 * f2(A) + Le^2 * f1(A) )


    # inverse muon decay
    elseif ( nu == Particle(14) ) && ( l == Particle(13) )
        E_th = calc_neutrino_prod_E_th( m_muon, m_electron)
        if E < E_th
            return 0 * cm^2
        end
        return σ0 * 4 * ( 
            1 - (m_muon^2 - m_electron^2)/s )^2 * f2(B)

    # inverse tau decay (?)
    elseif ( nu == Particle(16) ) && ( l == Particle(15) )
        m_tau = Particle("tau").mass.value * c^2
        E_th = calc_neutrino_prod_E_th( m_tau, m_electron)
        if E < E_th
            return 0 * cm^2
        end
        return σ0 * 4 * ( 
            1 - (m_tau^2 - m_electron^2)/s )^2 * f2(B)


    elseif ( nu == Particle(12) ) && ( l == Particle(11) )
        return σ0 * ( Re^2 * f1(A) + Le^2 * f2(A) + 
            4 * f2(B) + 4Le * f3( A, B )
        )

    elseif ( nu == Particle(-12) ) && ( l == Particle(11) )
        return σ0 * (
            Re^2 * f2(A) + Le^2 * f1(A) + 
            4/3 * ( (1 - B)^2 + C )^(-1) + 
            4Le * ( 1 - A ) * ( (1 - B)^2 + C )^(-1) * f4( A )
        )

    elseif ( nu == Particle(-12) ) && ( l == Particle(13) )

        E_th = calc_neutrino_prod_E_th( m_muon, m_electron)
        if E < E_th
            return 0 * cm^2
        end
    
        return σ0 * 4/3 * 
            ( 1 - (m_muon^2 - m_electron^2)/s )^2 * 
            ( (1 - B)^2 + C )^(-1)

    end

end

# function dσ_numuele_numuele( E, y )

#     # Eq. 15

#     Re = 2 * sin_sq_θW
#     Le = 2 * sin_sq_θW - 1 

#     s = 2 * m_electron * E
#     return G_fermi^2 * s / 4pi * 
#         (1 + s * y / M_Zboson^2 )^(-2) * 
#         ( Re^2 * (1-y)^2 - Le^2 )

# end

# function dσ_numubarele_numubarele( E, y )

#     # Eq. 15

#     Re = 2 * sin_sq_θW
#     Le = 2 * sin_sq_θW - 1 

#     s = 2 * m_electron * E
#     return G_fermi^2 * s / 4pi * 
#         (1 + s * y / M_Zboson^2 )^(-2) * 
#         ( Re^2 + Le^2 * (1-y)^2 )

# end
