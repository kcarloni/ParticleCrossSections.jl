
"""
    σ_inel_pp( Tp )

Total inelastic proton-proton collision cross-section,
as a function of the projectile proton energy `Tp`.
"""
σ_inel_pp( Tp ) = σ_inel_pp_Kafexhiu2014( Tp )

"""
    σ_inel_pnuc( Tp, At )

Total inelastic proton-nucleus collision cross-section,
as a function of the projectile proton energy `Tp`.
"""
function σ_inel_pnuc( Tp, At ) 

    if At == 1;
        return σ_inel_pp( Tp )
    end
    return σ_inel_pnuc_Kafexhiu2014( Tp, At )
end

"""
    σ_inel_pnuc( Tp, Ap, At, )

Total inelastic nucleus-nucleus collision cross-section,
as a function of the projectile energy per nucleon `Tp`.

Note: because it is a function of the energy per nucleon, 
this function is symmetric in the choice of target / projectile.
"""
function σ_inel_nucnuc( Tp, Ap, At ) 

    # projectile energy per nucleon is the same 
    # in the rest frame of either nucleon
    if Ap == 1
        return σ_inel_pnuc( Tp, At )
    elseif At == 1
        return σ_inel_pnuc( Tp, Ap )
    end

    return σ_inel_nucnuc_Kafexhiu2014( Tp, Ap, At )
end


"""
    w_AB( Tp, Ap, At )

The average number of wounded nucleons in a nucleus-nucleus collision, as a function of the projectile energy per nucleon `Tp`.

w_AB is given by: w_AB = (A σ_pB + B σ_pA)/σ_AB
"""
function w_AB( Tp, Ap, At )

    σ_t = σ_inel_pnuc( Tp, At )

    # target energy is the same in the rest frame 
    # of the projectile because they have the same mass
    σ_p = σ_inel_pnuc( Tp, Ap )
    σ_pt = σ_inel_nucnuc( Tp, Ap, At )

    return (Ap * σ_t + At * σ_p) / σ_pt 
end

"""
The secondary meson multiplicity enhancement factor,
as a function of the projectile energy per nucleon `Tp`,
according to the wounded nucleon model:
Białas+ (1976) [https://doi.org/10.1016/0550-3213(76)90329-1]

The wounded nucleon model assumes that the average meson production multiplicity in a nucleus-nucleus interaction is proportional to the number of nucleons that underwent at least one inelastic collision. The total number of mesons (e.g. π0) produced is then given by 
    ⟨n_AB (π0)⟩ = 1/2 w_AB ⟨n_pp(π0)⟩
"""
f_secondary_mult( Tp, Ap, At ) = 1/2 * w_AB( Tp, Ap, At )

# ==================================

# parametrizations from Sihver+ (1993)

    """
        σ_inel_nucnuc( Ap, At )

    Total inelastic nucleus-nucleus collision cross-section,
    which is energy-independent if the projectile energy per nucleon
    Tp ≥ 100 MeV/nucleon.

    Parametrization from Sihver+ (1993) [https://doi.org/10.1103/PhysRevC.47.1225]
    """
    function σ_inel_nucnuc_Sihver1993( Tp, Ap, At )
        
        if Tp <= 100MeV
            println( "Sihver1993 parametrization invalid for Tp < 100MeV/nucleon!")
            return NaN * mb
        end

        Ap13 = Ap^(1/3)
        At13 = At^(1/3)

        b0 = 1.581 - 0.876 * ( 1/Ap13 + 1/At13 )

        return 58.1mb * ( Ap13 + At13 - b0 * ( 1/Ap13 + 1/At13 ) )^2

    end

    """
        σ_inel_pnuc( At )

    Total inelastic proton-nucleus collision cross-section,
    which is energy-independent if the projectile (proton) energy Tp ≥ 200 MeV. 

    Parametrization from Sihver+ (1993) [https://doi.org/10.1103/PhysRevC.47.1225]
    """
    function σ_inel_pnuc_Sihver1993( Tp, At )

        if Tp <= 200MeV
            println( "Sihver1993 parametrization invalid for Tp < 200MeV!")
            return NaN * mb
        end

        At13 = At^(1/3)
        b0 = 2.247 - 0.915 * (1 + 1/At13 )

        return 58.1mb * ( 1 + At13 - b0 * ( 1 + 1/At13 ) )^2
    end

#

# parametrizations from Kafexhiu (2014)

    """
        σ_inel_pp( Tp )

    Total inelastic proton-proton collision cross-section,
    as a function of the projectile energy `Tp`.

    Parametrization from Kafexhiu (2014) [arXiv:1406.7369]
    """
    function σ_inel_pp_Kafexhiu2014( Tp )

        # Tp_th = 0.2797GeV
        Tp_th = 2m_pi0 + m_pi0^2/2m_proton

        r = Tp / Tp_th
        return mb * 
            ( 30.7 - 0.96 * log(r) + 0.18 * log(r)^2 ) *  
            ( 1 - r^(-1.9) )^3
    end

    """
        σ_inel_nucnuc( Tp, Ap, At )

    Total inelastic nucleus-nucleus collision cross-section,
    which is based on the parametrization from Sihver1993 `σ_inel_nucnuc_Sihver1993`
    but is enhanced logarithmically in projectile energy per nucleon `Tp`.
    """
    function σ_inel_nucnuc_Kafexhiu2014( Tp, Ap, At )

        Tp0 = 1TeV 
        return σ_inel_nucnuc_Sihver1993( Tp, Ap, At ) * ( 
            1 + log(  max( 1, σ_inel_pp_Kafexhiu2014(Tp)/σ_inel_pp_Kafexhiu2014(Tp0) ) )
        )
    end

    """
        σ_inel_nucnuc( Tp, At )

    Total inelastic proton-nucleus collision cross-section,
    which is based on the parametrization from Sihver1993 `σ_inel_pnuc_Sihver1993`
    but is enhanced logarithmically in projectile (proton) energy `Tp`.
    """
    function σ_inel_pnuc_Kafexhiu2014( Tp, At )

        Tp0 = 1TeV 
        return σ_inel_pnuc_Sihver1993( Tp, At ) * ( 
            1 + log(  max( 1, σ_inel_pp_Kafexhiu2014(Tp)/σ_inel_pp_Kafexhiu2014(Tp0) ) )
        )
    end

#


# # --------------------------------------

"""
    σ_inel_nucnuc_Tripathi1996( Tp, proj::Isotope, targ::Isotope )

Total inelastic nucleus-nucleus collision cross-section,
as a function of the projectile energy per nucleon `Tp`.

Parametrization originally from 
- Tripathi (1996) 
- Tripathi (1999) [https://doi.org/10.1016/0168-583X(96)00331-X],
"""
function σ_inel_nucnuc_Tripathi1996( Tp, proj::Isotope, targ::Isotope )

    # projectile is always the lighter isotope...
    if targ.A < proj.A
        A_targ, Z_targ = targ.A, targ.Z
        targ = isotopes[proj.Z, proj.A]
        proj = isotopes[Z_targ, A_targ]
    end

    Eproj = Tp * proj.A 

    # projectile energy per nucleon
    Tp_in_MeV = Tp / 1MeV

    # # center of mass energy for internal proton-proton collision
    # Ecm = sqrt( 2m_proton * Tp + 2m_proton^2 )

    # center of mass energy for *full system*
    m_proj = natural( proj.mass.val.val * u"u"; base=MeV )
    m_targ = natural( targ.mass.val.val * u"u"; base=MeV )
    Ecm = sqrt( 2 * m_targ * Eproj + m_targ^2 + m_proj^2 )
    # Ecm = E * Ap * At / (Ap + At) (??)
    Ecm_in_MeV = Ecm / 1MeV

    At13 = targ.A^(1/3)
    Ap13 = proj.A^(1/3)

    # radius of colliding ion
    r0 = 1.1 * u"fm" 

    # multiplication factor ( Xm in Tripathi99 )
    f = 1 

    E = Tp_in_MeV

    # --------------------
    # B = 

        # # rp, rt = hard sphere nuclear radius
        # rp = 1.29 * u"fm" * Ap13
        # rt = 1.29 * u"fm" * At13

        # rp, rt calculation from Geant4
        rp = 1.29 * 0.6 * 1.36u"fm" * Ap13 / r0
        rt = 1.29 * 0.6 * 1.36u"fm" * At13 / r0 

        # Coulomb barrier distance
        R = rp + rt + 1.2 * (Ap13 + At13) / Ecm_in_MeV^(1/3)

        # energy-dependent coulomb interaction barrier
        B = 1.44 * proj.Z * targ.Z / R 
    #

    # --------------------
    # δE = 

        # mass asymmetry term
        S = At13 * Ap13 / ( At13 + Ap13 )

        # D, T1 are reaction pair specific
        T1 = 40

        D = 1.75
        # proton-nucleus case
        if proj == isotopes[1,1];       D = 2.05
        # alpha-nucleus case
        elseif proj == isotopes[2,4];   D = 2.77 - (8e-3 * targ.A) + (1.8e-5 * targ.A^2) - 0.8 / (1 + exp((250-E)/75))
        end

        # transparency + pauli blocking
        CE = D * (1 + exp(-E/T1)) - 0.292 * exp( -E/792 ) * cos( 0.229 * E^0.453 )

        # correction term
        # - overlapping volume of colliding system
        # - energy dependence of nucleus transparency
        # - isospin effects
        δE = 1.85*S + 0.16*S*(Ecm_in_MeV)^(-1/3) - CE + 
            0.91 * (targ.A - 2*targ.Z)*proj.Z/(proj.A * targ.A)
    #   

    # --------------------
    # Rc = 

        # coulomb multiplier for light systems
        Rc = 1

        # from Tripathi99 Table 2
        if proj == isotopes[1,1]

            if     targ == isotopes[1,2];   Rc = 13.5
            elseif targ == isotopes[2,3];   Rc = 21 
            elseif targ == isotopes[2,4];   Rc = 27

            # Lithium
            elseif targ.Z == 3;             Rc = 2.2
            end

        elseif proj == isotopes[1,2]

            if     targ == isotopes[1,2];   Rc = 13.5
            elseif targ == isotopes[2,4];   Rc = 13.5
            
            # Carbon
            elseif targ.Z == 6 ;            Rc = 6.0
            end

        elseif proj == isotopes[2,4]

            # Tantalum
            if     targ.Z == 73;   Rc = 0.6

            # Gold
            elseif targ.Z == 79;   Rc = 0.6
            end

        end
    #

    # --------------------

    # modified Bradt-Peters form
    # B = energy-dependent coulomb barrier 
    # f = multiplication factor 
    σR = pi*r0^2 * (Ap13 + At13 + δE )^2 * ( 1- Rc * B / Ecm_in_MeV ) * f

end

