
# differential cross-sections for X production in nucleus-nucleus collisions. 

"""
    dσ_nucnuc( Eprim, Esec, prim, targ, sec )

Differential cross-section dσ/dEsec(Eprim, Esec) for production of `sec` with energy `Esec` in fixed-target collisions between `prim` with energy `Eprim` and `targ`.
"""
function dσ_nucnuc( Eprim, Esec, prim, targ, sec )

    dσ_AA = 0. * mb/GeV 
    try
    # first option: use AAfrag table
        dσ_AA = dσ_nucnuc_AAfrag( Eprim, Esec, prim, targ, sec )

    catch
    # # second recourse: wounded nucleon approximation 
    #     dσ_AA = dσ_nucnuc_woundednuc( Eprim, Esec, prim, targ, sec )
    end

    return dσ_AA
end    

# =========================================
# based on AAfrag tables: 

    function load_dxs_aafrag202()

        He = isotopes[2,4]
        C = isotopes[6,12]
        Al = isotopes[13, 26]
        Fe = isotopes[26, 56]

        # NamedTuple{ (:prim, :targ) }.(
        prim_targ_pairs = [
            ( p, p ),
            ( p, He ),
            ( He, p ),
            ( He, He ),
            ( C, p ),
            ( Al, p ),
            ( Fe, p ),
        ] 
        #)
        tup_strs = ["p_p_04", "p_He04", "He_p04", "HeHe04", "C_p_04", "Al_p04", "Fe_p04"]
        # tup_strs_lowE = ["p_p_04L", "p_He04L", "He_p04L", "HeHe04L"]

        secondaries = [
            [p, Particle(-2212)],
            [Particle("n"), Particle(-2112)],
            [Particle(12), Particle(-12), Particle(14), Particle(-14)],
            [Particle(11), Particle(-11)],
            [Particle(22)],
            # + antideuteron, antihelium-3
            # Isotope(-1,-2), Isotope(-2,-3)
        ]
        sec_strs = ["pap", "nan", "nu", "el", "gam"]
        all_secs = reduce(vcat, secondaries)

        # # for (ad = antideuteron, ah = antihelium-3 )
        # ["p_p_04", "p_He04", "He_p04", "HeHe04", "", "", "", "ap_p04", "apHe04"]

        fluxes = Matrix{Function}(undef, length(tup_strs), length(all_secs) )

        for (j, t_str) in enumerate( tup_strs )

            # index into all_secs 
            itot = 0

            for (i, s_str) in enumerate( sec_strs )
                for (sub_i, s) in enumerate( secondaries[i] )
                    
                    itot += 1
                    fname = s_str * "_" * t_str 
                    tupmat = load_aafrag_dxs_from_file( fname, sub_i )

                    if j <= 4
                        tupmat_low = load_aafrag_dxs_from_file( fname * "L", sub_i )
                        fluxes[j, itot] = make_aafrag_dxs_itp(
                            tupmat, tupmat_low 
                        )
                    else
                        fluxes[j, itot] = make_aafrag_dxs_itp(
                            tupmat
                        )
                    end

                end
            end
        end
        return (pt=prim_targ_pairs, s=all_secs), fluxes 
    end

    function load_aafrag_dxs_from_file( fname, ind )

        fpath = modpath * "/data/dsig_pp/aafrag202/"
        out = readdlm( fpath * fname  )

        # indices where each sub-table begins
        ind_start = findall( out[:,1] .== 0.0 )

        ind_start_end = collect( zip( ind_start, [ (ind_start[2:end] .- 1)..., size(out, 1) ] ) )

        i0, i1 = ind_start_end[ind]

        E_prim = out[(i0+1):i1, 1] * GeV
        T_sec = out[i0, 2:end] * GeV
        Esq_dσdE = out[(i0+1):i1, 2:end ] * mb * GeV

        rng_prim = try_make_range( log10.( E_prim/1GeV), 1e-1 )
        rng_sec  = try_make_range( log10.( T_sec/1GeV), 1e-1 )

        mat = Esq_dσdE/(1mb*GeV)
        mat[ mat .< 1e-50 ] .= exp10.(-100)
        mat .= log10.( mat )

        return (rng_prim, rng_sec), mat

        # itp = linear_interpolation( (rng_prim, rng_sec), Esq_dσdE,
        #     # extrapolation_bc=zero(eltype(Esq_dσdE)) 
        #     extrapolation_bc=Line()
        # )
        # return (E1, E2) -> itp( log10(E1/1GeV), log10(E2/1GeV) )
    end

    function make_aafrag_dxs_itp( tupmat )
        tup, mat = tupmat 
        itp = linear_interpolation(
            tup, mat, extrapolation_bc=Line() 
        )

        # eprim0, eprim1 = extrema(tup[1])
        # m0 = zero( eltype( mat ) )
        m0 = 0. * mb * GeV 

        eprim0 = minimum( exp10.(tup[1]) * GeV )
        esec0, esec1 = extrema( exp10.(tup[2]) * GeV )

        function dxs( eprim, esec )
            if !( esec0 < esec < esec1); return m0; end 
            if eprim < eprim0; return m0; end 
            return exp10(
                itp( log10(eprim/1GeV), log10(esec/1GeV) )
            ) * mb * GeV 
        end
        return dxs
    end

    function make_aafrag_dxs_itp( tupmat, tupmatL )
        tup, mat = tupmat 
        tupL, matL = tupmatL 

        itp = make_aafrag_dxs_itp( tupmat )
        itpL = make_aafrag_dxs_itp( tupmatL )

        # eprim0, eprim1 = extrema(tup[1])
        # m0 = zero( eltype( mat ) )
        m0 = 0. * mb * GeV 

        eprim0, eprim1 = extrema( exp10.(tupL[1]) * GeV )
        eprim2 = minimum( exp10.(tup[1]) * GeV )

        function dxs( eprim, esec )
            if eprim < eprim0; return m0; 
            elseif eprim < eprim1; return itpL( eprim, esec )
            elseif eprim > eprim2; return itp( eprim, esec )
            else
                dxsL = itpL( eprim, esec )
                dxs = itp( eprim, esec )

                return dxsL + (dxs - dxsL)/( eprim2 - eprim1) * (eprim - eprim1)
            end
        end
        return dxs
    end

    const aafrag_combos, dxsTables_aafrag = load_dxs_aafrag202()

    """
        dσ_nucnuc_AAfrag( Eprim, Esec, prim, targ, sec )

    Differential cross-section dσ/dEsec(Eprim, Esec) for production of `sec` with energy `Esec` in fixed-target collisions between `prim` with energy `Eprim` and `targ`.

    The cross-sections are based on interpolations of tabulated QGSJet-II-04m runs, provided by AAfrag (Ostapchenko+ 2019, 2021, 2023) arXiV:2110.00496
    """
    function dσ_nucnuc_AAfrag( Eprim, Esec, prim, targ, sec )

        j = findfirst( sec .== aafrag_combos.s )
        i = findfirst( Ref( (prim, targ) ) .== aafrag_combos.pt )
        
        # try swapping order (?) -- requires boosting secondary energy (?)
        i_swap = findfirst(  Ref( (targ, prim) ) .== aafrag_combos.pt )

        if isnothing(i) && isnothing(i_swap);    
            error("differential cross-section for (primary, target) combination ($(prim.name), $(targ.name)) combination is not tabulated by AAfrag!"); 

        elseif isnothing(i);
            error("not implemented -- differential cross-section for (primary, target) combination ($(prim.name), $(targ.name)) combination is tabulated for swapped configuration by AAfrag!"); 

        elseif isnothing(j); 
            error("differential cross-section for secondary $(sec.name) is not tabulated by AAfrag!"); 
        else
            # ? for whatever reason results don't match AAfrag unless I assume the tables are actually in mb, not mb GeV 
            return dxsTables_aafrag[i,j]( Eprim, Esec ) / Esec /1GeV 


            # return dxsTables_aafrag[i,j]( Eprim, Esec ) / Esec^2
        end
    end

    # """
    # dσ_pp_AAfrag( Eprim, Esec, sec )

    # Convenience function specific to p-p collisions.
    # """
    # function dσ_pp_AAfrag( Eprim, Esec, sec )

    #     j = findfirst( sec .== aafrag_combos.s )
    #     i = 1

    #     if isnothing(j); 
    #         error("differential cross-section for secondary $sec is not tabulated by AAfrag!"); 
    #     else
    #         return dxsTables_aafrag[i,j]( Eprim, Esec ) / Esec^2
    #     end
    # end




# =========================================
# using wounded nucleon approximation

# ?! off by a factor of 10 compared to AAfrag tabulated results ... ?!
function dσ_nucnuc_woundednuc( Eprim, Esec, prim, targ, sec )
    # dσ_pp = dσ_pp_AAfrag( Eprim/prim.A, Esec, sec )
    dσ_pp = dσ_nucnuc_AAfrag( Eprim/prim.A, Esec, p, p, sec )

    dσ_AA = dσ_pp * f_secondary_mult( Eprim/prim.A, prim.A, targ.A )
    return dσ_AA
end