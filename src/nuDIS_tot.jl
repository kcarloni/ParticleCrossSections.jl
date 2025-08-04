
function load_xs_CooperSarkar2011( typ, interaction )

    out = readdlm( modpath * "/data/sig_nuDIS/coopersarkar2011_" * typ * ".txt", skipstart=1 )

    ene = out[:,1] * GeV
    # xs_CC = out[:,2] * u"pb"
    # xs_NC = out[:,9] * u"pb"

    if interaction == "CC"
        xs = out[:,2] * u"pb"
    elseif interaction == "NC"
        xs = out[:,9] * u"pb"
    end

    itp = linear_interpolation(
        log10.(ene/1GeV), log10.(xs/1u"pb"), extrapolation_bc=Line()
    )

    return ( E -> u"pb" * exp10( itp( log10(E/GeV) ) ) )
end

const σ_tot_CCDIS_CS11_nu = load_xs_CooperSarkar2011("nu", "CC")
const σ_tot_CCDIS_CS11_nubar = load_xs_CooperSarkar2011("nubar", "CC")
const σ_tot_NCDIS_CS11_nu = load_xs_CooperSarkar2011("nu", "NC")
const σ_tot_NCDIS_CS11_nubar = load_xs_CooperSarkar2011("nubar", "NC")


function σ_tot_CCDIS( E, nu )
    if nu in Particle.( (12, 14, 16) )
        return σ_tot_CCDIS_CS11_nu(E)
    else
        return σ_tot_CCDIS_CS11_nubar(E)
    end
end

function σ_tot_NCDIS( E, nu )
    if nu in Particle.( (12, 14, 16) )
        return σ_tot_NCDIS_CS11_nu(E)
    else
        return σ_tot_NCDIS_CS11_nubar(E)
    end
end

# const modpath = "/Users/kiara/julia_packages/ParticleCrossSections.jl"
# load_xs_CooperSarkar2011( "nu" )


# ===============================

function load_xs_Gandhi1995(
    typ, interaction
)

    logE = 1:12
    if typ == "nu" && interaction == "CC"
        σ = [ 
            0.777e-37
            0.697e-36
            0.625e-35
            0.454e-34
            0.196e-33
            0.611e-33
            0.176e-32
            0.478e-32
            0.123e-31
            0.301e-31
            0.706e-31
            0.159e-31
        ] * cm^2
    elseif typ == "nu" && interaction == "NC"
        σ = [ 
            0.242e-37
            0.217e-36
            0.199e-35
            0.155e-34
            0.745e-34
            0.252e-33
            0.748e-33
            0.207e-32
            0.540e-32
            0.134e-31
            0.316e-31
            0.715e-31
        ] * cm^2
    elseif typ == "nubar" && interaction == "CC"
        σ = [ 
            0.368e-37
            0.349e-36
            0.338e-35
            0.292e-34
            0.162e-33
            0.582e-33
            0.174e-32
            0.477e-32
            0.123e-31
            0.301e-31
            0.706e-31
            0.159e-30
        ] * cm^2
    elseif typ == "nubar" && interaction == "NC"  
        σ = [ 
            0.130e-37
            0.122e-36
            0.120e-35
            0.106e-34
            0.631e-34
            0.241e-33
            0.742e-33
            0.207e-32
            0.540e-32
            0.134e-31
            0.316e-31
            0.715e-31
        ] * cm^2 
    end

    logσ = log10.( σ/1mb )
    itp = linear_interpolation( logE, logσ )

    return (E -> exp10( itp( log10.(E/1GeV) )) * 1mb )
end

const σ_tot_CCDIS_G95_nu    = load_xs_Gandhi1995("nu", "CC")
const σ_tot_CCDIS_G95_nubar = load_xs_Gandhi1995("nubar", "CC")
const σ_tot_NCDIS_G95_nu    = load_xs_Gandhi1995("nu", "NC")
const σ_tot_NCDIS_G95_nubar = load_xs_Gandhi1995("nubar", "NC")

function σ_tot_CCDIS_G95( E, nu; approx=false )

    is_nu = nu in Particle.( (12, 14, 16) )

    if approx
        if is_nu
            return 2.69e-36 * cm^2 * (E / 1GeV)^(0.402)
        else
            return 2.53e-36 * cm^2 * (E / 1GeV)^(0.404)
        end
    else
        if is_nu
            return σ_tot_CCDIS_G95_nu(E)
        else
            return σ_tot_CCDIS_G95_nubar(E)
        end
    end
end

function σ_tot_NCDIS_G95( E, nu; approx=false )

    is_nu = nu in Particle.( (12, 14, 16) )

    if approx
        if is_nu
            return 1.06e-36 * cm^2 * (E / 1GeV)^(0.408)
        else
            return 0.98e-36 * cm^2 * (E / 1GeV)^(0.410)
        end
    else
        if is_nu
            return σ_tot_NCDIS_G95_nu(E)
        else
            return σ_tot_NCDIS_G95_nu(E)
        end
    end
end

# function σ_tot_CCDIS_Gandhi1995( E, nu )
#     logE = 1:12

#     logσ = log10.( σ/1mb )
#     itp = 

# end