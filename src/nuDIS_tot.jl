
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