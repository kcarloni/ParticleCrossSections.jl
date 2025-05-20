
module ParticleCrossSections

const modpath = pkgdir(@__MODULE__) 

using AstroParticleUnits
using ParticleTable

using DelimitedFiles
using Interpolations

# using PyCall
# const aafragpy = PyNULL()
# function __init__()
#     copy!( aafragpy, pyimport("aafragpy") )
# end

function try_make_range( x, rtol=1e-5 )

    xr = range( x[1], x[end], length(x) )
    xmax = maximum( x )

    if all( abs.( x .- xr ) .< xmax * rtol ); return xr 
    else; error("steps are not approximately equal")
    end
end

using IsotopeTable 
function Base.getproperty(obj::Isotope, sym::Symbol)
    if sym === :Z; return obj.atomic_number
    elseif sym === :A; return obj.mass_number
    else # fallback to getfield
        return getfield(obj, sym)
    end
end
Base.broadcastable( obj::Isotope) = Ref(obj)
export isotopes

function Base.getproperty(obj::Particle, sym::Symbol)
    if sym === :Z; return ParticleTable.Corpuscles.Z(obj)
    elseif sym === :A; return ParticleTable.Corpuscles.A(obj)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end

# local scope constants for particles -- shorter code
const p = Particle("p")


include("pp_inel_tot.jl")
export σ_inel_pp, σ_inel_pnuc, σ_inel_nucnuc
export w_AB, f_secondary_mult

include("pp_inel_diff.jl")
export dσ_nucnuc

include("nuele_tot.jl")
export σ_tot_nuele
export calc_neutrino_prod_E_th

include("IBD_tot.jl")
export σ_tot_IBD

include("nuDIS_tot.jl")
export σ_tot_CCDIS, σ_tot_NCDIS

end # module ParticleCrossSections
