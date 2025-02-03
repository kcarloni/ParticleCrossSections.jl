
module ParticleCrossSections

const modpath = pkgdir(@__MODULE__) 

using AstroParticleUnits
# using DelimitedFiles

using IsotopeTable 
function Base.getproperty(obj::Isotope, sym::Symbol)
    if sym === :Z; return obj.atomic_number
    elseif sym === :A; return obj.mass_number
    else # fallback to getfield
        return getfield(obj, sym)
    end
end



export isotopes

include("pp_inel.jl")
export σ_inel_pp, σ_inel_pnuc, σ_inel_nucnuc
export w_AB, f_secondary_mult


end # module ParticleCrossSections
