
module ParticleCrossSections

const modpath = pkgdir(@__MODULE__) 

using AstroParticleUnits
using ParticleTable: Particle

using DelimitedFiles
using Interpolations

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

include("nuele.jl")
export σ_tot_nuele
export calc_neutrino_prod_E_th

include("nuDIS.jl")
export σ_tot_CCDIS, σ_tot_NCDIS

end # module ParticleCrossSections
