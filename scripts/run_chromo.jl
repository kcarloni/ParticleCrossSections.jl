
# proj_dir = "/Users/kiara/harvard_grad/research_projects/project_TANDEM/"
# using Pkg
# Pkg.activate( proj_dir )

using DelimitedFiles
using PyCall
chromo = pyimport("chromo")

log_E   = parse(Float64, ARGS[1])
proj_A, proj_Z  = parse.(Int, ARGS[2:3])
targ_A, targ_Z  = parse.(Int, ARGS[4:5])

outfile = ARGS[6]
modelname = ARGS[7]

k = chromo.kinematics.FixedTarget(
    proj_A * exp10(log_E) * chromo.constants.GeV, 
    (proj_A, proj_Z), (targ_A, targ_Z)
)

g = getproperty( chromo.models, Symbol(modelname) )( k )

# out = g_DPM.cross_section() 
# xs = out.prod

# :total, :elastic, :inelastic, :prod, 
# :b_elastic, :quasielastic,
# :coherent...

cols = [
    "log[ projectile_energy / GeV ]",
    "xs_production [mb]",
]
vals = [
    log_E,
    g._inel_or_prod_cross_section
]

mkpath( dirname(outfile) )

delim = ','
if !isfile(outfile)
    println("creating file $outfile")
    open( outfile, "w") do io 
        writedlm( io, reshape(cols, 1, :), delim; )
    end
end

open( outfile, "a") do io 
    writedlm( io, reshape(vals, 1, :), delim )
end

