
fpath = "/Users/kiara/julia_packages/ParticleCrossSections.jl/data/"
fname = "crxsecs_fragmentation_Evoli2019_cumulative_modified.dat"


out = readdlm( fpath * fname, '\t', skipstart=1)

parent = [
    isotopes[Z,A] for (Z,A) in eachcol(
        parse.(Int, reduce( hcat, split.( 
            reduce( hcat, split.(out[:,1]) )[1,:], '.' ) ))
    )
]
daughter = [
    isotopes[Z,A] for (Z,A) in eachcol(
        parse.(Int, reduce( hcat, split.( 
            reduce( hcat, split.(out[:,1]) )[2,:], '.' ) ))
    )
]

energy = out[:, 2] * MeV
xsec = out[:, 3] * mb 
xsec_err = out[:,4] * mb

Table( 
    parent=parent,
    daughter=daughter,
    energy=energy,
    xsec=xsec,
    xsec_err=xsec_err
)