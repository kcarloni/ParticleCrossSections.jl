

begin
# UHE specific
    fig = Figure()
    ax = Axis( fig[1,1], yscale=log10, xscale=log10 )

    u_x = GeV
    u_y = cm^2

    x = exp10.(2:0.05:12) * GeV

    for (nu, l) in [
        ( 14, 11 ),
        ( -14, 11 ),
        ( 14, 13 ),
        ( 12, 11 ),
        (-12, 11),
        (-12, 13)
    ]

        nu, l = Particle(nu), Particle(l)
        y = σ_tot_nuele_Gandhi1995.( x, nu, l )

        lines!( ax, x/u_x, y/u_y, 
            label=nu.glyph*" → "*l.glyph,
            # label = (@lstring nu.latex * L"\to" * l.glyph),
            linestyle=(:dash,:dense)
        )
    end

    axislegend(ax)
    fig
end

begin
# full energy range
    u_x = 1GeV
    u_y = 1cm^2

    fig = Figure()
    ax = Axis( fig[1,1], yscale=log10, xscale=log10,
        xlabel="neutrino energy $(unit(u_x))",
        ylabel="cross-section $(unit(u_y))"
    )

    x = exp10.(-3:0.05:12) * GeV

    for (i, nu) in enumerate( Particle.( [12, -12, 14, -14, 16, -16] ) )

        y = σ_tot_nuele.( x, nu )
        c = distinct_sequential[12][i+2] #Cycled(i)

        if nu.pdgid.value < 0;
            ls=(:dash, :dense)
        else;
            ls=:solid
        end

        lines!( ax, x/u_x, y/u_y, 
            label="$(nu.glyph)",
            color=c,
            linestyle=ls
        )
    end
    vlines( 11GeV/u_x )

    axislegend(ax, position=:lt)
    save( "figures/xs_nuele.pdf", fig )
    fig
end

