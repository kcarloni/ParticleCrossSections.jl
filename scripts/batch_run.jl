
modpath = "/Users/kiara/julia_packages/ParticleCrossSections.jl"

begin
rng_log_E = 0:0.2:7
proj_A, proj_Z = 12, 6
targ_A, targ_Z = 1,1

modelname = "DpmjetIII193"
# modelname = "QGSJetII04"
# modelname = "UrQMD34"
# modelname = "EposLHC"
    
    outfile = modpath * "/data/sig_inel_nucnuc/" * modelname * "/"
    outfile = outfile * "proj$(proj_A).$(proj_Z)-targ$(targ_A).$(targ_Z).txt"

    rm( outfile, force=true )

    @time for log_E in 0:0.2:7
        cmd = `julia $modpath/scripts/run_chromo.jl $log_E $proj_A $proj_Z $targ_A $targ_Z $outfile $modelname`
        try
            run(cmd)
        catch
            println("failed $log_E")
        end
    end
end