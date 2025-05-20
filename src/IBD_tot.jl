

"""
from Giunti + Kim, eq. (12.17)
"""
function σ_tot_IBD( E )

    if E < 1.806MeV; return 0 * cm^2; end

    E_ele = E - 1.293MeV
    p_ele = sqrt( E_ele^2 - m_electron^2 )
    return 9.56e-44 * cm^2 * (E_ele * p_ele/1MeV^2)
end
