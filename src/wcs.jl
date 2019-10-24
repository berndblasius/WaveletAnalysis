##  Wavelet cross spectrum in julia
##  Author Bernd Blasius (bernd.blasius@gmail.com)

using DSP, Statistics

include("wavelet.jl")
include("tools.jl")


function wcs(x1,x2,dt)
# calculate the wavelet cross spectrum (amplitude and phase) and
# wavelet coherence of two time series
    # x1, x2 : the two vectors of the time series
    # dt: sampling time interval
    # this function assumes that both time series cover the same time range
    # and are sampled equally
    # no checks are made to ensure this!!

    pad = 1       # pad the time series with zeros
    nvoices = 32   # total number of voices
    noctave = 4   # total number of octaves
    dj = 1 ./ nvoices
    s0 = 2*dt     # start at a scale of 2*dt
    j1 = noctave * nvoices
    mother = "MORLET"
    k0 = 6
    
    # do the wavelet transform  
    # period and scale should be returned the same
    wave1, period, scale, coi = wavelet(x1,dt,pad,dj,s0,j1,mother,k0)
    wave2, period, scale, coi = wavelet(x1,dt,pad,dj,s0,j1,mother,k0)
    
    power1 = abs2.(wave1)   # wavelet power spetrum
    power2 = abs2.(wave2)   # wavelet power spetrum

    # calculate wavelet cross spectrum (CWS)
    wcs = wave1 .* conj.(wave2)

    # smoothing
    spower1 = smooth(power1,dt,nvoices,period)
    spower2 = smooth(power2,dt,nvoices,period)
    swcs    = smooth(wcs,sw,dt,nvoices,period)  

    amplitude = abs.(swcs)   # amplitude of wavelet cross spectrum
    phase     = angle.(swcs) # phase of wavelet cross spectrum

    # wavelet coherence
    wco = swcs_abs ./ (sqrt.(spower1) .* sqrt.(spower2))

    return spower1, spower2, amplitude, phase, wco
end