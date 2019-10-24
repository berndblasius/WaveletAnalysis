##  Continuous wavelet analysis in julia
##  Author Bernd Blasius (bernd.blasius@gmail.com)

# for a more involved example, including plotting, wavelet power spectrum
# and significance levels, see test/cwt_test_elnino.jl 


using Statistics
include("wavelet.jl")

# continuoous wavelet transform (cwt)
function cwt(x,dt)
    # x : vector of time series
    # dt: sampling time interval
    pad = 1       # pad the time series with zeros
    nvoices = 32   # total number of voices
    noctave = 4   # total number of octaves
    dj = 1 ./ nvoices
    s0 = 2*dt     # start at a scale of 2*dt
    #n1 = length(y)
    #j1=trunc(Int,(log(n1*dt/s0)/log(2.0))/dj)
    j1 = noctave * nvoices
    mother = "MORLET"
    k0 = 6
    
    x = x .- mean(x)    # normalize
    x = x ./ std(x)
    variance = var(x)

    wave, period, scale, coi = wavelet(x,dt,pad,dj,s0,j1,mother,k0)
    power = abs2.(wave)   # wavelet power spetrum

    # Global wavelet spectrum 
    global_ws = variance*(mean(sqrt.(power),dims=2))   # time-average over all times

    return wave, period, power, global_ws, x, variance, scale, coi, nvoices
end

