##  Test continuous wavelet transform and wavelet power spectrum
##  with the example of the elnino time series
##  again closely following Torrence & Compo
## 
##  Author Bernd Blasius (bernd.blasius@gmail.com)

module W

using CSV, DSP, PyPlot, Statistics

include("../src/wavelet.jl")
include("../src/tools.jl")


function test_elnino()
  # read in the data
  df = CSV.read("../data/nino3data.csv"; comment="#") 
  t   = df[:,1]  # time (year)
  sst = df[:,2]  # el nino signal (sea surface temperature SST)
  x   = df[:,3]  # el nino temperature anomaly

  figure("elnino time series")
  clf()
  subplot(2,1,1)
  plot(t,sst)
  xlabel("time")
  ylabel("SST")
  subplot(2,1,2)
  plot(t,x)
  xlabel("time")
  ylabel("anomaly")


  dt = mean(diff(t))   # actually we would need perfectly regular sampling intervals
  pad = 1       # pad the time series with zeros
  nvoices = 32  # total number of voices
  noctave = 5   # total number of octaves
  dj = 1 ./ nvoices
  s0 = 1.0      # starting scale at 1 years
  j1 = noctave * nvoices
  mother = "MORLET"
  k0 = 6
  
  x = x .- mean(x)    # normalize
  x = x ./ std(x)
  variance = var(x)

  wave, period, scale, coi = wavelet(x,dt,pad,dj,s0,j1,mother,k0)   # cwt

  power = abs2.(wave)                         # wavelet power
  global_ws = variance*(mean(power,dims=2))   # time-average --> global wavelet spectrum
  spower = smooth(power,dt,nvoices,period)    # smoothing
  
  # Significance levels: (variance=1 for the normalized SST)
  lag1 = 0.72    # lag-1 autocorrelation for red noise background
  signif,fft_theor = wave_signif(x,dt,scale,0,lag1)
  sig95 = signif .* ones(length(x))'  # expand signif --> (J+1)x(N) array
  sig95 = power ./ sig95              # where ratio > 1, power is significant

  # significance level for global spectrum
  dof = length(x) .- scale  # the -scale corrects for padding at edges
  global_signif,fft_theor = wave_signif(variance,dt,scale,1,lag1,-1,dof)

  # plotting scale range
  lminp = log2(minimum(period))
  lmaxp = log2(maximum(period))
  lrange = (trunc(Int, lminp) :  trunc(Int, lmaxp))

  figure("wavelet power spectrum")
  clf()
  pcolormesh(t, log2.(period),spower, cmap=PyPlot.cm.jet, alpha = 0.5)
  plot(t,log2.(coi),"k")           # cone-of-influence, anything "below" is dubious
  contour(t, log2.(period),sig95, [-99,1], colors="black")   # 95% significance value
  ylabel("log2(period length)")
  xlabel("time")
  ylim(lmaxp,lminp)
  yticks(lrange, 2 .^ lrange)

  figure("global spectrum")
  clf()
  plot(log2.(period),log2.(global_ws))
  plot(log2.(period),log2.(global_signif),"k.")  # 95% significance line
  xlabel("log2(period)")
  ylabel("log2(power)")
  xlim(lmaxp,lminp)
  xticks(lrange, 2 .^ lrange)
end

test_elnino()

end # module
