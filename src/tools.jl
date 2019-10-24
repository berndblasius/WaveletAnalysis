## Some auxilliary functions, used for continuous wavelet analysis
## Author Bernd Blasius (bernd.blasius@gmail.com)



function scramble(x,n)
# sample n random values from time series x
# to create surrogate data
    len = length(x)
    x[rand(1:len,n)]
end

function smooth(x,dt,nvoice,scale,sw=0.6)
    # this function needs julia package DSP.jl
    xf = similar(x)  # smoothed data
    nscales, ntimes = size(x)
 
    # smooth in time, by filtering with a Gaussian
    for i=1:nscales
      s = scale[i]
      b = [exp(-0.5*(j*dt/s)^2)/(s*2*pi)  for j=1:ntimes]
      xf[i,:] = filtfilt(b,x[i,:])  
    end
 
    # smooth in scale according to Si et al (2008)
    # window sw=0.6 corresponds to decorrelation length of Morlet wavelet
    sw_tot = trunc(Int,sw*nvoice)
    b = (1/sw_tot)*ones(sw_tot)     # flat (rectangular) window -> moving average
    for i=1:ntimes
      xf[:,i] = filtfilt(b,xf[:,i])
    end
 
    xf
 end