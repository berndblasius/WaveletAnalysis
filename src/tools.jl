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


# some functions for circular Statistics

function unwrap(v,vmin,vmax)
  # unwrap phase angle (phase in multiples of pi)
    if v < vmin
        v += 2.0
    elseif v > vmax
        v -= 2.0
    end
  v
end

function circstats(ph)
    # ph is a vector if phases in units of 2 pi
    n = length(ph)

    c = sum(cos.(ph)) / n
    s = sum(sin.(ph)) / n

    R = sqrt(c*c + s*s)       # order parameter
    var = 1.0 - R             # cicurlar variance
    theta = atan(s,c)         # mean direction
    std = sqrt(-2. * log(R))  # circular standard deviation
    theta, var, std
end


function get_maxInd(wcs_abs, wco, phs, n, ind1, ind2)
    # find maximal cross-spectrum at each time instance
    # within a range of indices (ind1, ind2) - corresponding to large power

    wcs_max_ind = zeros(Int,n)
    wco_max  = zeros(n)
    angle_max = zeros(n)
    for i=1:n
        ind = argmax(wcs_abs[ind1:ind2,i]) + ind1 - 1
        wcs_max_ind[i] = ind
        wco_max[i]   = wco[ind,i]
        angle_max[i] = phs[ind,i] / pi   # we always measure phase in multiples of pi
    end
    angle_max = unwrap.(angle_max, -0.5, 1.5)   # restrict phases to range -0.5 .. 1.5

    wcs_max_ind, angle_max, wco_max
end


function phasehist(angles, indices)
    # calculate histogram of phase angles (given in multiples of pi)
    delta = 0.01
    angle_hist = fit(Histogram, angles[indices],closed=:right, -0.5:delta:1.5)
    h_angle = angle_hist.weights
    h_bins = pi* (angle_hist.edges[1][2:end])

    #smoothing of histogram
    width = 25
    filt = (1/width)*ones(width)   # filter window
    ll = length(h_angle)
    h_angle1 = [h_angle; h_angle ; h_angle]  # wrap the histogram to consider periodic boundaries
    h_angle1 = filtfilt(filt,h_angle1)
    h_angle = h_angle1[ll+1:2*ll]
    
    # normalization of histogram
    h_angle = h_angle / sum(h_angle)

    h_angle, h_bins
end
