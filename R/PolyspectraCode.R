#' Calculate DFT
#'
#' @param x Time Series Vector
#' @param lambda Fourier Frequency
#'
#' @return Scalar Value, denoting the Discrete Fourier Transform of x at lambda
#' @export
#'
#' @examples


discrete.Fourier = function(x, lambda){
  s = 0
  for(i in 0:(length(x)-1)){
    z = complex(1, cos(lambda*i), -sin(lambda*i))
    s = s + x[i+1]*z
  }
  return(s)
}


#' Obtain an estimate of polyspectra
#'
#' @param x Time Series Vector
#' @param lambda Fourier Frequency. The kength of lambda indicate the order of
#' the polyspectra to be computed
#'
#' @return An estimate of the Polyspectra
#' @export
#'
#' @examples


polyspectra.Estimate = function(x, lambda){
  ll = length(lambda)
  pl = 1
  for(i in 1:ll){
    pl = pl*discrete.Fourier(x, lambda[i])
  }
  pl = pl*discrete.Fourier(x, - sum(lambda))
  #pl = pl/length(x)
  pl = pl*((2*pi)^(-ll))/length(x)
  return(pl)
}



#' Estimate of spectral mean
#'
#' @param x Time Series Vector
#' @param g Weight Function
#'
#' @return
#' @export
#'
#' @examples

spectral.Mean.Estimate = function(x, g, breaks){
  lambda = 2*pi*seq(-breaks, breaks, 1)/(2*breaks)
  s = 0
  for(i in 1:length(lambda)){
    s = s + polyspectra.Estimate(x, lambda[i])*g(lambda[i])
  }
  s = s / (length(x))
  return(Re(s))
}

#' Estimate of Bispectral mean
#'
#' @param x Time series vector
#' @param g Weight Function
#'
#' @return
#' @export
#'
#' @examples

bispectral.Mean.Estimate = function(x, g, breaks){
  lambda = 2*pi*seq(-breaks, breaks, 1)/(2*breaks)  s = 0
  sc = NULL
  for(i in 1:length(lambda)){
    for(j in 1:length(lambda)){
      f = polyspectra.Estimate(x, c(lambda[i], lambda[j]))*
        g(c(lambda[i], lambda[j]))
      s = s + f
      sc = c(sc, f)
    }
  }
  s = s / (length(x)^2)
  return(Re(s))
}


#' Theoretical spectral mean
#'
#' @param f Spectral Density
#' @param g Weight Function
#' @param breaks Granularity
#'
#' @return
#' @export
#'
#' @examples

spectral.Mean.Theory = function(f, g, breaks){
  lambda = 2*pi*seq(-floor(breaks/2), floor(breaks/2), 1)/breaks
  s = 0
  for(i in 1:length(lambda)){
    s = s + f(lambda[i])*g(lambda[i])
  }
  s = s / (breaks)
  return(Re(s))
}

#' Theoretical Bispectral mean
#'
#' @param f Bispectra
#' @param g Weight Function
#' @param breaks
#'
#' @return
#' @export
#'
#' @examples

bispectral.Mean.Theory = function(f, g, breaks){
  lambda = 2*pi*seq(-floor(breaks/2), floor(breaks/2), 1)/breaks
  s = 0
  for(i in 1:length(lambda)){
    for(j in 1:length(lambda)){
      s = s + f(c(lambda[i], lambda[j]))*g(c(lambda[i], lambda[j]))
    }
  }
  s = s / (breaks^2)
  return(Re(s))
}

#' Theoretical trispectral mean
#'
#' @param f Trispectra
#' @param g Weight Function
#' @param breaks
#'
#' @return
#' @export
#'
#' @examples
trispectral.Mean.Theory = function(f, g, breaks){
  s = 0
  lambda = 2*pi*seq(-floor(breaks/2), floor(breaks/2), 1)/breaks
  for(i in 1:length(lambda)){
    for(j in 1:length(lambda)){
      for(k in 1:length(lambda)){
        s = s + f(c(lambda[i], lambda[j], lambda[k]))*
          g(c(lambda[i], lambda[j], lambda[k]))
      }
    }
  }
  s = s / (breaks^3)
  return(Re(s))
}



#' Theoretical Tetraspectral mean
#'
#' @param f Tetraspectra
#' @param g Weight Function
#' @param breaks
#'
#' @return
#' @export
#'
#' @examples

tetraspectral.Mean.Theory = function(f, g, breaks){
  s = 0
  lambda = 2*pi*seq(-floor(breaks/2), floor(breaks/2), 1)/breaks
  for(i in 1:length(lambda)){
    for(j in 1:length(lambda)){
      for(k in 1:length(lambda)){
        for(l in 1:length(lambda)){
          s = s + f(c(lambda[i], lambda[j], lambda[k], lambda[l]))*
            g(c(lambda[i], lambda[j], lambda[k], lambda[l]))
        }
      }
    }
  }
  s = s / (breaks^4)
  return(Re(s))
}




#' Theoretical Asymptotic Variance of bispectral mean
#'
#' @param ps Spectra
#' @param bs Bispectra
#' @param ts Trispectra
#' @param qs Tetraspectra
#' @param g Weight Function
#' @param breaks Granularity
#'
#' @return
#' @export
#'
#' @examples

bispectralAsympVariance = function(ps, bs, ts, qs, g, breaks){
  #
  # if(type == "separable"){
  #   g1 = function(y){g(y[1])*g(y[2])}
  #   g2 = function(y){g(y[1])*g(y[2])*g(y[3])}
  #   y = 6*(3 * spectral.Mean.Theory(ps, g, breaks)*
  #            trispectral.Mean.Theory(ts, g2, breaks)
  #          + 3 * bispectral.Mean.Theory(bs, g1, breaks)^2 +
  #            spectral.Mean.Theory(ps,g, breaks)^3)
  # }else{
  lambda = 2*pi*seq(-floor(breaks/2), floor(breaks/2), 1)/breaks

  s1 = 0
  s2 = 0
  s3 = 0
  for(i in 1:length(lambda)){
    for(j in 1:length(lambda)){

      s3 = s3 + g(c(lambda[i], lambda[j]))*
        g(c(-lambda[i], -lambda[j]))*given_spectra(lambda[i])*
        given_spectra(lambda[j])*given_spectra(-lambda[i] - lambda[j])


      for(k in 1:length(lambda)){
        s1 = s1 + g(c(lambda[i], lambda[j]))*
          g(c(-lambda[i] - lambda[j], lambda[k]))*
          given_bispectra(c(lambda[i], lambda[j]))*
          given_bispectra(c(-lambda[i] - lambda[j], lambda[k]))


        s2 = s2 + g(c(lambda[i], lambda[j]))*
          g(c(-lambda[i] - lambda[j] - lambda[k], lambda[k]))*
          given_trispectra(c(lambda[i], lambda[j], lambda[k]))*
          given_spectra(-lambda[i] - lambda[j])

      }
    }
  }

  s1 = (s1/(breaks^3))*9*2
  s2 = (s2/(breaks^3))*9*2
  s3 = (s3/(breaks^2))*6

  y = s1 + s2 + s3



  # }


  return(Re(y))

}



