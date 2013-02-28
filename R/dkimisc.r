rotateKurtosis <- function( evec, KT, i = 1, j = 1, k = 1, l = 1) {
  
  if ( ( i < 1) | ( i > 3)) stop( "rotateKurtosis: index i out of range")
  if ( ( j < 1) | ( j > 3)) stop( "rotateKurtosis: index j out of range")
  if ( ( k < 1) | ( k > 3)) stop( "rotateKurtosis: index k out of range")  
  if ( ( l < 1) | ( l > 3)) stop( "rotateKurtosis: index l out of range")

  if ( length( dim( evec)) != 3) stop( "rotateKurtosis: dimension of direction array is not 3")
  if ( dim( evec)[1] != 3) stop( "rotateKurtosis: length of direction vector is not 3")
  if ( dim( evec)[2] != 3) stop( "rotateKurtosis: number of direction vectors is not 3")

  if ( length( dim( KT)) != 2) stop( "rotateKurtosis: dimension of kurtosis array is not 2")
  if ( dim( KT)[1] != 15) stop( "rotateKurtosis: kurtosis tensor does not have 15 elements")
  if ( dim( KT)[2] != dim( evec)[3]) stop( "rotateKurtosis: number of direction vectors does not match number of kurtosis tensors")
  
  nvox <- dim( KT)[2]
  
  ## we create the full symmetric kurtosis tensor for each voxel ...
  W <- defineKurtosisTensor( KT)
  ## ..., i.e., dim(W) == c( 3, 3, 3, 3, nvox) 
  
  Wtilde <- rep( 0, nvox)
  
  for ( ii in 1:3) {
    
    for ( jj in 1:3) {
      
      for ( kk in 1:3) {
        
        for ( ll in 1:3) {
          
          ## I am not sure about the order of i and ii, should it be evec[ i, ii, n]? Same for j, k, l
          Wtilde <- Wtilde + evec[ ii, i, ] * evec[ jj, j, ] * evec[ kk, k, ] * evec[ ll, l, ] * W[ ii, jj, kk, ll, ]   
          
        }
        
      }
      
    }
    
  }
  
  invisible( Wtilde)
  
}


defineKurtosisTensor <- function( DK) {
  
  W <- array( 0, dim = c( 3, 3, 3, 3, dim( DK)[2]))
  
  W[ 1, 1, 1, 1, ] <- DK[ 1, ]
  W[ 2, 2, 2, 2, ] <- DK[ 2, ]
  W[ 3, 3, 3, 3, ] <- DK[ 3, ]
  
  W[ 1, 1, 1, 2, ] <- DK[ 4, ]
  W[ 1, 1, 2, 1, ] <- DK[ 4, ]
  W[ 1, 2, 1, 1, ] <- DK[ 4, ]
  W[ 2, 1, 1, 1, ] <- DK[ 4, ]
  
  W[ 1, 1, 1, 3, ] <- DK[ 5, ]
  W[ 1, 1, 3, 1, ] <- DK[ 5, ]
  W[ 1, 3, 1, 1, ] <- DK[ 5, ]
  W[ 3, 1, 1, 1, ] <- DK[ 5, ]
  
  W[ 2, 2, 2, 1, ] <- DK[ 6, ]
  W[ 2, 2, 1, 2, ] <- DK[ 6, ]
  W[ 2, 1, 2, 2, ] <- DK[ 6, ]
  W[ 1, 2, 2, 2, ] <- DK[ 6, ]
  
  W[ 2, 2, 2, 3, ] <- DK[ 7, ]
  W[ 2, 2, 3, 2, ] <- DK[ 7, ]
  W[ 2, 3, 2, 2, ] <- DK[ 7, ]
  W[ 3, 2, 2, 2, ] <- DK[ 7, ]
  
  W[ 3, 3, 3, 1, ] <- DK[ 8, ]
  W[ 3, 3, 1, 3, ] <- DK[ 8, ]
  W[ 3, 1, 3, 3, ] <- DK[ 8, ]
  W[ 1, 3, 3, 3, ] <- DK[ 8, ]
  
  W[ 3, 3, 3, 2, ] <- DK[ 9, ]
  W[ 3, 3, 2, 3, ] <- DK[ 9, ]
  W[ 3, 2, 3, 3, ] <- DK[ 9, ]
  W[ 2, 3, 3, 3, ] <- DK[ 9, ]
  
  W[ 1, 1, 2, 2, ] <- DK[ 10, ]
  W[ 1, 2, 1, 2, ] <- DK[ 10, ]
  W[ 1, 2, 2, 1, ] <- DK[ 10, ]
  W[ 2, 1, 2, 1, ] <- DK[ 10, ]
  W[ 2, 1, 1, 2, ] <- DK[ 10, ]
  W[ 2, 2, 1, 1, ] <- DK[ 10, ]
  
  W[ 1, 1, 3, 3, ] <- DK[ 11, ]
  W[ 1, 3, 1, 3, ] <- DK[ 11, ]
  W[ 1, 3, 3, 1, ] <- DK[ 11, ]
  W[ 3, 1, 3, 1, ] <- DK[ 11, ]
  W[ 3, 1, 1, 3, ] <- DK[ 11, ]
  W[ 3, 3, 1, 1, ] <- DK[ 11, ]
  
  W[ 2, 2, 3, 3, ] <- DK[ 12, ]
  W[ 2, 3, 2, 3, ] <- DK[ 12, ]
  W[ 2, 3, 3, 2, ] <- DK[ 12, ]
  W[ 3, 2, 3, 2, ] <- DK[ 12, ]
  W[ 3, 2, 2, 3, ] <- DK[ 12, ]
  W[ 3, 3, 2, 2, ] <- DK[ 12, ]
  
  W[ 1, 1, 2, 3, ] <- DK[ 13, ]
  W[ 1, 1, 3, 2, ] <- DK[ 13, ]
  W[ 3, 1, 1, 2, ] <- DK[ 13, ]
  W[ 2, 1, 1, 3, ] <- DK[ 13, ]
  W[ 2, 3, 1, 1, ] <- DK[ 13, ]
  W[ 3, 2, 1, 1, ] <- DK[ 13, ]
  W[ 1, 3, 1, 2, ] <- DK[ 13, ]
  W[ 1, 2, 1, 3, ] <- DK[ 13, ]
  W[ 3, 1, 2, 1, ] <- DK[ 13, ]
  W[ 2, 1, 3, 1, ] <- DK[ 13, ]
  W[ 1, 2, 3, 1, ] <- DK[ 13, ]
  W[ 1, 3, 2, 1, ] <- DK[ 13, ]
  
  W[ 2, 2, 1, 3, ] <- DK[ 14, ]
  W[ 2, 2, 3, 1, ] <- DK[ 14, ]
  W[ 3, 2, 2, 1, ] <- DK[ 14, ]
  W[ 1, 2, 2, 3, ] <- DK[ 14, ]
  W[ 1, 3, 2, 2, ] <- DK[ 14, ]
  W[ 3, 1, 2, 2, ] <- DK[ 14, ]
  W[ 3, 2, 1, 2, ] <- DK[ 14, ]
  W[ 2, 3, 2, 1, ] <- DK[ 14, ]
  W[ 1, 2, 3, 2, ] <- DK[ 14, ]
  W[ 2, 1, 2, 3, ] <- DK[ 14, ]
  W[ 2, 1, 3, 2, ] <- DK[ 14, ]
  W[ 2, 3, 1, 2, ] <- DK[ 14, ]
  
  W[ 3, 3, 1, 2, ] <- DK[ 15, ]
  W[ 3, 3, 2, 1, ] <- DK[ 15, ]
  W[ 2, 3, 3, 1, ] <- DK[ 15, ]
  W[ 1, 3, 3, 2, ] <- DK[ 15, ]
  W[ 1, 2, 3, 3, ] <- DK[ 15, ]
  W[ 2, 1, 3, 3, ] <- DK[ 15, ]
  W[ 1, 3, 2, 3, ] <- DK[ 15, ]
  W[ 2, 3, 1, 3, ] <- DK[ 15, ]
  W[ 3, 2, 3, 1, ] <- DK[ 15, ]
  W[ 3, 1, 3, 2, ] <- DK[ 15, ]
  W[ 3, 1, 2, 3, ] <- DK[ 15, ]
  W[ 3, 2, 1, 3, ] <- DK[ 15, ]
  
  invisible( W)
}


## for efficiency we should combine both functions
kurtosisFunctionF1 <- function( l1, l2, l3) {
  
  require( gsl)
  ## Tabesh et al. Eq. [28]
  ## this function is defined without MD^2!!
  ( ellint_RF( l1/l2, l1/l3, 1) * sqrt( l2*l3) / l1 + ellint_RD( l1/l2, l1/l3, 1) * ( 3* l1^2 - l1*l2 - l1*l3 - l2*l3) / (3*l1*sqrt(l2*l3)) - 1) / 2 / ( l1-l2) / ( l1-l3)

  ## consider removable singularities!!
  # ind1: ((l1 == l2) | (l1 == l3)) & !(l2 == l3))
  # kurtosisFunctionF2( l2, l1, l1) / 2 
  # kurtosisFunctionF2( l3, l1, l1) / 2 

  # ind2: (l1 == l2 == l3)
  # 1/5
}

kurtosisFunctionF2 <- function( l1, l2, l3) {
  
  require( gsl)
  ## Tabesh et al. Eq. [28]
  ## this function is defined without MD^2!!
  3 * ( ellint_RF( l1/l2, l1/l3, 1) * (l2+l3) / sqrt(l2*l3) + ellint_RD( l1/l2, l1/l3, 1) * (2*l1-l2-l3) / 3/sqrt(l2*l3) - 2 ) / (l2 - l3) / (l2 - l3)

  ## consider removable singularities!!
  # ind1: (!((l1 == l2) | (l1 == l3))) & (l2 == l3))
  # alpha(x) = 1/sqrt(abs(x)) * atan(sqrt(abs(x)))
  # 6 (l1+2*l3)^2/144/l3^2/(l1-l3)^2 *(l3*(l1+2*l3)+ l1*(l1-4*l3)*alpha(1-l1/l3))
  
  # ind2: (l1 == l2 == l3)
  # 6/15
  
}

pseudoinverseSVD <- function( xxx) {
  svdresult <- svd( xxx)
  svdresult$v %*% diag( 1 / svdresult$d) %*% t( svdresult$u)
}