# Operations for GF(s) with s being prime number
add_GFs <- function(x, y, s){
  return((x + y)%%s)
}

minus_GFs <- function(x, y, s){
  return((x - y)%%s)
}

multiply_GFs <-function(x, y, s){
  return((x * y)%%s)
}

divide_GFs <- function(x, y, s){
  
  for (j in 1:(s - 1)) {
    if(multiply_GFs(j, y, s) == 1){
      return(multiply_GFs(x, j, s))
    }
  }
}

# Perform polynomial remainder on GF (s)
polynomial_mod <- function(poly_a, poly_b, s) {
  #poly_a and poly_b are coefficient vectors of the polynomials 
  if (all(poly_b == 0)) {
    stop("Modulus polynomial cannot be zero")
  }
  
  remainder <- poly_a
  divisor <- poly_b
  
  while (any(remainder != 0) && length(remainder) >= length(divisor)) {
    quotient_leading_term <- divide_GFs(remainder[1], divisor[1], s)  
    
    divisor_shifted <- c(divisor, rep(0, length(remainder) - length(divisor)))
    product <- (quotient_leading_term * divisor_shifted) %% s
    
    remainder <- (remainder - product) %% s
    
    first_nonzero_index <- min(which(remainder != 0))
    
    remainder <- remainder[first_nonzero_index:length(remainder)]
  }
  
  return(remainder)
}


# A function to define the multiplications on galois fied 
multiply_GF <- function(a, b, poly_primitive, s){
  # a and b are coefficient vectors of the polynomials as elements in GF(s^p)
  # poly_primitive is the coefficient vector of the primitive polymonial
  if(length(poly_primitive) == 1){
    return(multiply_GFs(a, b, s))
  }
  
  a <- as.array(a)
  b <- as.array(b)
  n <- length(a) + length(b) - 1
  c <- numeric(n)
  
  for (i in seq_along(a)) {
    
    for (j in seq_along(b)) {
      c[i + j - 1] <- c[i + j - 1] + a[i] * b[j]
    }
  }
  c <- c %% s
  if(all(c == 0)){
    c <- 0
  }else{
    nonzero_indices <- which(c != 0)
    c <- c[nonzero_indices[1]:length(c)]
  }
  reminder <- polynomial_mod(c, poly_primitive, s)
  
  return(reminder)
}

add_GF <- function(a,b,s){
  longest <- max(c(length(a),length(b)))
  if( length(a) < longest){
    a <- c(rep(0,longest-length(a)),a)
  }
  if( length(b) < longest){
    b <- c(rep(0,longest-length(b)),b)
  }
  c <- (a + b) %% s
  
  return(c)
}



# A function to construct the GSOA(s^p,K*(s^p-1)/(s-1),s^q,t), which is a part of mutiplication table
GSOA_opt <- function(polyprimitive, s, K = (s - 1), q = (length(polyprimitive)-1), add = 0){
  # poly_primitive is the coefficient vector of the primitive polynomial
  # 1 <= K <= s-1, q <= degree(primitive polynomial)
  # add is a shift element, represented by a vector 
  p <- length(polyprimitive) - 1
  z <- rev(c(s^(0:(q-1)),rep(0,p-q)))
  G <- matrix(0, nrow = s^p, ncol = (s^p - 1)/(s - 1))
  result <- NULL
  M <- as.matrix(expand.grid(rep(list(0:(s - 1)),p)))

  for (k in 1:K) {
    A <- NULL
    
    for (i in 1:nrow(M)) {
      nonzero_indices <- which(M[i, ] != 0)
      if (length(nonzero_indices) > 0) {
        first_nonzero_value <- M[i, nonzero_indices[1]]
        if (first_nonzero_value == k) {
          A <- rbind(A, M[i, ])
        }
      }
    }

    for (i in 1:nrow(M)) {
      x <- M[i, ]
      nonzero_indices <- which(x != 0)
      if (length(nonzero_indices) > 0){
        x <- x[nonzero_indices[1]:length(x)]
      }else{
        x <- 0
      }
      
      for (j in 1:nrow(A)) {
        y <- as.array(A[j, ])
        nonzero_indices <- which(y != 0)
        y <- y[nonzero_indices[1]:length(y)]
        a <- multiply_GF(x, y, polyprimitive, s)
        a <- add_GF(a,add,s)
        a <- c(rep(0, p - length(a)), a)
        G[i, j] <- sum(a * z)
      }
    }
    
    result <- cbind(result, G)
  }

  return(result)
}
