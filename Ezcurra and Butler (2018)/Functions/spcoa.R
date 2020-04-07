spcoa <- function (D,
                   Correction = "none",
                   DoStress = FALSE,
                   DoR2likeratio = FALSE) {
  
  ### This function performs Primary Coordinates Analysis on a provided distance
  ### matrix. It can also correct for negative eigenvalues and show the stress
  ### and R^2-like ratio of the analysis.
  
  ### SET UP VARIABLES.
  
  # Check the proper class of distance matrix is provided.
  if (!is(D, "dist")) {
    stop("A distance matrix of class 'dist' must be provided.")
  }
  
  # Transform the distance matrix into a proper matrix.
  D <- as.matrix(D)
  
  # Get the number of rows of the distance matrix.
  n <- nrow(D)
  
  # Get the threshold for accurate results.
  epsilon <- sqrt(.Machine$double.eps)
  
  # Get the index for the correction.
  CORRECTION <- c("none", "lingoes", "cailliez")
  Correction <- pmatch(x = Correction, table = CORRECTION, nomatch = NA)
  
  # Stop function if method is not correctly specified.
  if (is.na(Correction)) {
    stop("Invalid correction method.")
  }
  
  ### UNCORRECTED PCoA.
  
  # Make centering matrix 'Cn'.
  Cn <- diag(x = 1, nrow = n, ncol = n) -
    (matrix(data = 1, nrow = n, ncol = n) / n)
  
  # Create double-centered matrix 'B'.
  B <- -0.5 * (Cn %*% (D^2) %*% Cn)
  
  # Obtain the eigenvalues and eigenvectors of B.
  E <- eigen(B)
  
  # Transform values smaller than the threshold into proper 0s.
  E$values[abs(E$values) < epsilon] <- 0
  E$vectors[abs(E$vectors) < epsilon] <- 0
  
  ### IS CORRECTION NECESSARY?
  
  if (any(E$values < 0)) {
    
    ## No correction.
    
    if (Correction == 1) {
      
      # No correction note.
      COR <- "Negative eigenvalues found. No correction was applied."
      
    }
    
    ## Lingoes correction.
    
    if (Correction == 2) {
      
      # Set the correction value as the absolute value of the min eigenvalue.
      c1 <- abs(min(E$values))
      # Create a new B matrix, applying the correction.
      B <- (-0.5 * (D^2)) - c1
      # Restore the diagonal.
      diag(B) <- 0
      # Center the B matrix.
      B <- Cn %*% B %*% Cn
      
      # Lingoes note.
      COR <- list(method = "Negative eigenvalues found. Lingoes correction applied.",
                  value = c1)
      
    }
    
    ## Cailliez correction.
    
    if (Correction == 3) {
      
      # Create the delta2 matrix.
      delta2 <- Cn %*% (-0.5 * D) %*% Cn
      # Create the special matrix to calculate c2.
      c2matrix <- rbind(cbind(matrix(data = 0, nrow = n, ncol = n),
                              diag(x = -1, nrow = n, ncol = n)),
                        cbind(2 * B, -4 * delta2))
      # Calculate c2.
      c2 <- max(Re(eigen(x = c2matrix,
                         symmetric = FALSE,
                         only.values = TRUE)$values))
      # Calculate the new B matrix.
      B <- -0.5 * (D + c2)^2
      # Restore the diagonal.
      diag(B) <- 0
      # Centre the matrix.
      B <- Cn %*% B %*% Cn
      
      # Cailliez note.
      COR <- list(method = "Negative eigenvalues found. Cailliez correction applied.",
                  value = c2)
      
    }
    
    ## Recalculate the eigenvalues.
    
    # Obtain the eigenvalues and eigenvectors of B.
    E <- eigen(B)
    
    # Keep the real part if any complex numbers are present.
    E$values <- Re(E$values)
    E$vectors <- Re(E$vectors)
    
    # Transform values smaller than the threshold into proper 0s.
    E$values[abs(E$values) < epsilon] <- 0
    E$vectors[abs(E$vectors) < epsilon] <- 0
    
    if (any(E$values < 0)) {
      
      warning("Negative eigenvalues found!")
      
    }
    
  } else {
    
    COR <- "No negative eigenvalues were found. No correction applied."
    
  }
  
  # Get a vector with the indeces of the eigenvalues greater than 0.
  ind <- 1:sum(E$values > 0)
  
  # Calculate the final vectors by multiplying the eigenvectors
  # by the sqrt of the eigenvalues (in a diagonal matrix).
  X <- E$vectors[, ind] %*% sqrt(diag(E$values[ind]))
  
  rownames(X) <- rownames(D)
  
  # Calculate the variance explained per dimension.
  exvar <- E$values[ind] / sum(E$values[ind])
  # Give names to this vector.
  names(exvar) <- paste("dim", 1:ncol(X), sep = "_")
  
  # Calculate cumulative variance.
  cumvar <- cumsum(exvar)
  # Give names to this vector.
  names(cumvar) <- paste(1:ncol(X), "dims", sep = "_")
  
  ### STRESS CALCULATION.
  
  # If stress in wanted...
  if (DoStress) {
    
    # Set up a vector.
    stress <- rep(0, ncol(X))
    
    # For 1 to all available dimensions...
    for (i in 1:ncol(X)) {
      
      # Create a distance matrix for the new vectors with the corresponding dimensions.
      Y <- as.matrix(dist(X[, 1:i], "euc"))
      
      # Set up a vector 'v' for the squared differences between the original distances
      # and the results of the PCoA.
      v <- rep(0, (n - 1))
      
      # For every vector...
      for (j in 1:(n - 1)) {
        
        # Calculate the Dij - PCoA(Dij)
        v[j] <-  D[j, (j + 1)] - Y[j, (j + 1)]
        
      }
      
      # Transform values below threshold into proper 0s.
      v[v < epsilon] <- 0L
      
      # Square 'em.
      v <- v ^ 2
      
      # Calculate the stress for the corresponding dimensions.
      stress[i] <- sqrt(sum(v) / (sum(Y ^ 2)))
      
    }
    
    # Set the names for the vector.
    names(stress) <- paste(1:ncol(X), "dims", sep = "_")
    
  } else {
    
    # If stress is not required, report it.
    stress <- "Stress was not required."
    
  }
  
  ### R2 LIKE RATIO CALCULATION.
  
  # If R2 like ratio is wanted...
  if (DoR2likeratio) {
    
    # Set up a vector for the values.
    R2lr <- rep(0, ncol(X))
    
    # For every number of dimensions.
    for (i in 1:ncol(X)) {
      
      # Calculate the R2 like ratio.
      R2lr[i] <- sum(E$values[1:i]) / sum(E$values)
      
    }
    
    # Set the names for the vectors.
    names(R2lr) <- paste(1:ncol(X), "dims", sep = "_")
    
  } else {
    
    R2lr <- "R^2-like ratio was not required."
    
  }
  
  ### RETURN.
  
  return(list("correction" = COR,
              "results" = list("values" = E$values[ind],
                               "vectors" = X,
                               "all_values" = E$values),
              "quality" = list("cumvar" = cumvar,
                               "expvar" = exvar,
                               "stress" = stress,
                               "R2" = R2lr)))
  
}
