ComparableTaxa <- function (CladisticMatrix, Silent = FALSE) {
  
  require("Claddis")
  
  n_matrices <- length(CladisticMatrix) - 1
  
  res <- vector(mode = "list", length = n_matrices)
  
  perc <- vector(mode = "numeric", length = n_matrices)
  
  names(res) <- names(CladisticMatrix[-1])
  
  for (i in (seq_len(n_matrices) + 1)) {
    
    mtx <- CladisticMatrix[[i]]$Matrix
    
    logmtx <- matrix(data = TRUE,
                     nrow = nrow(mtx),
                     ncol = ncol(mtx),
                     dimnames = dimnames(mtx))
    
    logmtx[which(is.na(mtx) | mtx == "")] <- FALSE
    
    combs <- combn(x = nrow(mtx), m = 2)
    
    comp_mtx <- matrix(NA, nrow = ncol(combs), ncol = ncol(mtx))
    
    for (j in seq_len(ncol(combs))) {
      
      comp_mtx[j, ] <- logmtx[combs[1, j], ] & logmtx[combs[2, j], ]
      
    }
    
    res[[i - 1]] <- sum(rowSums(comp_mtx) == 0)
    
    perc[i - 1] <- round((res[[i - 1]] / nrow(comp_mtx)) * 100, digits = 2)
    
  }
  
  if (!Silent) {
    
    if (n_matrices == 1) {
      
      cat(n_matrices, "matrix has been found.\n\n")
      
    } else {
      
      cat(n_matrices, "matrices have been found.\n\n")
      
    }
    
    for (i in seq_len(n_matrices)) {
      
      if (res[[i]] == 0) {
        
        cat(paste0(names(res)[i], ":"),
            "All taxon pairs are comparable for at least one character.")
        
      } else {
        
        cat(paste0(names(res)[i], ":"),
            res[[i]], "taxon pairs cannot be compared for any character",
            paste0("(", perc[i], "%)."))
        
      }
      
    }
    
  }
  
  return(invisible(res))
  
}