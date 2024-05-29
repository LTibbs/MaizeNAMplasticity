# Genomic prediction functions --------------------------------------------

# Laura Tibbs Cortes
# June 9, 2020

# make FUNCTIONS to do the repetitive parts of genomic prediction (not just giant loops anymore)


# function to split observations into folds for kfold cross-validation:
fold.maker <- function(num.obs, kfold) { # provide the total number of observations you are splitting and the number of folds
  new.order <- sample(1:num.obs, num.obs) # randomly permute observation order
  Q <- num.obs%/%kfold # find whole-number quotient (number of obs in some folds)
  R <- num.obs%%kfold # find remainder (how many folds need an extra obs)
  
  folds <- vector("list", length = kfold) # make list to hold folds
  for (i in 1:kfold) {
    if (i <=R) { # for the larger folds
      # print(((i-1)*(Q+1)+1):(i*(Q+1)))
      folds[[i]] <- new.order[c(((i-1)*(Q+1)+1):(i*(Q+1)))]
    }
    else if (i>R) { # for the smaller folds
      # print(((R*(Q+1))+(i-R-1)*Q+1):((R*(Q+1))+(i-R)*Q))
      folds[[i]] <- new.order[c((R*(Q+1))+(i-R-1)*Q+1):((R*(Q+1))+(i-R)*Q)] # identify the next fold
    } else {warning("Error in fold creation.")}
  }

return(folds)
}

# make folds given observation ids rather than number of obs
fold.maker.2 <- function(obs.ids, kfold) { # provide the ids of observations you are splitting
  new.order <- sample(obs.ids, length(obs.ids))
  Q <- length(obs.ids)%/%kfold # find whole-number quotient (number of obs in some folds)
  R <- length(obs.ids)%%kfold # find remainder (how many folds need an extra obs)  
  
  folds <- vector("list", length = kfold) # make list to hold folds
  for (i in 1:kfold) {
    if (i <=R) { # for the larger folds
      # print(((i-1)*(Q+1)+1):(i*(Q+1)))
      folds[[i]] <- new.order[c(((i-1)*(Q+1)+1):(i*(Q+1)))]
    }
    else if (i>R) { # for the smaller folds
      # print(((R*(Q+1))+(i-R-1)*Q+1):((R*(Q+1))+(i-R)*Q))
      folds[[i]] <- new.order[c((R*(Q+1))+(i-R-1)*Q+1):((R*(Q+1))+(i-R)*Q)] # identify the next fold
    } else {warning("Error in fold creation.")}
  }

return(folds)
}
