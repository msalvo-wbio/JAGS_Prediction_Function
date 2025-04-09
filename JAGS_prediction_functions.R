#Creator: Marco Salvo
#Last Updated: 04/01/2025
#R Version Updated on: 4.4.2

library(reshape2)

# The function for predictions
jags.predict <- function(intercepts,
                         slope.coefficients,
                         exp.cov.values,
                         credible.intervals = c(0.025, 0.5, 0.975),
                         res = 100,
                         link.function,
                         lognormal.sd,
                         covariate.type = "continuous"
){
  #extract the number of iterations from your simulated intercepts
  iters <- length(intercepts)
  
  ### Setup Data:
  
  #contin
  if(covariate.type == "continuous"){
    #create sequences of your covariate values
    for(ixi in 1:length(exp.cov.values)){ #for every covariate value
      if("numeric" %in% class(exp.cov.values[[ixi]])){ # if it is numeric:
        if(length(exp.cov.values[[ixi]]) == 1){ #if expected value is length 1, repeat it up the times of res
          exp.cov.values[[ixi]] <- rep(exp.cov.values[[ixi]], times = res)
        } else if(length(exp.cov.values[[ixi]]) == 2){ #if expected values are length 2 (a range), create a sequence of times of res
          exp.cov.values[[ixi]] <- seq(exp.cov.values[[ixi]][1], 
                                       exp.cov.values[[ixi]][2],
                                       length.out = res)
        } else { #if a covariate value is not length 1 or 2, give an error
          stop("Expected covariate values must be of length 1 or 2")
        }
        
      }
    }
  } else if(covariate.type == "categorical"){
    
    #check if the name of each list item is "category"
    is.category <- sapply(exp.cov.values, function(x) x == "category")
    
    #for every covariate value
    for(ixi in 1:length(exp.cov.values)){ 
      if(is.category[ixi] == T){ #if it is categorical
        exp.cov.values[[ixi]] <- rep(0, times = length(exp.cov.values[is.category]) + 1) #create a sequence of 0's which is 1+length of the number of categories
      } else if(is.category[ixi] == F){ #else if it's continuous
        if(length(exp.cov.values[[ixi]]) == 1){
          exp.cov.values[[ixi]] <- rep(exp.cov.values[[ixi]], 
                                       times = length(exp.cov.values[is.category]) + 1) #repeat the means the same number of times as the number of categories
        } else { #no varying values within a categorical prediction (yet)
          stop("Current version does not account for interactions, set continuous covariates to their mean value")
        }
      }
    }
    
    #Fill in 1's in every 0 sequence
    for(i in seq_along(exp.cov.values)){
      if(i == 1){pos <- 2} #set position index to 2 for the first value
      
      if (all(exp.cov.values[[i]] == 0)) {  # Check if the vector is all zeros
        if(i != 1){pos <- pos + 1} #add a 1 to the position index
        exp.cov.values[[i]][pos] <- 1 #set the position value to 1
      }
    }
    
    #set res
    res <- length(exp.cov.values[is.category]) + 1
    
  } else {stop("covariate.type must be one of 'continuous' or 'categorical'.")}
  
  
  #create a matrix for:
  pred.mat <- matrix(NA, iters, res) #predicted values
  pred.cred <- matrix(NA, res, length(credible.intervals)) #predicted median and credible intervals
  
  # Get your estimates
  for(jlv in 1:res){ # for all resolutions:
    
    parameter.values <- list() #initialize list to store parameter values in
    for(ixi in 1:length(exp.cov.values)){ #for all covariates in the model:
      
      if("numeric" %in% class(exp.cov.values[[ixi]])){ #different indexing for vectors and matrices
        parameter.values[[ixi]] <- slope.coefficients[[ixi]] * exp.cov.values[[ixi]][jlv]#get the product of the slope and the covariate value
      } else if("matrix" %in% class(exp.cov.values[[ixi]])){
        parameter.values[[ixi]] <- slope.coefficients[[ixi]] * exp.cov.values[[ixi]][,jlv]#get the product of the slope and the covariate value
      } else {
        stop("expected values must be a numeric vector or a expected.response.matrix matrix derived from this function")
      }
      
      
      if(ixi == 1){ #for the first value, add the intercepts and the first covariate term
        pred.mat[,jlv] <- intercepts + parameter.values[[ixi]]
      } else { #for all other values, add them to the current matrix column to get the complete expected value based on a normal distribution
        pred.mat[,jlv] <- pred.mat[,jlv] + parameter.values[[ixi]]
      } 
      
    }
    
    #If the model is not following a normal distribution:
    if(link.function == "log"){ #for log links
      pred.mat[,jlv] <- exp(pred.mat[,jlv])
    } else if(link.function == "lognormal"){ #for log-normal links
      pred.mat[,jlv] <- rlnorm(iters, pred.mat[,jlv], lognormal.sd)
    } else if(link.function == "logit"){ # for logit link
      pred.mat[,jlv] <- exp(pred.mat[,jlv])/(1+exp(pred.mat[,jlv]))
    }
    
    #get predicted credible intervals
    pred.cred[jlv,] <- quantile(pred.mat[,jlv], credible.intervals)
  }
  
  #melted observations for graphing
  melted.obs <- melt(pred.mat); names(melted.obs) <- c('iteration','cov','response')
  
  return(list(expected.cov.values = exp.cov.values,
              expected.response.matrix = pred.mat,
              expected.quantiles = pred.cred,
              plotting.data = melted.obs))
  
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# The function for plotting predictions

plot.jags.predictions <- 
  function(plotting.data, 
           expected.covariate.values,
           quantiles,
           x.lab,
           y.lab,
           covariate.type = "continuous",
           category.names){
    
    if(covariate.type == "continuous"){
      
      smoothScatter(plotting.data$response ~ expected.covariate.values[plotting.data$cov],
                    las = 1, nrpoints = 0,
                    ylab = y.lab,
                    xlab = x.lab)
      lines(quantiles[,2] ~ expected.covariate.values, lty = 1, lwd = 3, col = 'white')
      lines(quantiles[,1] ~ expected.covariate.values, lty = 2, lwd = 3, col = 'white')
      lines(quantiles[,3] ~ expected.covariate.values, lty = 2, lwd = 3, col = 'white')
      
    } else if(covariate.type == "categorical"){
      
      plot.df <- 
        data.frame(category = category.names,
                   cov = 1:nrow(quantiles),
                   low.CI = quantiles[,1],
                   median = quantiles[,2],
                   high.CI = quantiles[,3])
      
      ggplot() + 
        geom_violin(data = plotting.data, aes(as.factor(cov), response), bw = 0.01, fill = "gray90") + 
        geom_linerange(data = plot.df, aes(x = as.factor(cov), ymin = low.CI, ymax = high.CI), size = 0.7) +
        geom_point(data = plot.df, aes(as.factor(cov), median), color = "red", size = 2.5) +
        theme_bw() +
        scale_x_discrete(name = waiver(), labels = category.names) +
        ylab("Probability of Selection") +
        xlab("")
      
    } else {
      stop("covariate.type must be one of 'continuous' or 'categorical'")
    }
    
  }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
