#'tree_logL
#'
#'Calculate the loglikelihood of a Cox-Ingersoll-Ross Diffusion Process in the Tree of Life
#'
#'  This function caluclates the log likelihood.
#'
#'@param tr An object of class \code{phylo} from the ape package.  In this version the CIR
#'process parameters alpha, mu, and sigma for each branch are within vectors in the same order
#'as the edge (branch) labelling.
#'@param tipdata A numeric vector containing tip values in the same order as the tip labels in \code{tr$tip.label}.
#'@param lst A list of time series projects for each simulated path equal to the length of
#'the number of branches in the \code{tr} object.
#'@param alpha Set to NULL if alpha is to be estimated,
#'otherwise set to a numeric value or a numeric vector specifying
#'the value of the parameter for all the branches/edges.
#'In the latter case, the values must be specifying in the same order
#'as the edges in the \code{tr} object.
#'@param mu As \code{alpha}.
#'@param sigma As \code{alpha}.
#'@param model A list containing drift, diffusion, and the partial differentiation of diffusion as quoted
#'expressions using method quote. For the Euler scheme
#'the drift coefficient as \code{drift}, the diffusion coefficient as
#'\code{diffusion}, and the partial differentiation of
#'\code{diffusion} by \code{x} as \code{dx_diffusion} is required.
#'@param method Specified as either "euler" or "milstein."
#'@param ... Not used.
#'
#'@return logL, an integer.
#'
#'@export


tree_logL <- function(tr, tipdata, lst, alpha, mu, sigma, model,
              method, ...) {

tipdata <- as.numeric(tipdata)
n_tips  <- length(tr$tip.label)
rt_node <- n_tips + 1
logL <- NULL
rt_node_dist <- ape::dist.nodes(tr)[rt_node, ]

logL_edges <- function (node, tr, tipdata, lst, alpha, mu, sigma, model) {
  #The daughter nodes are those that at first branch directily from rt_node, node is updated for subsequent splits
  daughters <- tr$edge[which(tr$edge[, 1] == node), 2]
  
  for (ind_d in 1:2) {
    edge <- which((tr$edge[, 1] == node) & (tr$edge[, 2] == daughters[ind_d]))
    theta <- c(alpha[edge], mu[edge], sigma[edge])

    if (daughters[ind_d] > n_tips) {
      logL[edge] <<- logl_fn(X = lst[[edge]], theta = theta,
                            model = model, log = TRUE,
                            method = method)
     
    } else {
      logL[edge] <<- logl_fn(X = lst[[edge]], theta = theta,
                          model = model, log = TRUE, method = method) +
                     dc_fn(x = tipdata[daughters[ind_d]],
                          t = rt_node_dist[daughters[ind_d]],
                          x0 = lst[[edge]][length(lst[[edge]])],
                          t0 = tsp(lst[[edge]])[2],
                          theta = theta,
                          model = model,
                          log = TRUE,
                          method = method)
    }
    
    
    
  
    # Reset the node to the "new root" 
    if (daughters[ind_d] > n_tips) {
      logL_edges(daughters[ind_d], tr, tipdata, lst, alpha, mu, sigma, model)
    }
  }
}

logL_edges(rt_node, tr, tipdata, lst, alpha, mu, sigma, model)
# Log Likelihood sums the likelihoods rather than taking the product which is more efficient
return(sum(logL))
}
