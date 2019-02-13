#'phylo_sde
#'
#'Simulate a Cox-Ingersoll-Ross Diffusion Process in the Tree of Life
#'
#'@param tr An object of class \code{phylo} from the ape package.  In this version the CIR
#'process parameters alpha, mu, and sigma for each branch are within vectors in the same order
#'as the edge (branch) labelling.
#'@param rt_value Value at the root of \code{tr}.
#'@param N Data imputation frequency.
#'@param theta Matrix of parameter values for each edge of the tree.
#'@param model A list containing drift, diffusion, and the partial differentiation of diffusion as quoted
#'expressions using method quote. For the Euler scheme
#'the drift coefficient as \code{drift}, the diffusion coefficient as
#'\code{diffusion}, and the partial differentiation of
#'\code{diffusion} by \code{x} as \code{dx_diffusion} is required.
#'See the Examples.
#'@param method Specified as either "euler" or "milstein."
#'@param ... Not used.
#'
#'@return A list of time series projects for each simulated path equal to the length of
#'the number of branches in the \code{tr} object.
#'
#'@export


phylo_sde <- function(tr, rt_value, N, theta, model, method, ...) {

  lst <- list()
  n_tips <- length(tr$tip.label)
  rt_node <- n_tips + 1

  #only relevant if specified
  dotslist <- list(...)
  if ("pred.corr" %in% names(dotslist)) {
    pred.corr <- dotslist$pred.corr
  } else {
    pred.corr <- FALSE
  }

  if (method == "milstein") {
    pred.corr <- TRUE
  }

  if (pred.corr) {
    if (!exists("dx_diffusion", model)) {
      model$dx_diffusion <- D(model$diffusion, "x")
    }
  } else {
      model$dx_diffusion <- quote(NULL)
  }

  sde_edges <- function(tr, node, X0, t0) {
    
    # Preceeding nodes
    daughters <- tr$edge[which(tr$edge[, 1] == node), 2]
    for (d_ind in 1:2) {
      
      # Traverse tree through daughter nodes
      edge <- which((tr$edge[, 1] == node) & (tr$edge[, 2] == daughters[d_ind]))
      
     
      drift <- as.expression(force(eval(substitute(substitute(e,
                              list(alpha = theta[edge, "alpha"],
                                   mu = theta[edge, "mu"],
                                   sigma = theta[edge, "sigma"])),
                              list(e = model$drift)))))
      
     
      diffusion <- as.expression(force(eval(substitute(substitute(e,
                              list(alpha = theta[edge, "alpha"],
                                   mu = theta[edge, "mu"],
                                   sigma = theta[edge, "sigma"])),
                              list(e = model$diffusion)))))
      
        
      diffusion_x <- as.expression(force(eval(substitute(substitute(e,
                              list(alpha = theta[edge, "alpha"],
                                   mu = theta[edge, "mu"],
                                   sigma = theta[edge, "sigma"])),
                              list(e = model$dx_diffusion)))))
     
      n_steps <- tr$edge.length[edge] * N
      
      tE <- t0 + tr$edge.length[edge]
      
      lst[[edge]] <<- sde::sde.sim(X0 = X0, t0 = t0, T = tE, N = n_steps,
                                   method = method,
                                   drift = drift,
                                   sigma = diffusion,
                                   sigma.x = diffusion_x,
                                   pred.corr = pred.corr)
      tE <- tsp(lst[[edge]])[2]
      if (daughters[d_ind] > n_tips) {
        sde_edges(tr, daughters[d_ind], lst[[edge]][n_steps + 1], tE)
      }
    }
  }
  sde_edges(tr, rt_node, X0 = rt_value, t0 = 0)

  # Remove tip values (we have observed tip values)
  node_len <- ape::node.depth.edgelength(tr)
  for (nthtip in 1:n_tips) {
    nEdge <- ape::which.edge(tr, tr$tip.label[nthtip])
    ntsp <- tsp(lst[[nEdge]])
    # If the end time for the node edges is less than T - 1/N,
    # the last sample is removed from the simulated data,
    # if there is more than one sample.
    if (ntsp[2] > (node_len[nthtip] - 1 / N)) {
      if (length(lst[[nEdge]]) > 1) {
        if (length(lst[[nEdge]]) == 2) {
          lst[[nEdge]] <- ts(lst[[nEdge]][1], start = ntsp[1], end = ntsp[1],
                              frequency = N)
        } else {
          lst[[nEdge]] <- window(lst[[nEdge]], start = ntsp[1],
                                 end = ntsp[2] - 1 / N)
        }
      }
    }
  }
  return(lst)
}

