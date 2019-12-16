
fitadj  <- function(adj_norm, thres) {

  size <- function(x, i) {
    return(dim(x)[i])
  }

  zeros <- function(...) {
    dots <- list(...)
    dims <- as.integer( unlist(dots))
    return(array(0, dims))
  }

  K = length(thres);

  a = size(adj_norm$zy);

  q = a[1]
  p = a[2]
  L = a[3]

  adj_thres = zeros(p+q, p+q, L, K);

  for (k in 1:K) {

    adj_thres[1:q, 1:q, ,k] = (adj_norm$zz > thres[k]);

    adj_thres[(q+1):(p+q), (q+1):( p+q), , k] = (adj_norm$yy > thres[k]);

    adj_thres[1:q, (q+1):(p+q), ,k] = (adj_norm$zy >thres[k]);

    for (l in 1:L) {

      adj_thres[(q+1):(p+q), 1:q,l,k] = t(adj_thres[1:q, (q+1):(p+q),l,k]);

    }

  }


  return(adj_thres)

}

