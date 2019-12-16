
edgenorm = function (fitlistpost) {

  transind <- function(dim, j, k) {

    l = rep(0, length(k))


    if (j > dim || max(k) > dim) {

      stop('Incorrect Input')

    } else {

      for (i in 1:length(k)) {

        if (j < k[i]) {

          l[i] = (j-1)*(2*dim-j)/2+k[i]-j;

        } else if (j > k[i]) {

          l[i] = (k[i]-1)*(2*dim-k[i])/2+j-k[i];

        }

      }

    }
    return(l)
  }

  reshape <- function(A, ...) {
    nargs <- length(dots <- list(...))
    if (nargs == 1)   dims <- as.integer(dots[[1]])
    else { dims <-   as.integer(unlist(dots)) }
    return(array(as.vector(A), dims))
  }

  size <- function(x, i) {
    return(dim(x)[i])
  }

  zeros <- function(...) {
    dots <- list(...)
    dims <- as.integer( unlist(dots))
    return(array(0, dims))
  }

  a = size(fitlistpost$phi_j);
  p = a[1]
  p = a[2]
  q = a[3]
  L = a[4]

  adj_norm = list('zz'= zeros(q,q,L),
                  'zy' = zeros(q,p,L),
                  'yy' = zeros(p,p,L));

  for (j in  1:(q-1)) {
    for (k in  (j+1):q) {
      ind = transind(q,j,k);
      adj_norm$zz[j,k,] = fitlistpost$lambda_jk[ind,]^2;
      adj_norm$zz[k,j,] = adj_norm$zz[j,k,];
    }
  }

  for (j in 1:q) {
    for (k in 1:p) {
      par = zeros(p+1,L);
      par[1,] = fitlistpost$ita_j[k,j,];
      par[2:(p+1),] = reshape(fitlistpost$phi_j[,k,j,],p,L);
      adj_norm$zy[j,k,] = apply(par^2, 2, sum)
    }
  }

  for (j in  1:(p-1)) {
    for (k in  (j+1):p) {
      par = zeros(1+q, L);
      par[1,]  = fitlistpost$phi_0[j,k,];
      par[2:(q+1),] = reshape(fitlistpost$phi_j[j,k,,],q,L);
      adj_norm$yy[j,k,] = apply(par^2, 2, sum);
      adj_norm$yy[k,j,] = adj_norm$yy[j,k,];
    }
  }
  return(adj_norm)
}
