
pargen <- function(adjmat, p, q, a, b, c) {

  posdef <- function(m, force) {

    p <- dim(m)[1]

    if (! all(m == t(m))) {

      stop('input matrix is not symmetric');

    } else if (sum(sum(abs(m)))==0) {

      mnew <- m;

    } else if (force == 0) {

      a <- min(eigen(m)$values);

      if( a >= 0) {

        mnew <- m;

      } else if(a < 0) {

        mnew <- m - 5/4*a*diag(1,p);

      }

    }  else if (force == 1) {

      a <- min(eigen(m)$values);

      if(a > 0) {

        mnew <- m;

      } else if(a < 0) {

        mnew <- m - 5/4*a*diag(1,p);

      } else {

        mnew <- m + abs(mean(apply(m, 2, mean)))/5 * diag(1,p);

      }

    }

    return(mnew)

  }

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

  if ((p+q) != nrow(adjmat)) {
    stop('dimension dismatch')
  }

  adj = list('zz' , 'zy' ,  'yy' )
  adj[['zz']] = adjmat[1:q,1:q]
  adj[['zy']] = adjmat[1:q, (q+1):(q+p)]
  adj[['yy']] = adjmat[(q+1):(q+p), (q+1):(q+p)]


  parlist = list('lambda_j' ,'lambda_jk','ita_0','ita_j', 'phi_0',  'phi_j')

  parlist[['lambda_j']] = c*sign(matrix(stats::runif(q*1,-1,1),q,1)) * matrix(stats::runif(q*1,0.9,1.1),q,1)

  parlist[['lambda_jk']] = c*sign(matrix(stats::runif(q*(q-1)/2*1,-1,1),q*(q-1)/2,1))*matrix(stats::runif(q*(q-1)/2*1,0.9,1.1),q*(q-1)/2,1)

  parlist[['ita_0']] = a*sign(matrix(stats::runif(p*1,-1,1),p,1)) * matrix(stats::runif(p*1,0.9,1.1),p,1)

  parlist[['ita_j']] =  a*sign(matrix(stats::runif(p*q,-1,1),p,q)) * matrix(stats::runif(p*q,0.9,1.1),p,q);


  tmp = b*matrix(stats::runif(p*p,0.9,1.1),p,p)

  tmp[lower.tri(tmp)] = 0

  diag(tmp) = 0


  parlist[['phi_0']] = tmp + t(tmp);

  parlist[['phi_j']] = array(0, c(p,p,q));

  for (j in 1:q) {


    tmp = b*matrix(stats::runif(p*p,0.9,1.1),p,p)

    tmp[lower.tri(tmp)] = 0

    diag(tmp) = 0

    parlist[['phi_j']][,,j] = tmp + t(tmp);


  }

  for (j in 1:(q-1)) {

    for (k in (j+1):q) {

      if(adj[['zz']][j,k]==0) {

        ind = transind(q,j,k);

        parlist[['lambda_jk']][ind] = 0;

      }

    }

  }


  for (j in  1:q) {

    for (k in 1:p) {

      if (adj[['zy']][j,k]==0) {

        parlist[['ita_j']][k,j] = 0;

        parlist[['phi_j']][, k, j] = 0;

        parlist[['phi_j']][k, , j] = 0;


      }
    }

  }

  for (j in 1:(p-1)) {
    for (k in (j+1):p) {

      if (adj[['yy']][j,k]==0) {

        parlist[['phi_0']][j,k] = 0;
        parlist[['phi_0']][k,j] = 0;
        parlist[['phi_j']][j,k,] = 0;
        parlist[['phi_j']][k,j,] = 0;
      }
    }
  }

  parlist[['phi_0']] = posdef(parlist[['phi_0']],1);

  for (j in 1:q) {

    parlist[['phi_j']][,,j] = posdef(parlist[['phi_j']][,,j],0);
    parlist[['phi_0']] = parlist[['phi_0']] + diag(diag(parlist[['phi_j']][,,j]));
    parlist[['phi_j']][,,j] = parlist[['phi_j']][,,j] - diag(diag(parlist[['phi_j']][,,j]));

  }

  return(parlist)

}


