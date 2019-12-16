
datagen <- function(parlist, n) {

  zprob <- function(parlist) {

    pq = dim(parlist[['ita_j']])
    p = pq[1]
    q = pq[2]
    prob = as.matrix(rep(0,2^q));


    cparlist = list('g' = as.matrix(rep(0,2^q)),
                    'h' = matrix(0,p,2^q),
                    'K' =  array(0,c(p,p,2^q)))


    z <- lapply(0:(2^q-1),function(x){binaryLogic::as.binary(x,littleEndian=TRUE)})

    maxz <- max(sapply(z, length))

    z <- lapply(z, function(x) {
      if(length(x) < maxz) return(c(x,rep(0, maxz-length(x))))
      return(x)
    })

    z <- do.call('rbind',z)


    for (j in 1:2^q) {

      ztmp = z[j, ];

      zz = ztmp %*% t(ztmp);

      zztmp =  zz[lower.tri(zz)]

      cparlist[['g']][j] = (ztmp %*%parlist[['lambda_j']] + zztmp %*% parlist[['lambda_jk']])[1,1];

      cparlist[['h']][,j] = parlist[['ita_0']] + parlist[['ita_j']] %*% as.matrix(ztmp);

      z3 = array(t(kronecker(ztmp, matrix(1,p,p))), c(p,p,q))

      cparlist[['K']][ , ,j] = parlist[['phi_0']] + apply(z3 * parlist[['phi_j']],c(1,2), sum);

      prob[j] = Matrix::det(cparlist[['K']][ , ,j])^(-0.5) * exp(cparlist[['g']][j]+ ( t(cparlist[['h']][,j]) %*%  solve(cparlist[['K']][ ,,j]) %*% as.matrix(cparlist[['h']][ ,j]/2) )[1,1]);

    }

    prob = prob/sum(prob)

    return(list(prob = prob,
                cparlist = cparlist))

  }

  randp <- function(P, ...) {

    varargin <- unlist(list(...))

    nargin <- length(varargin) + 1

    if(nargin < 2) {
      stop("incorrect dimensions")
    }


    e = tryCatch({

      X = array(stats::runif(prod(varargin) , 0, 1) , varargin)

    }, error = function(e) e)


    if(any(P < 0)) {

      stop('All probabilities should be 0 or larger.');

    }

    if(is.null(P) || sum(P) == 0) {

      stop('All zero probabilities');

    } else {

      numbers <- c(1:length(P))[P!=0]
      prob <- P[P!=0]/sum(P[P!=0])

      X <- array(sample(numbers, size = prod(varargin), replace = TRUE, prob = prob) , varargin)

    }

    return(X)

  }

  pq = dim(parlist[['ita_j']])
  p = pq[1]
  q = pq[2]


  res = zprob(parlist);

  prob = res$prob
  cparlist = res$cparlist

  a = 0;

  b = 0;

  while (a==0 || b==1) {

    ztmp = randp(prob,n,1);

    z <- lapply((ztmp-1)[,1],function(x){binaryLogic::as.binary(x,littleEndian=TRUE, n = q)})

    z <- matrix(as.integer(unlist(z)), ncol = q, byrow = TRUE)



    message('binary z generation finished')

    a = min(apply(z,2,mean));
    b = max(apply(z,2,mean));

  }
  y = matrix(0,n,p);

  for (i in 1:n) {

    sigma = solve(cparlist[['K']][ , , ztmp[i,1]]);

    mu = solve(cparlist[['K']][ , , ztmp[i,1]] ) %*% cparlist[['h']][ ,ztmp[i,1]];

    y[i, ] = MASS::mvrnorm(n = 1, mu, sigma);

  }

  message('data generation finished')

  return(list(z = z,
              y = y,
              prob = prob,
              cparlist = cparlist))

}
