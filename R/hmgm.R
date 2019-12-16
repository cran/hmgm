
hmgm<- function(z,y,tune1,tune2,method,kappa,penalty1 = NULL, penalty2 = NULL) {


  size <- function(x, i) {
    return(dim(x)[i])
  }

  p = dim(y)[2]
  q = dim(z)[2]
  q1 = q*(q-1)/2
  L = length(tune1)

  zeros <- function(...) {
    dots <- list(...)
    dims <- as.integer( unlist(dots))
    return(array(0, dims))
  }

  fitlist = list('lambda_j'= zeros(q,L),
                 'lambda_jk' = zeros(q1,2,L),
                 'ita_0' = zeros(p,L),
                 'ita_j' = zeros(p,q,2,L),
                 'phi_0' = zeros(p,p,2,L),
                 'phi_j' = zeros(p,p,q,3,L))

  message(sprintf('fitting logistic regression'))

  repmat <- function(A, ...) {

    nargs <- length(dots <- list(...))
    if (nargs == 1)   dims <- as.integer(dots[[1]])
    else { dims <-   as.integer( unlist(dots)) }
    if (length(dims) == 1) {
      dims[2] <- dims[1]
    }
    if (dims[length(dims)] == 1) {
      return( t(kronecker(array(1, rev(dims)), A)) )
    } else {
      return( kronecker(array(1, dims), A))
    }

  }

  reshape <- function(A, ...) {
    nargs <- length(dots <- list(...))
    if (nargs == 1)   dims <- as.integer(dots[[1]])
    else { dims <-   as.integer(unlist(dots)) }
    return(array(as.vector(A), dims))
  }

  design = function(z, y, type, j) {

    n = nrow(z)
    p = ncol(y)
    q = ncol(z)

    p1 = (p-1)*p/2;

    if (type == 'z') {

      x = zeros(n,q-1+p+p1);
      tmp = z;
      tmp = tmp[,-j]
      x[,1:(q-1)] = tmp
      x[ ,q:(q+p-1)] = y
      ind1 = repmat(1:p,p,1);
      ind1 = ind1[lower.tri(ind1)]
      ind2 = t(repmat(1:p,p,1))
      ind2 = ind2[lower.tri(ind2)]
      yy = y[,ind1] * y[ ,ind2]
      x[,(q+p):(q+p+p1-1)] = yy


    } else if (type == 'y') {

      x = zeros(n,p+p*q-1)
      x[,1:q] = z;
      tmp = y[,-j]
      x[,(q+1):(q+p-1)] = tmp

      a1 = reshape(repmat(t(z), p-1,1),n,(p-1)*q)
      a2 = repmat(tmp, 1,q)

      x[,(q+p):ncol(x)] = a1 * a2

    }

    return(x)

  }

  penalty_fac_flex = function(type,p,q,kappa) {

    p1 = (p-1)*p/2;

    if (type == "z") {

      w = numeric(p+q+p1-1);

      w[1:(q-1)] = 1*kappa;

      w[q:(q+p-1)] = 1;

      w[(q+p):length(w)] = 2;

    } else if (type == 'y') {

      w = numeric(p-1+p*q);

      w[1:q] = 1;

      w[(q+1):(q+p-1)] = 1;

      w[(q+p):length(w)] = 2;

    }

    return(w)

  }

  storage = function(fitlist, a0, beta, j, sigma2 = NULL) {

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

    p = dim(fitlist[['ita_0']])[1]

    q = dim(fitlist[['lambda_j']])[1]

    p1 = (p-1)*p/2

    L = dim(beta)[2]

    start = 0

    if (is.null(sigma2)) {

      if (j > 1) {

        ind1 = transind(q,j,1:(j-1))

      }

      if (j < q) {

        ind2 = transind(q,j, (j+1):q)

      }

      fitlist[['lambda_j']][j,] = a0

      if (j > 1) {

        fitlist[['lambda_jk']][ind1,2,] = beta[(start+1):(start+(j-1)), ]

      }

      start = start + j - 1

      if( j < q) {

        fitlist[['lambda_jk']][ind2,1,] = beta[(start+1):(start+(q-j)),]

      }

      start = start + q - j


      fitlist[['ita_j']][ ,j,1,] = beta[(start+1):(start+p), ]

      start = start+p

      tmp = matrix(1,p,p)

      tmp[!lower.tri(tmp)] = 0

      tmp2 = array(dim = c(p,p,L))

      for( i in 1:L) {

        tmp2[, ,i] = tmp
      }

      tmp = tmp2


      tmp[tmp==1] = beta[(start+1):(start+p1), ]


      for (l in 1:L) {

        fitlist[['phi_j']][,,j,1,l] = tmp[,,l]+ t(tmp[,,l]) - diag(diag(tmp[,,l]));

      }

    } else if (!is.null(sigma2) ) {

      fitlist[['ita_0']][j, ]= a0/sigma2

      beta = sweep(beta,2,sigma2, FUN = "/")

      fitlist[['phi_0']][j,j,1, ] = 1/sigma2

      fitlist[['phi_0']][j,j,2, ] = 1/sigma2

      fitlist[['ita_j']][j, ,2, ] = beta[(start+1):(start+q), ]

      start = start+q

      if(j > 1) {

        fitlist[['phi_0']][j,1:(j-1),2, ] = beta[(start+1):(start+j-1), ]
        fitlist[['phi_0']][1:(j-1),j,2, ] = beta[(start+1):(start+j-1), ]

      }

      start = start+j-1

      if (j < p) {

        fitlist[['phi_0']][j,(j+1):p, 1, ] = beta[(start+1):(start+p-j), ]
        fitlist[['phi_0']][(j+1):p, j,1, ] = beta[(start+1):(start+p-j), ]

      }

      start = start + p - j

      for (k in 1:q) {

        if( j > 1) {

          fitlist[['phi_j']][j,1:(j-1),k,3, ] = beta[(start+1):(start+j-1), ]
          fitlist[['phi_j']][1:(j-1),j,k,3, ] = beta[(start+1):(start+j-1), ]

        }

        start = start+j-1

        if( j < p) {

          fitlist[['phi_j']][j,(j+1):p,k,2, ] =  beta[(start+1):(start+p-j), ]
          fitlist[['phi_j']][(j+1):p, j,k,2, ] =  beta[(start+1):(start+p-j), ]

        }

        start = start+p-j

      }

    }

    return(fitlist)

  }

  for (j in 1:q) {

    message(sprintf('logistic regression %d',j))
    ytmp = (z[,j]>0)+1;
    xtmp = design(z,y,'z',j)


    if(is.null(penalty1)){
      penaltyfac = penalty_fac_flex('z',p,q,kappa)
    }
    else{
      penaltyfac <- penalty1[,j]
    }

    options.type = "covariance"

    if (size(xtmp,2) > size(xtmp,1)) {

      options.type = 'naive';

    }

    fit = glmnet::glmnet(xtmp, ytmp, family = 'binomial',
                 penalty.factor = penaltyfac,
                 maxit=  3000,
                 type.gaussian = options.type,
                 lambda= tune1,
                 standardize = FALSE)

    fit$beta <- as.matrix(fit$beta)

    if (size(fit$beta,2)<L) {

      beta_tmp = cbind(fit$beta , repmat(fit$beta[,ncol(fit$beta)],1,L-size(fit$beta,2)) );
      a0_tmp = c(fit$a0 , repmat(fit$a0[length(fit$a0)],1,L-length(fit$a0)));

    } else {

      beta_tmp = fit$beta;
      a0_tmp = fit$a0;

    }



    fitlist = storage(fitlist, a0_tmp , beta_tmp, j)

  }

  message(sprintf('fitting linear regression') )


  for (j in 1:p) {

    message(sprintf('least square regression %d',j));

    ytmp = y[,j];

    xtmp = design(z,y,'y',j);

    if(is.null(penalty2)){

      penalty_fac_flex = function(type,p,q,kappa) {

        p1 = (p-1)*p/2;

        if (type == "z") {

          w = numeric(p+q+p1-1);

          w[1:(q-1)] = 1*kappa;

          w[q:(q+p-1)] = 1;

          w[(q+p):length(w)] = 2;

        } else if (type == 'y') {

          w = numeric(p-1+p*q);

          w[1:q] = 1;

          w[(q+1):(q+p-1)] = 1;

          w[(q+p):length(w)] = 2;

        }

        return(w)

      }
    }
    else{
      penaltyfac <- penalty2[,j]
    }


    options.type = "covariance"

    if (size(xtmp,2) > size(xtmp,1)) {
      options.type = 'naive';
    }

    fit = glmnet::glmnet(xtmp, ytmp, 'gaussian',
                 penalty.factor = penaltyfac,
                 maxit=  3000,
                 type.gaussian = options.type,
                 lambda= tune2,
                 standardize = FALSE)

    fit$beta <- as.matrix(fit$beta)

    sigma2_tmp = (1-fit$dev)*fit$nulldev;

    if (size(fit$beta,2)<L) {

      beta_tmp = cbind(fit$beta , repmat(fit$beta[,ncol(fit$beta)],1,L-size(fit$beta,2)) );
      a0_tmp = c(fit$a0 , repmat(fit$a0[length(fit$a0)],L-length(fit$a0),1));
      sigma2 = c(sigma2_tmp, repmat(sigma2_tmp[length(sigma2_tmp)],L-length(sigma2_tmp),1));
    } else {

      beta_tmp = fit$beta;
      a0_tmp = fit$a0;
      sigma2 = sigma2_tmp;
    }

    fitlist = storage(fitlist, a0_tmp , beta_tmp, j, sigma2);

  }

  fitpost = function(fitlist, method) {


    p = nrow(fitlist[['ita_0']]);

    q = dim(fitlist[['lambda_j']])[1];

    L = dim(fitlist[['lambda_j']])[2];


    q1 = q*(q-1)/2;

    fitlist_post = list('lambda_j'= matrix(0,q,L),
                        'lambda_jk' = matrix(0,q1,L),
                        'ita_0' = matrix(0, p, L),
                        'ita_j' = array(0, c(p, q, L)),
                        'phi_0' = array(0, c(p,p,L)),
                        'phi_j' = array(0, c(p,p,q,L)));

    fitlist_post[['lambda_j']] = fitlist[['lambda_j']];
    fitlist_post[['ita_0']] = fitlist[['ita_0']];

    if (method == 'max') {

      res = abs(fitlist[['lambda_jk']])

      a = apply(res, c(1,3), max)
      ind = apply(res, c(1,3), which.max)

      s = dim(fitlist[['lambda_jk']]);

      a1 = repmat(1:s[1],1,L)
      a2 = ind
      a3 = reshape(repmat(1:L,s[1],1),L*s[1],1)


      aa = cbind(as.vector(a1) ,as.vector(a2), as.vector(a3))

      res = NULL
      for(i in 1:nrow(aa)) {

        res[i] = nat::sub2ind(s, aa[i,])
      }

      ind1 =  res

      fitlist_post[['lambda_jk']] = reshape(matrix(fitlist[['lambda_jk']][ind1],ncol=1),q1,L);

      res = abs(fitlist[['ita_j']])

      a = apply(res, c(1,2,4), max)
      ind = apply(res, c(1,2,4), which.max)

      s = dim(fitlist[['ita_j']]);


      a1 =  as.vector(repmat(1:s[1],1,L*s[2]))
      a2 = as.vector(t(repmat(reshape(repmat(1:s[2],s[1],1),s[1]*s[2],1),L,1)))
      a3 = as.vector(ind)
      a4 =   as.vector(reshape(repmat(1:L,s[1]*s[2],1),L*s[1]*s[2],1))


      aa = cbind(a1,a2,a3,a4)

      res = NULL
      for(i in 1:nrow(aa)) {

        res[i] = nat::sub2ind(s, aa[i,])
      }

      ind1 =  res

      fitlist_post[['ita_j']] = reshape(array(fitlist[['ita_j']][ind1]),p,q,L);

      res = abs(fitlist[['phi_0']])

      a = apply(res, c(1,2,4), max)
      ind = apply(res, c(1,2,4), which.max)

      s = dim(fitlist[['phi_0']]);


      a1 =  as.vector(repmat(1:s[1],1,L*s[2]))
      a2 = as.vector(t(repmat(reshape(repmat(1:s[2],s[1],1),s[1]*s[2],1),L,1)))
      a3 = as.vector(ind)
      a4 = as.vector(reshape(repmat(1:L,s[1]*s[2],1),L*s[1]*s[2],1))


      aa = cbind(a1,a2,a3,a4)

      res = NULL
      for(i in 1:nrow(aa)) {

        res[i] = nat::sub2ind(s, aa[i,])
      }

      ind1 =  res

      fitlist_post[['phi_0']] = reshape(array(fitlist[['phi_0']][ind1]) ,c(p,p,L));


      res = abs(fitlist[['phi_j']])

      a = apply(res, c(1,2,3,5), max)
      ind = apply(res, c(1,2,3,5), which.max)

      s = dim(fitlist[['phi_j']]);


      a1 =  as.vector(repmat(1:s[1],1,L*s[2]*s[3]))
      a2 = matrix( t(repmat(reshape(repmat(1:s[2],s[1],1),s[1]*s[2],1),s[3]*L,1)),  ncol = 1);
      a3 = as.vector(t(repmat(reshape(repmat(1:s[3],s[1]*s[2],1),s[1]*s[2]*s[3],1),L,1)));
      a4 = as.vector(ind)
      a5 = as.vector(reshape(repmat(1:L,s[1]*s[2]*s[3],1),s[1]*s[2]*s[3]*L,1));


      aa = cbind(a1,a2,a3,a4,a5)

      res = NULL
      for(i in 1:nrow(aa)) {

        res[i] = nat::sub2ind(s, aa[i,])
      }

      ind =  res

      fitlist_post[['phi_j']] = reshape(array(fitlist[['phi_j']][ind]),p,p,q,L);

    }  else if (method == 'min') {



      res = abs(fitlist[['lambda_jk']])

      a = apply(res, c(1,3), min)
      ind = apply(res, c(1,3), which.min)

      s = dim(fitlist[['lambda_jk']]);

      a1 = repmat(1:s[1],1,L)
      a2 = ind
      a3 = reshape(repmat(1:L,s[1],1),L*s[1],1)


      aa = cbind(as.vector(a1) ,as.vector(a2), as.vector(a3))

      res = NULL
      for(i in 1:nrow(aa)) {

        res[i] = nat::sub2ind(s, aa[i,])
      }

      ind1 =  res

      fitlist_post[['lambda_jk']] = reshape(matrix(fitlist[['lambda_jk']][ind1],ncol=1),q1,L);

      res = abs(fitlist[['ita_j']])

      a = apply(res, c(1,2,4), min)
      ind = apply(res, c(1,2,4), which.min)

      s = dim(fitlist[['ita_j']]);


      a1 =  as.vector(repmat(1:s[1],1,L*s[2]))
      a2 = as.vector(t(repmat(reshape(repmat(1:s[2],s[1],1),s[1]*s[2],1),L,1)))
      a3 = as.vector(ind)
      a4 =   as.vector(reshape(repmat(1:L,s[1]*s[2],1),L*s[1]*s[2],1))


      aa = cbind(a1,a2,a3,a4)

      res = NULL
      for(i in 1:nrow(aa)) {

        res[i] = nat::sub2ind(s, aa[i,])
      }

      ind1 =  res

      fitlist_post[['ita_j']] = reshape(array(fitlist[['ita_j']][ind1]),p,q,L);

      res = abs(fitlist[['phi_0']])

      a = apply(res, c(1,2,4), min)
      ind = apply(res, c(1,2,4),which.min)

      s = dim(fitlist[['phi_0']]);


      a1 =  as.vector(repmat(1:s[1],1,L*s[2]))
      a2 = as.vector(t(repmat(reshape(repmat(1:s[2],s[1],1),s[1]*s[2],1),L,1)))
      a3 = as.vector(ind)
      a4 = as.vector(reshape(repmat(1:L,s[1]*s[2],1),L*s[1]*s[2],1))


      aa = cbind(a1,a2,a3,a4)

      res = NULL
      for(i in 1:nrow(aa)) {

        res[i] = nat::sub2ind(s, aa[i,])
      }

      ind1 =  res

      fitlist_post[['phi_0']] = reshape(array(fitlist[['phi_0']][ind1]) ,c(p,p,L));


      res = abs(fitlist[['phi_j']])

      a = apply(res, c(1,2,3,5), min)
      ind = apply(res, c(1,2,3,5), which.min)

      s = dim(fitlist[['phi_j']]);


      a1 =  as.vector(repmat(1:s[1],1,L*s[2]*s[3]))
      a2 = matrix( t(repmat(reshape(repmat(1:s[2],s[1],1),s[1]*s[2],1),s[3]*L,1)),  ncol = 1);
      a3 = as.vector(t(repmat(reshape(repmat(1:s[3],s[1]*s[2],1),s[1]*s[2]*s[3],1),L,1)));
      a4 = as.vector(ind)
      a5 = as.vector(reshape(repmat(1:L,s[1]*s[2]*s[3],1),s[1]*s[2]*s[3]*L,1));


      aa = cbind(a1,a2,a3,a4,a5)

      res = NULL
      for(i in 1:nrow(aa)) {

        res[i] = nat::sub2ind(s, aa[i,])
      }

      ind =  res

      fitlist_post[['phi_j']] = reshape(array(fitlist[['phi_j']][ind]),p,p,q,L);


    }

    return(fitlist_post)
  }

  fitlist_post = fitpost(fitlist, method);


  return(list(fitlist_post = fitlist_post,
              fitlist = fitlist ))

}
