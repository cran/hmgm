\name{hmgm}

\alias{hmgm}

\title{
    High-dimensional Mixed Graphical Models Estimation
}

\description{
The main function for high-dimensional Mixed Graphical Models estimation.
}

\usage{
hmgm(z,y,tune1,tune2,method,kappa,penalty1=NULL,penalty2=NULL)
}

\arguments{
  \item{z}{
    z is a n x q discrete data matrix  (n is the sample size and q is the number of discrete variables).
}
 \item{y}{
    y is a n x p continous data matrix  (n is the sample size and p is the number of continous variables).
}
 \item{tune1}{
   Tuning parameter vector for logistic regression (rho in the orginal paper).
}
 \item{tune2}{
   Tuning parameter vector for linear regression (chi in the orginal paper).
}
 \item{method}{
  Can only be max or min, which implies the function takes the maximum or minimum of absolute values as the final estimate.
}
 \item{kappa}{
  tuning parameters for lambda.
}
  \item{penalty1}{
     Penalty for logistics regression. The default penalty is weighted lasso penalty.
     See details at formulation (10) in High-dimensional Mixed Graphical Models.
}
   \item{penalty2}{
    Penalty for linear regression. The default penalty is weighted lasso penalty. See details atformulation (11) in High-dimensional Mixed Graphical Models.
}
}

\details{The graph structure is estimated by maximizing the conditional likelihood of one variable given the rest. We focus on the conditional log-likelihood of each variable and fit separate regressions to estimate the parameters, much in the spirit of the neighborhood selection approach proposed by Meinshausen-Buhlmann for the Gaussian graphical model and by Ravikumar for the Ising model. We incorporating a group lasso penalty, approximated by a weighted lasso penalty for computational efficiency.
}

\value{
 The function returns is a structure of fitted parameters path, the notations are the same as the paper.

 \item{fitlist_post}{the fitted parameter path by taking the maximum or minimum absolute values with signs}
 \item{fitlist }{The original fitlist}
}

\references{
Jie Cheng, Tianxi Li, Elizaveta Levina, and Ji Zhu.(2017)
   \emph{ High-dimensional Mixed Graphical Models. Journal of Computational and Graphical Statistics 26.2: 367-378},
   \url{https://arxiv.org/pdf/1304.2810.pdf}\cr
Simon, N., Friedman, J., Hastie,T., Tibshirani, R. (2011)
\emph{Regularization Paths for Cox's ProportionalHazards Model via Coordinate Descent, Journal of Statistical Software, Vol.39(5) 1-13},
\url{https://www.jstatsoft.org/v39/i05/}\cr
Meinshausen, N. and Buhlmann, P. (2006)
  \emph{High dimensional graphs and variable selection with the lasso, Annals of Statistics, 34, 1436–1462.},
  \url{https://arxiv.org/pdf/math/0608017.pdf}\cr
Ravikumar, P., Wainwright, M., and Lafferty, J. (2010)
\emph{ High-dimensionalIsing model selection using l1-regularized logistic regression,Annals of
Statistics, 38, 1287–1319.},
\url{https://arxiv.org/pdf/1010.0311.pdf} \cr
Liu, H., Han, F., Yuan, M., Lafferty, J., and Wasserman, L. (2012)
\emph{ High dimensional semiparametric Gaussian copula graphical models, Annals
of Statistics, 40, 2293–2326.},
\url{https://arxiv.org/pdf/1202.2169.pdf}\cr
Zhao, P., Rocha, G., and Yu, B. (2009)
\emph{ The composite absolute penalties family for grouped and hierarchical variable selection, The Annals of
Statistics, 3468–3497.},
\url{https://arxiv.org/pdf/0909.0411.pdf}
}


\author{
     Mingyu Qi, Tianxi Li
}

\seealso{
\code{\link{datagen}}
}
\examples{
n = 100
p = 20
q = 10
a = 1
b=  2
c = 1


adj = matrix(0, p+q, p+q)
adj[10:16, 10:16] = 1
adj[1:5, 1:5] = 1
adj[25:30, 25:30] = 1
adj = adj-diag(diag(adj))

parlist = pargen(adj, p, q, a, b,c)

mydata = datagen(parlist, n)

z = mydata$z

y = mydata$y

tune1 = tune2 = 0.1

kappa = 0.1

## parameter estimation

fit = hmgm(z, y, tune1, tune2, 'max', kappa)

}
