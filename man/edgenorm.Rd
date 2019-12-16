\name{edgenorm}
\alias{edgenorm}

\title{
Calculate the group L2 norm for each pair of edges
}

\description{
Function to calculate the group L2 norm for each pair of edges
}

\usage{
edgenorm(fitlistpost)
}
\arguments{
  \item{fitlistpost}{
     The fitted parameter path
}
}

\value{
 The function returns a list of group L2 norm for each pair of edges

 \item{zz }{Group L2 norm for each pair of edges connecting binary variables}
 \item{zy }{Group L2 norm for each pair of edges connecting binary variables and continous variables}
 \item{yy }{Group L2 norm for each pair of edges connecting continous variables}
}

\references{
Jie Cheng, Tianxi Li, Elizaveta Levina, and Ji Zhu. (2017)
    \emph{High-dimensional Mixed Graphical Models. Journal of Computational and Graphical Statistics 26.2: 367-378},
    \url{https://arxiv.org/pdf/1304.2810.pdf}
}

\author{
Mingyu Qi, Tianxi Li
}

\seealso{
  \code{\link{hmgm}}
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

##calculate the group L2 norm for each pair of edges

fitlist_post = fit$fitlist_post
adj_norm = edgenorm(fitlist_post)

}
