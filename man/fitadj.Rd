\name{fitadj}
\alias{fitadj}

\title{
Obtain the adjascent matrix by thresholding the adj norm matrix
}

\description{
Function to obtain the adjascent matrix by thresholding the adj norm matrix
}

\usage{
fitadj(adj_norm, thres)
}

\arguments{
  \item{adj_norm}{
      A structure with adj norm matrix zz zy yy
}
   \item{thres}{
     Length of thresholding vector
}
}

\value{
 The function returns a 4-dimentional array to record the adj matrix.
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
  \code{\link{edgenorm}}
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

fit = hmgm(z, y, tune1, tune2, 'max',kappa)

#calculate the group L2 norm for each pair of edges

fitlist_post = fit$fitlist_post
adj_norm = edgenorm(fitlist_post)

adj_lambda = fitadj(adj_norm, 0)
}
