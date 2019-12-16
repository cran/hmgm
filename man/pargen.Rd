\name{pargen}

\alias{pargen}

\title{
Generating parameters according to the graph
}

\description{
The function generates parameters for different types of edges based on the graph.
}

\usage{
pargen(adjmat, p, q, a, b, c)
}

\arguments{
  \item{adjmat}{
           A  m x m adjacency matrix (m is the number of total variables). The program automatically check whether the matrix is symmetric and  positive.
}
   \item{p}{
            The number of continous variables.
}
    \item{q}{
              The number of binary variables.
}
    \item{a}{
 Control overall magnitude of the non-zero parameters for edges connecting continuous variables.
}
     \item{b}{
 Control overall magnitude of the non-zero parameters for edges connecting binary and continuous variables.
}
      \item{c}{
 Control overall magnitude of the non-zero parameters for edges connecting binary variables.
}
}

\details{
    In order to generate simulation data, first generate the parameters. Once the adjacency matrix is given, we set all parameters corresponding to absent edges to 0. For the non-zero parameters, we set \code{lambda}j, \code{lambda}jk, \code{eta}j to be positive or negative with equal probability and the absolute value of each
non-zero \code{eta}j is drawn from the uniform distribution on the interval (0.9a, 1.1a) and each non-zero \code{lambda}j or \code{lambda}jk is from (0.9c,1.1c).
The program makes sure that all the probability values are not negative.
}

\value{
The function returns a paramter list.
}

\references{
Jie Cheng, Tianxi Li, Elizaveta Levina, and Ji Zhu. (2017)
    \emph{ High-dimensional Mixed Graphical Models. Journal of Computational and Graphical Statistics 26.2: 367-378},
    \url{https://arxiv.org/pdf/1304.2810.pdf}
}


\author{
    Mingyu Qi, Tianxi Li
}

\seealso{
 \code{\link{datagen}}
}
\examples{
## set controlling parameters
p = 20
q = 10
a = 1
b = 2
c = 1

# set adjacency matrix
adj = matrix(0, p+q, p+q)
adj[10:16, 10:16] = 1
adj[1:5, 1:5] = 1
adj[25:30, 25:30] = 1
adj = adj-diag(diag(adj))

#generate list
parlist = pargen(adj, p, q, a, b,c)

}
