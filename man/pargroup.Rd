\name{pargroup}
\alias{pargroup}
\title{
    Function to partition overlapping groups into non-overlapping groups
}
\description{
    Function to identify all regions where groups intersect,
    make all variables in each overlapping region into a new group.
}
\usage{
    pargroup(A)
}

\arguments{
  \item{A}{ An n x p matrix represents the relationship between variables and groups. (n is the number of groups and p is the number of variables)
}
}
\details{
     In order to partition groups, we propose a method based on Gaussian-Jordan elimination with pivot on \code{A} to get a reduced row echelon form matrix. Then we use the reduced row echelon form matrix to determine groups. This method can obtain an accurate result as well as reduce computational complexity in R.
}
\value{
    A m x p matrix which represents the relationship between variables and groups after partitioning.
}

\references{
Jie Cheng, Tianxi Li, Elizaveta Levina, and Ji Zhu. (2017)
    \emph{ High-dimensional Mixed Graphical Models. Journal of Computational and Graphical Statistics 26.2: 367-378},
    \url{https://arxiv.org/pdf/1304.2810.pdf}
}

\author{Mingyu Qi, Tianxi Li}


\examples{
## Set an overlap group
A<-rbind(c(1,1,1,0,0), c(0,1,1,1,1))


## Use pargroup() to partion this overlap group to non-overlap group

G = pargroup(A)
}
