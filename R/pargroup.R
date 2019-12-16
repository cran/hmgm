
pargroup <- function (A) {

  nr <- nrow(A)
  nc <- ncol(A)
  tol <- 2.220446e-16
  r <- 1
  for (i in 1:nc) {
    pivot <- which.max(abs(A[r:nr, i]))
    pivot <- r + pivot - 1
    m <- abs(A[pivot, i])
    if (m <= tol) {
      A[r:nr, i] <- 0
    }
    else {
      A[c(pivot, r), i:nc] <- A[c(r, pivot), i:nc]
      A[r, i:nc] <- A[r, i:nc]/A[r, i]
      if (r == 1) {
        ridx <- c((r + 1):nr)
      }
      else if (r == nr) {
        ridx <- c(1:(r - 1))
      }
      else {
        ridx <- c(1:(r - 1), (r + 1):nr)
      }
      A[ridx, i:nc] <- A[ridx, i:nc] - A[ridx, i, drop = FALSE] %*%
        A[r, i:nc, drop = FALSE]
      if (r == nr)
        break
      r <- r + 1
    }
  }
  A[abs(A) < tol] <- 0

  pos <- c()
  for(i in 1:nc) {
    if(!all(A[,i] %in% c(0,1)) | sum(A[,i] == 1) > 1) {
      pos <- c(pos, i)
    }
  }
  A2 <- A
  while(length(pos) > 1) {
    A2[,pos[1]] <- 0
    temp <- rep(0, ncol(A2))
    temp[pos[1]] <- 1
    pos2 <- c(1)
    for(j in 2:length(pos)) {
      if(all(A[,pos[j]] == A[,pos[1]])) {
        temp[pos[j]] <- 1
        A2[,pos[j]] <- 0
        pos2 <- c(pos2, j)
      }
    }
    A2 <- rbind(A2, temp)
    pos <- setdiff(pos, pos[pos2])
  }

  if(length(pos) == 1) {
    temp <- rep(0, ncol(A2))
    temp[pos[1]] <- 1
    A2[,pos[1]] <- 0
    A2 <- rbind(A2, temp)
  }

  rownames(A2)<- NULL
  A2
}
