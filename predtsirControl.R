
predtsirControl <- function (times, births, beta, alpha, S0, I0, nsim, stochastic, controlStart,  controlEnd, betachange) 
{
  I.mat <- S.mat <- matrix(NA, length(times), nsim)
  alpha <- alpha[1]
  if (length(beta) < length(times)) {
    beta <- rep(beta, length(times))[1:length(times)]
  }
  if (length(births) < length(times)) {
    births <- rep(births, length(times))
  }
  for (n in 1:nsim) {
    S <- I <- rep(NA, length(times))
    S[1] <- round(S0)
    I[1] <- round(I0)
    for (t in 2:length(times)) {
      lambda <- min(S[t - 1], unname(beta[t - 1] * S[t - 
                                                       1] * (I[t - 1])^alpha))
      if(t >= controlStart & t <= controlEnd){ lambda <- min(S[t - 1], unname(betachange*beta[t - 1] * S[t - 
                                                                                              1] * (I[t - 1])^alpha))}
      if (stochastic) {
        I[t] <- rnbinom(n = 1, mu = lambda, size = I[t - 
                                                       1] + 1e-10)
      }
      else {
        I[t] <- lambda
      }
      S[t] <- max(S[t - 1] + births[t - 1] - I[t], 0)
    }
    I.mat[, n] <- I
    S.mat[, n] <- S
  }
  I.mat <- data.frame(I.mat)
  S.mat <- data.frame(S.mat)
  I.error <- apply(as.matrix(I.mat), 1, function(x) {
    mean(x) + c(-1.96, 1.96) * sd(x)/sqrt(length(x))
  })
  I.mat$high <- I.error[2, ]
  I.mat$low <- I.error[1, ]
  S.error <- apply(as.matrix(S.mat), 1, function(x) {
    mean(x) + c(-1.96, 1.96) * sd(x)/sqrt(length(x))
  })
  S.mat$high <- S.error[2, ]
  S.mat$low <- S.error[1, ]
  I.mat$mean <- rowMeans(I.mat[, 1:nsim], na.rm = T)
  S.mat$mean <- rowMeans(S.mat[, 1:nsim], na.rm = T)
  I.mat$time <- times
  S.mat$time <- times
  return(list(I = I.mat, S = S.mat))
}
