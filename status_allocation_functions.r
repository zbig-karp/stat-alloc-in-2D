# The three types of status allocation
mar <- function(x) {
  x <- x[rowSums(x) > 0, colSums(x) > 0]
  l <- rowSums(x) %o% colSums(x) / sum(x)
  dimnames(l) <- dimnames(x)
  m <- matrix(0, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
  for (i in 1:nrow(m)) {
    for (j in 1:ncol(m)) m[i,j] <- min(rowSums(x)[i] - sum(m[i,1:j]), colSums(x)[j] - sum(m[1:i,j]))
  }
  list(Actual = x, Meritocratic = m, Lottery = l)
}

# Constant mixing model, the Pythagorean estimation
merit1 <- function(dat) {
  
  # Function which returns the distance between the observed and the model-predicted allocations
  mc <- function(dat, alpha) {
    aprime <- alpha * dat$Meritocratic + (1 - alpha) * dat$Lottery
    norm(x = dat$Actual - aprime, type = "F")^2
    }
  
  # Minimising the distance between the observed and the model-predicted allocations
  out <- optim(par = 0.9, 
               fn = mc, 
               dat = dat,
               method = "L-BFGS-B", lower = 0, upper = 1,
               hessian = TRUE)
  
  est_table <- tibble(
    Estimate = out$par,
    `S.E.` = sqrt(1/out$hessian[1, 1]),
    `z test` = Estimate/`S.E.`,
    `p value` = pnorm(`z test`, lower.tail = F)
  )
  
  est_table <- est_table %>%
    mutate(`p value` = cut(`p value`, breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf), labels = c("< 0.001", "< 0.0$", "< 0.05", "< 0.1", "n.s.")))
  
  aprime <- out$par * dat$Meritocratic + (1 - out$par) * dat$Lottery
  
  list(`Mixing coefficient` = as.data.frame(est_table),
       `Model-predicted allocation` = aprime,
       `Pearson's X^2` = sum((dat$Actual - aprime)^2/aprime),
       `Distance` = sum(abs(prop.table(dat$Actual) - prop.table(aprime)))/2)
  
}

# Differential mixing model, the Pythagorean estimation
merit2 <- function(dat) {
  
  # Function which returns the distance between the observed and the model-predicted allocations
  mc <- function(dat, alpha) {
    # Row and column marginals
    row_mrg <- rowSums(dat$Actual)
    col_mrg <- colSums(dat$Actual)
    
    # "Target" proportions
    target <- vector(mode = "double", length = length(col_mrg))
    for (i in 1:length(target)) {
      target[i] <- (col_mrg[i] - sum(alpha * dat$Meritocratic[, i]))/sum((1 - alpha) * row_mrg)
    }
    
    # "Compound" allocation
    aprime <- dat$Actual
    for (i in 1:nrow(aprime)) {
      for (j in 1:ncol(aprime)) 
        aprime[i, j] <- (alpha[i] * dat$Meritocratic[i, j] + (1 - alpha[i]) * target[j] * row_mrg[i])
    }
    
    fake <- 99999
    ll <- norm(x = dat$Actual - aprime, type = "F")^2
    
    if (is.finite(ll)) {
      fake <- ll
    } else {
      ll <- fake
    }
    ll
  }
  
  # Minimising the the distance between the observed and the model-predicted allocations
  out <- optim(par = rep(0.9, times = nrow(dat$Actual)), 
               fn = mc, 
               dat = dat,
               method = "L-BFGS-B", 
               lower = rep(0, times = nrow(dat$Actual)), 
               upper = rep(1, times = nrow(dat$Actual)),
               hessian = TRUE)
  
  est_table <- tibble(
    Estimate = out$par,
    `S.E.` =  ifelse(diag(out$hessian) <= 0, NA, sqrt(diag(solve(out$hessian)))),
    `z test` = Estimate/`S.E.`,
    `p value` = pnorm(`z test`, lower.tail = F)
  )
  
  est_table <- est_table %>%
    mutate(`p value` = cut(`p value`, breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf), labels = c("< 0.001", "< 0.0$", "< 0.05", "< 0.1", "n.s.")))
  
  # Calculating expected counts
  row_mrg <- rowSums(dat$Actual)
  col_mrg <- colSums(dat$Actual)
  alpha <- est_table$Estimate
  
  target <- vector(mode = "double", length = length(col_mrg))
  for (i in 1:length(target)) {
    target[i] <- (col_mrg[i] - sum(alpha * dat$Meritocratic[, i]))/sum((1 - alpha) * row_mrg)
  }
  
  aprime <- dat$Actual
  for (i in 1:nrow(aprime)) {
    for (j in 1:ncol(aprime)) 
      aprime[i, j] <- alpha[i] * dat$Meritocratic[i, j] + (1 - alpha[i]) * target[j] * row_mrg[i]
  }
  
  list(`Mixing coefficient` = as.data.frame(est_table),
       `Model-predicted allocation` = aprime,
       `Pearson's X^2` = sum((dat$Actual - aprime)^2/aprime),
       `Distance` = sum(abs(prop.table(dat$Actual) - prop.table(aprime)))/2)
  
}

# Constant mixing model, the ML estimation
merit3 <- function(dat) {
  
  # Function which returns the log-likelihood for the compound matrix
  mc <- function(dat, alpha) {
    aprime <- (alpha * dat$Meritocratic + (1 - alpha) * dat$Lottery)/sum(dat$Actual)
    
    fake <- 99999
    ll <- -sum(dat$Actual * log(aprime))
    
    if (is.finite(ll)) {
      fake <- ll
    } else {
      ll <- fake
    }
    ll
  }
  
  # Maximising the likelihood with respect to the mixing coefficient
  out <- optim(par = 0.5, 
               fn = mc, 
               dat = dat,
               method = "L-BFGS-B", lower = 0, upper = 1,
               hessian = TRUE)
  
  est_table <- tibble(
    Estimate = out$par,
    `S.E.` = sqrt(1/out$hessian[1, 1]),
    `z test` = Estimate/`S.E.`,
    `p value` = pnorm(`z test`, lower.tail = F)
  )
  
  est_table <- est_table %>%
    mutate(`p value` = cut(`p value`, breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf), labels = c("< 0.001", "< 0.0$", "< 0.05", "< 0.1", "n.s.")))
  
  aprime <- out$par * dat$Meritocratic + (1 - out$par) * dat$Lottery
  
  list(`Mixing coefficient` = as.data.frame(est_table),
       `Model-predicted allocation` = aprime,
       `Pearson's X^2` = sum((dat$Actual - aprime)^2/aprime),
       `Distance` = sum(abs(prop.table(dat$Actual) - prop.table(aprime)))/2)
  
}

# Differential mixing model, the ML estimation
merit4 <- function(dat) {
  
  # Function which returns the log-likelihood for the compound matrix
  mc <- function(dat, alpha) {
    # Row and column marginals
    row_mrg <- rowSums(dat$Actual)
    col_mrg <- colSums(dat$Actual)
    
    # "Target" proportions
    target <- vector(mode = "double", length = length(col_mrg))
    for (i in 1:length(target)) {
      target[i] <- (col_mrg[i] - sum(alpha * dat$Meritocratic[, i]))/sum((1 - alpha) * row_mrg)
    }
    
    # "Compound" allocation
    aprime <- dat$Actual
    for (i in 1:nrow(aprime)) {
      for (j in 1:ncol(aprime)) 
        aprime[i, j] <- (alpha[i] * dat$Meritocratic[i, j] + (1 - alpha[i]) * target[j] * row_mrg[i])/sum(dat$Actual)
    }
    
    fake <- 99999
    ll <- -sum(dat$Actual * log(aprime))
    
    if (is.finite(ll)) {
      fake <- ll
    } else {
      ll <- fake
    }
    ll
  }
  
  # Maximising the likelihood with respect to the mixing coefficient
  out <- optim(par = rep(0.5, times = nrow(dat$Actual)), 
               fn = mc, 
               dat = dat,
               method = "L-BFGS-B", 
               lower = rep(0, times = nrow(dat$Actual)), 
               upper = rep(1, times = nrow(dat$Actual)),
               hessian = TRUE)
  
  est_table <- tibble(
    Estimate = out$par,
    `S.E.` =  ifelse(diag(out$hessian) <= 0, NA, sqrt(diag(solve(out$hessian)))),
    `z test` = Estimate/`S.E.`,
    `p value` = pnorm(`z test`, lower.tail = F)
  )
  
  est_table <- est_table %>%
    mutate(`p value` = cut(`p value`, breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf), labels = c("< 0.001", "< 0.0$", "< 0.05", "< 0.1", "n.s.")))
  
  # Calculating expected counts
  row_mrg <- rowSums(dat$Actual)
  col_mrg <- colSums(dat$Actual)
  alpha <- est_table$Estimate
  
  target <- vector(mode = "double", length = length(col_mrg))
  for (i in 1:length(target)) {
    target[i] <- (col_mrg[i] - sum(alpha * dat$Meritocratic[, i]))/sum((1 - alpha) * row_mrg)
  }
  
  aprime <- dat$Actual
  for (i in 1:nrow(aprime)) {
    for (j in 1:ncol(aprime)) 
      aprime[i, j] <- alpha[i] * dat$Meritocratic[i, j] + (1 - alpha[i]) * target[j] * row_mrg[i]
  }
  
  list(`Mixing coefficient` = as.data.frame(est_table),
       `Model-predicted allocation` = aprime,
       `Pearson's X^2` = sum((dat$Actual - aprime)^2/aprime),
       `Distance` = sum(abs(prop.table(dat$Actual) - prop.table(aprime)))/2)
  
}