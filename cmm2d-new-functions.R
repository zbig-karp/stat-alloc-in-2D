# The constant mixing model for 1 merit dimension ----

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
       `Dissimilarity index` = sum(abs(dat$Actual - aprime))/(2 * sum(dat$Actual)),
       `Distance (Frobenius)` = norm(dat$Actual - aprime, type = "F")/sum(dat$Actual))
  
}

# The new status allocation function for two merit dimensions ----

stall2d <- function(d, dim1, dim2, status, f) {
  
  # Changing the names of the original variables
  
  d <- d[, c(f, dim1, dim2, status)]
  names(d) <- c("f", "dim1", "dim2", "status")
  
  # A function which returns the meritocratic allocation
  mall <- function(x) {
    m <- matrix(0, nrow = nrow(x), ncol = ncol(x))
    dimnames(m) <- dimnames(x)
    for (i in 1:nrow(m)) {
      for (j in 1:ncol(m)) m[i, j] <- 
          min(
            rowSums(x)[i] - sum(m[i, 1:j]), 
            colSums(x)[j] - sum(m[1:i, j])
          )
    }
    m
  }
  
  # Reordering the rows of d by levels of the merit dimensions
  d <- d |>
    arrange(dim1, dim2)
  
  # Cross-classification of the merit dimensions
  t1 <- xtabs(f ~ dim1 + dim2, data = d)
  
  # Status allocation table collapsed over the categories of dim2
  t2 <- xtabs(f ~ dim1 + status, data = d)
  
  # Lottery allocation based on t2
  t3 <- outer(rowSums(t2), colSums(t2), "*")/sum(t2)
  
  # Meritocratic allocation based on t2
  t4 <- mall(t3)
  
  # Joint origin category
  d <- d |>
    unite(col = "origin", dim1, dim2, sep = "-", remove = FALSE) |>
    mutate(origin = fct_inorder(origin))
  
  # Observed status allocation with joint origin categories
  t5 <- xtabs(f ~ origin + status, data = d)
  
  ll <- vector(mode = "list", length = nrow(t1))
  for (i in 1:length(ll)) ll[[i]] <- outer(t1[i, ], t3[i, ], "*")/sum(t1[i, ])
  
  lm <- map(ll, mall)
  
  ml <- vector(mode = "list", length = nrow(t1))
  for (i in 1:length(ml)) ml[[i]] <- outer(t1[i, ], t4[i, ], "*")/sum(t1[i, ])
  
  mm <- map(ml, mall)
  
  out <- list(mm, ml, lm, ll) |>
    map(~do.call(rbind, .))
  
  for (i in 1:length(out)) rownames(out[[i]]) <- rownames(t5)
  
  c(list(t5), out)
}

# Constant mixing in two dimensions ----

cmm2d <- function(dat) {
  mixcoef <- function(dat, mc) {
    a <- mc[1]
    b <- mc[2]
    pred <- (1 - a) * (1 - b) * dat[[5]] + (1 - a) * b * dat[[4]] + a * (1 - b) * dat[[3]] + 
      a * b * dat[[2]]
    norm(dat[[1]] - pred, type = "F")
  }
  
  # Minimising the distance between the observed and the model-predicted allocations
  out <- optim(par = c(0.1, 0.1),
               fn = mixcoef,
               dat = dat,
               method = "L-BFGS-B",
               lower = c(0, 0), upper = c(1, 1),
               hessian = TRUE)
  
  est_table <- tibble(
    Estimate = out$par,
    `S.E.` = sqrt(1/out$hessian[1, 1]),
    `z test` = Estimate/`S.E.`,
    `p value` = pnorm(`z test`, lower.tail = F)
  )
  
  a <- out$par[1]
  b <- out$par[2]
  pred <- (1 - a) * (1 - b) * dat[[5]] + (1 - a) * b * dat[[4]] + a * (1 - b) * dat[[3]] + 
    a * b * dat[[2]]
  
  list(`Mixing coefficient` = as.data.frame(est_table),
       `Model-predicted allocation` = pred,
       `Dissimilarity index` = sum(abs(dat[[1]] - pred))/(2 * sum(dat[[1]])),
       `Distance (Frobenius)` = norm(dat[[1]] - pred, type = "F")/sum(dat[[1]]))
}

# Constant mixing in two dimensions: an extension ----

cmm2d_ext <- function(dat) {
  mixcoef <- function(dat, mc) {
    a <- mc[1]
    bM <- mc[2]
    bL <- mc[3]
    pred <- (1 - a) * (1 - bL) * dat[[5]] + (1 - a) * bL * dat[[4]] + a * (1 - bM) * dat[[3]] + 
      a * bM * dat[[2]]
    norm(dat[[1]] - pred, type = "F")
  }
  
  # Minimising the distance between the observed and the model-predicted allocations
  out <- optim(par = c(0.1, 0.1, 0.1),
               fn = mixcoef,
               dat = dat,
               method = "L-BFGS-B",
               lower = c(0, 0, 0), upper = c(1, 1, 1),
               hessian = TRUE)
  
  est_table <- tibble(
    Estimate = out$par,
    `S.E.` = sqrt(1/diag(out$hessian)),
    `z test` = Estimate/`S.E.`,
    `p value` = pnorm(`z test`, lower.tail = F)
  )
  
  a <- out$par[1]
  bM <- out$par[2]
  bL <- out$par[3]
  pred <- (1 - a) * (1 - bL) * dat[[5]] + (1 - a) * bL * dat[[4]] + a * (1 - bM) * dat[[3]] + 
    a * bM * dat[[2]]
  
  list(`Mixing coefficient` = as.data.frame(est_table),
       `Model-predicted allocation` = pred,
       `Dissimilarity index` = sum(abs(dat[[1]] - pred))/(2 * sum(dat[[1]])),
       `Distance (Frobenius)` = norm(dat[[1]] - pred, type = "F")/sum(dat[[1]]))
}

## Differential mixing in two dimensions ----

dmm2d <- function(dat) {
  
  # Auxiliary objects ----
  # Joint categories of the merit characteristics
  mcats <- rownames(dat[[1]])
  mcats <- str_split(string = mcats, pattern = "-", simplify = TRUE)
  # Categories of the primary characteristic
  mcat1 <- unique(mcats[, 1])
  # Categories of the secondary characteristic
  mcat2 <- unique(mcats[, 2])
  # The number of categories of the primary characteristic
  k1 <- length(mcat1)
  # The number of categories of the secondary characteristic
  k2 <- length(mcat2)
  
  # The function to be optimised over ----
  mixcoef <- function(dat, mc) {
    a <- mc[1:k1]
    b <- mc[(k1 + 1):(k1 + k2)]
    
    # Adjustment proportions
    target <- rep(0, times = ncol(dat[[1]]))
    target <- (colSums(dat[[1]]) - 
                 colSums(as.double(outer(b, a)) * dat[[2]]) - 
                 colSums(as.double(outer(1 - b, a)) * dat[[3]]) - 
                 colSums(as.double(outer(b, 1 - a)) * dat[[4]]))/
      sum(as.double(outer(1 - b, 1 - a)) * rowSums(dat[[5]]))
    
    # Predicted counts
    pred <- as.double(t(outer(a, b))) * dat[[2]] +
      as.double(t(outer(a, 1 - b))) * dat[[3]] + 
      as.double(t(outer(1 - a, b))) * dat[[4]] + 
      as.double(t(outer(1 - a, 1- b))) * outer(rowSums(dat[[5]]), target)
    
    # The distance from prediction to observation
    ll <- norm(pred - dat[[1]], type = "F")
    fake <- 99999999
    if (is.finite(ll)) {
      fake <- ll
    } else {
      ll <- fake
    }
    ll
  }
  
  out <- optim(par = rep(0.5, k1 + k2),
               fn = mixcoef,
               dat = dat,
               method = "L-BFGS-B",
               lower = rep(0, k1 + k2), 
               upper = rep(1, k1 + k2),
               hessian = TRUE)
  
  # A table with mixing coefficients
  est_table <- tibble(
    term = paste(c(rep("Dim 1", k1), rep("Dim 2", k2)), c(mcat1, mcat2), sep = ": "),
    estimate = out$par,
    std.error = sqrt(diag(1/out$hessian)),
    `z-statistic` = estimate/std.error,
    `p-value` = 2 * pnorm(abs(`z-statistic`), lower.tail = FALSE)
  )
  
  # Adjustment proportions
  a <- out$par[1:k1]
  b <- out$par[(k1 + 1):(k1 + k2)]
  target <- rep(0, times = ncol(dat[[1]]))
  target <- (colSums(dat[[1]]) - 
               colSums(as.double(outer(b, a)) * dat[[2]]) - 
               colSums(as.double(outer(1 - b, a)) * dat[[3]]) - 
               colSums(as.double(outer(b, 1 - a)) * dat[[4]]))/
    sum(as.double(outer(1 - b, 1 - a)) * rowSums(dat[[5]]))
  
  
  # A table with adjustment proportions
  adj_prop <- tibble(
    status = colnames(dat[[1]]),
    adj_prop = target
  )
  
  # Predicted counts
  pred <- as.double(t(outer(a, b))) * dat[[2]] +
    as.double(t(outer(a, 1 - b))) * dat[[3]] + 
    as.double(t(outer(1 - a, b))) * dat[[4]] + 
    as.double(t(outer(1 - a, 1- b))) * outer(rowSums(dat[[5]]), target)
  
  # Goodness of fit
  delta <- sum(abs(pred - dat[[1]]))/(2 * sum(dat[[1]]))
  
  # Distance
  distance <- out$value
  
  list("Mixing coefficient" = est_table, "Adjustment proportions" = target,
       "Dissimilarity index" = delta, "Minimised distance" = distance,
       "Predicted counts" = pred)
}

## Differential mixing in two dimensions, maximum likelihood estimation ----

dmm2d_ml <- function(dat) {
  
  # Auxiliary objects ----
  # Joint categories of the merit characteristics
  mcats <- rownames(dat[[1]])
  mcats <- str_split(string = mcats, pattern = "-", simplify = TRUE)
  # Categories of the primary characteristic
  mcat1 <- unique(mcats[, 1])
  # Categories of the secondary characteristic
  mcat2 <- unique(mcats[, 2])
  # The number of categories of the primary characteristic
  k1 <- length(mcat1)
  # The number of categories of the secondary characteristic
  k2 <- length(mcat2)
  
  # The function to be optimised over ----
  mixcoef <- function(dat, mc) {
    a <- mc[1:k1]
    b <- mc[(k1 + 1):(k1 + k2)]
    
    # Adjustment proportions
    target <- rep(0, times = ncol(dat[[1]]))
    target <- (colSums(dat[[1]]) - 
                 colSums(as.double(outer(b, a)) * dat[[2]]) - 
                 colSums(as.double(outer(1 - b, a)) * dat[[3]]) - 
                 colSums(as.double(outer(b, 1 - a)) * dat[[4]]))/
      sum(as.double(outer(1 - b, 1 - a)) * rowSums(dat[[5]]))
    
    # Predicted counts
    pred <- as.double(t(outer(a, b))) * dat[[2]] +
      as.double(t(outer(a, 1 - b))) * dat[[3]] + 
      as.double(t(outer(1 - a, b))) * dat[[4]] + 
      as.double(t(outer(1 - a, 1- b))) * outer(rowSums(dat[[5]]), target)
    
    # The distance from prediction to observation
    ll <- sum(dat[[1]] * log(pred))
    fake <- 99999999
    if (is.finite(ll)) {
      fake <- ll
    } else {
      ll <- fake
    }
    ll
  }
  
  out <- optim(par = rep(0.5, k1 + k2),
               fn = mixcoef,
               dat = dat,
               method = "L-BFGS-B",
               lower = rep(0, k1 + k2), 
               upper = rep(1, k1 + k2),
               hessian = TRUE)
  
  # A table with mixing coefficients
  est_table <- tibble(
    term = paste(c(rep("Dim 1", k1), rep("Dim 2", k2)), c(mcat1, mcat2), sep = ": "),
    estimate = out$par,
    std.error = sqrt(diag(1/out$hessian)),
    `z-statistic` = estimate/std.error,
    `p-value` = 2 * pnorm(abs(`z-statistic`), lower.tail = FALSE)
  )
  
  # Adjustment proportions
  a <- out$par[1:k1]
  b <- out$par[(k1 + 1):(k1 + k2)]
  target <- rep(0, times = ncol(dat[[1]]))
  target <- (colSums(dat[[1]]) - 
               colSums(as.double(outer(b, a)) * dat[[2]]) - 
               colSums(as.double(outer(1 - b, a)) * dat[[3]]) - 
               colSums(as.double(outer(b, 1 - a)) * dat[[4]]))/
    sum(as.double(outer(1 - b, 1 - a)) * rowSums(dat[[5]]))
  
  
  # A table with adjustment proportions
  adj_prop <- tibble(
    status = colnames(dat[[1]]),
    adj_prop = target
  )
  
  # Predicted counts
  pred <- as.double(t(outer(a, b))) * dat[[2]] +
    as.double(t(outer(a, 1 - b))) * dat[[3]] + 
    as.double(t(outer(1 - a, b))) * dat[[4]] + 
    as.double(t(outer(1 - a, 1- b))) * outer(rowSums(dat[[5]]), target)
  
  # Goodness of fit
  delta <- sum(abs(pred - dat[[1]]))/(2 * sum(dat[[1]]))
  
  # Distance
  distance <- out$value
  
  list("Mixing coefficients" = est_table, "Adjustment proportions" = target,
       "Dissimilarity index" = delta, "Minimised distance" = distance,
       "Predicted counts" = pred)
}
