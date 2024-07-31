# Single-dimensional model of status allocation: a function returning a list with three status allocation matrices: actual (or observed), meritocratic, and lottery-based

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

# Single-dimensional model of status-allocation: a function fitting the constant mixing model to the data and estimating the constant mixing coefficient. The function returns a list with four elements: (a) a tibble with the estimated coefficient and relevant statistical tests; (b) a matrix representing the model-predicted status allocation; (c) the index of dissimilarity between the observed and model-predicted status allocations; (d) the minimised distance (i.e., the Frobenius norm of the difference) between the observed and model-predicted status allocations ----

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

# Two-dimensional analysis: a function returning a list with five status allocation matrices -- the observed one and the four reference allocations: M_xM_y, M_xL_y, L_xM_y, and L_xL_y, respectively. The function takes five arguments: ----

# - d, or the data frame containing the data
# - dim1, or the name of the variable in `d` corresponding to the primary merit characteristic
# - dim2, or the name of the variable in `d` corresponding to the secondary merit characteristic
# - status, or the name of the variable in `d` corresponding to the destination status
# - f, or the name of the variable in `d` corresponding to the counts of unique combinations of the values of dim1, dim2, and status

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

# Two-dimensional model of status allocation: a function fitting the model to the data and estimating the mixing coefficients along both merit dimensions. he function returns a list with four elements: (a) a tibble with the estimated coefficients and relevant statistical tests; (b) a matrix representing the model-predicted status allocation; (c) the index of dissimilarity between the observed and model-predicted status allocations; (d) the minimised distance (i.e., the Frobenius norm of the difference) between the observed and model-predicted status allocations----

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
    `S.E.` = sqrt(diag(1/out$hessian)),
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