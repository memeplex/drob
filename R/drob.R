library(drc)
library(DEoptimR)

#' The 4-parameter logistic function (4PL)
#'
#' `fpl` is a list containing three elements related to the FPL model, as
#' required by the `model` parameter of the `drob` function:
#' - `fun`: the FPL function itself. `fpl$fun` takes as arguments a vector
#'   of values `x` and a vector of 4 parameters `t` .
#' - `grad`: the gradient of the 4PL function, also expressed as a function
#'   of `x` and `t`.
#' - `init`: a function that computes initial parameters and search bounds.
#'
#' @export
fpl <- list(
  fun = function(x, t) t[4] + (t[1] - t[4]) / (1 + (x / t[3])^t[2]),
  grad = function(x, t) {
    a <- (x / t[3])^t[2]
    b <- ((t[1] - t[4]) / (1 + a)^2) * a
    c(
      1 / (1 + a),
      -b * max(log(x / t[3]), -100),
      b * t[2] / t[3],
      a / (1 + a)
    )
  },
  init = function(x, y, scale = 3) {
    coef <- summary(drm(y ~ x, fct = LL.4()))$coefficients
    idx <- if (coef[1] >= 0) c(3, 1, 4, 2) else c(2, 1, 4, 3)
    t <- unname(coef[idx, 1])
    se <- unname(coef[idx, 2])
    t[2] <- abs(t[2])
    extend <- function(x, dir) x + dir * scale * abs(x)
    init <- list(t = t, se = se, lower = extend(t, -1), upper = extend(t, +1))
    init$lower[c(2, 3)] <- 1e-100
    init
  }
)

#' Return bisquare and its derivatives
#'
#' `bisquare` computes bisquare (aka Tukey's biweight) function for a given
#' cutoff point. It also computes its first two derivatives. All three functions
#' are returned as elements of a list, as required by the `step_1` parameter of
#' the `drob` function.
#'
#' @param k The cutoff point below and above which rho evaluates to 1.
#'
#' @return A list with three elements:
#' - `rho`: The bisquare function
#' - `psi`: The first derivative of rho
#' - `dpsi`: The second derivative of rho
#'
#' @export
bisquare <- function(k) {
  f <- function(r, a, b = 0) ifelse(abs(r) <= k, a, b)
  list(
    rho = function(r) f(r, 1 - (1 - (r / k)^2)^3, 1),
    psi = function(r) f(r, 6 * r * (1 - (r / k)^2)^2 / k^2),
    dpsi = function(r) f(r, 6 * (1 - (r / k)^2) * (1 - 5 * (r / k)^2) / k^2)
  )
}

#' Compute M-estimate of scale
#'
#' `m_scale` computes an M-estimate of scale for a given rho function using
#' one-dimensional root finding.
#'
#' @param r A vector of values (typically residuals) which scale is to be
#'   computed.
#' @param rho The rho function that defines the M-estimate of scale.
#' @param extend The interval for root finding will extend from 0 to `extend`
#'   times the median absolute value of `r`. Defaults to 5.
#'
#' @return The `rho`-based M-estimate of scale for `r`.
#'
#' @export
m_scale <- function(r, rho, extend = 5) {
  f <- function(s) mean(rho(r / s)) - 0.5
  m <- median(abs(r))
  uniroot(f, lower = 0, upper = extend * m, check.conv = TRUE)$root
}

#' @export
drob <- function( # nolint
  x, y,
  model = "fpl",
  step_1 = "sbi",
  step_2 = "sbi",
  step_3 = "mbi",
  lts_q = 0.7,
  mbi_k = 3.44,
  sbi_k = 1.55,
  sli_k = 1.48,
  de_args = identity,
  qn_args = identity,
  qn_gr = FALSE,
  ms_extend = 5,
  init_extend = 5
) {
  select <- function(arg, ...) if (is.character(arg)) switch(arg, ...) else arg
  mbi <- bisquare(mbi_k)
  sbi <- bisquare(sbi_k)
  model <- select(
    model,
    "fpl" = fpl,
    stop("Invalid model '", model, "'")
  )
  init <- model$init(x, y, init_extend)
  lower <- init$lower
  upper <- init$upper

  # Step 1 (t0)

  loss <- select(
    step_1,
    "lts" = function(r) mean((r^2)[r <= quantile(r, lts_q)]),
    "lms" = function(r) median(r^2),
    "ml1" = function(r) mean(abs(r)),
    "sl1" = function(r) median(abs(r)),
    "mbi" = { s <- mad(y); function(r) mean(mbi$rho(r / s)) }, # nolint
    "sbi" = function(r) m_scale(r, sbi$rho, ms_extend),
    stop("Step 1: invalid value '", step_1, "'")
  )
  grid <- matrix(runif(1000, lower, upper), nrow = length(lower))
  mloss <- median(abs(apply(grid, 2, loss)))
  args <- de_args(list(
    fn = function(t) loss(y - model$fun(x, t)),
    lower = lower, upper = upper, fnscale = mloss, tol = 1e-6
  ))
  res <- do.call(JDEoptim, args)
  if (res$convergence == 1) stop("Step 1: optimizer failed")
  t0 <- res$par

  # Step 2 (s, sx)

  scale <- select(
    step_2,
    "sl1" = function(r) median(abs(r)) * sli_k,
    "sbi" = function(r) m_scale(r, sbi$rho, ms_extend),
    stop("Step 2: invalid value '", step_2, "'")
  )
  idx <- as.factor(x)
  s <- tapply(y - model$fun(x, t0), idx, scale)
  sx <- s[idx]

  # Step 3 (t)

  loss <- select(
    step_3,
    "mbi" = mbi,
    stop("Step 3: invalid value '", step_3, "'")
  )
  args <- list(
    fn = function(t) mean(loss$rho((y - model$fun(x, t)) / sx)),
    gr = if (qn_gr && !is.null(model$grad) && !is.null(loss$psi)) function(t) {
      f <- function(r, x) r * model$grad(x, t)
      rowMeans(mapply(f, -loss$psi((y - model$fun(x, t)) / sx) / sx, x))
    },
    par = t0,
    control = list(parscale = t0, trace = 0)
  )
  optimize <- function(...) {
    res <- do.call(optim, qn_args(append(args, list(...))))
    if (res$convergence %in% c(51, 52)) {
      stop("Step 3: optimizer failed with '", res$message, "'")
    }
    res
  }
  res <- try(optimize(method = "L-BFGS-B", lower = lower, upper = upper))
  if (inherits(res, "try-error")) res <- try(optimize(method = "BFGS"))
  if (inherits(res, "try-error")) res <- optimize(method = "CG")
  t <- res$par

  # Standard errors (se)

  se <- if (!is.null(model$grad) && !is.null(loss$psi) && !is.null(loss$dpsi)) {
    r <- (y - model$fun(x, t)) / sx
    mp <- mean(loss$psi(r)^2)
    md <- mean(loss$dpsi(r))^2
    mg <- Reduce("+", mapply(function(x, s) {
      g <- model$grad(x, t)
      outer(g, g) / s^2
    }, x, sx, SIMPLIFY = FALSE)) / length(x)
    sqrt(((mp / md) * diag(ginv(mg))) / length(x))
  }

  list(t = t, t0 = t0, s = s, se = se, init = init, loss = res$value)
}


# https://cran.r-project.org/web/packages/robustbase/vignettes/psi_functions.pdf
# The constant k for 95% efficiency of the regression estimator is 4.685 and
# the constant for a breakdown point of 0.5 of the S-estimator is 1.548

# p.35, rho(r) = I(|r|>1), d = 0.5 => median(abs(r))
# p.130
# When computing an initial estimate β0 for the MM-estimate,
# the bisquare ρ works quite well and we recommend its use for this purpose.
