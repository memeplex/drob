#' The 4-parameter logistic model (4PL)
#'
#' `fpl` is a list containing three elements related to the 4-parameter logistic
#' model, as required by the `model` parameter of the `drob` function:
#' - `fun`: the 4PL-function itself. It takes as arguments a vector of values
#'   `x` and a vector of 4 parameters `t` .
#' - `grad`: the gradient of the 4PL-function, also expressed as a function
#'   of `x` and `t`.
#' - `init`: a function that computes initial parameters and search bounds.
#'
#' @examples
#'
#' t <- c(10, 10, 40, 20)
#' x <- 1:100
#' eta <- fpl$fun(x, t)
#' y <- eta + rnorm(100)
#' est <- fpl$init(x, y)
#' plot(x, y)
#' lines(x, eta)
#' lines(x, fpl$fun(x, est$t), lty = "dashed")
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
  init = function(x, y, extend = 15, eps1 = 1e-6, eps2 = 1e-12) {
    ux <- unique(x)
    uy <- tapply(y, x, mean)
    r <- range(uy)
    p <- (uy - r[2]) / (r[1] - r[2])
    p <- (1 - 2 * eps1) * p[ux != 0] + eps1
    b <- coef(lm(log(p / (1 - p)) ~ log(ux[ux != 0])))
    t0 <- list(t1 = r[1], t2 = -b[2], t3 = exp(-b[1] / b[2]), t4 = r[2])
    w <- 1 / as.vector(tapply(y, x, var)[as.factor(x)])
    b <- summary(nls(
      y ~ t4 + (t1 - t4) / (1 + (x / t3)^t2), start = t0, weights = w
    ))$coefficients
    i <- if (b[2, 1] < 0) c(4, 2, 3, 1) else 1:4
    t <- unname(b[i, 1])
    t[2] <- abs(t[2])
    se <- unname(b[i, 2])
    lower <- t - extend * se
    upper <- t + extend * se
    lower[2] <- 0
    lower[3] <- eps2 * min(ux[ux != 0])
    list(t = t, se = se, lower = lower, upper = upper)
  }
)

#' Return bisquare and its derivatives
#'
#' This computes bisquare (aka Tukey's biweight) function for a given cutoff
#' point. It also computes its first two derivatives. All three functions
#' are returned as elements of a list with names `rho`, `psi` and `dpsi`,
#' as required by the `step_1` parameter of the `drob` function.
#'
#' @param k The cutoff point below and above which the bisquare function
#'   evaluates to 1.
#'
#' @return A list with three elements:
#' - `rho`: The bisquare function.
#' - `psi`: The first derivative of the bisquare function.
#' - `dpsi`: The second derivative of the bisquare function.
#'
#' @examples
#'
#' x <- seq(-3, 3, length.out = 100)
#' bi2 <- bisquare(2)
#' plot(x, bi2$rho(x), type = "l", ylim = c(-1.5, 1.5))
#' lines(x, bi2$psi(x), col = 2)
#' lines(x, bi2$dpsi(x), col = 3)
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
#' `m_scale` computes an M-estimate of scale for a given rho-function using
#' a one-dimensional root finding routine.
#'
#' @param r A vector of values (typically residuals) whose scale is to be
#'   computed.
#' @param rho The rho-function that defines the M-estimate of scale.
#' @param extend The interval for root finding will extend from 0 to `extend`
#'   times the median absolute value of `r`. Defaults to 5.
#'
#' @return An M-estimate of scale for `r` based on function `rho`.
#'
#' @examples
#' r <- rnorm(1000, 0, 3)
#' rho <- bisquare(1.55)$rho
#' m_scale(r, rho)
#'
#' @export
m_scale <- function(r, rho, extend = 5) {
  f <- function(s) mean(rho(r / s)) - 0.5
  m <- median(abs(r))
  uniroot(f, lower = 0, upper = extend * m, check.conv = TRUE)$root
}

#' Compute robust estimates of dose-response model parameters
#'
#' @description
#'
#' `drob` computes an M-estimate of location from data `x` and `y`, given a
#' potentially nonlinear and heteroscedastic model. It implements a 3-step
#' procedure in the spirit of MM-estimation:
#' - Step 1: heuristically computes an initial location estimate `t0`.
#' - Step 2: computes scale estimates `s` for each dose, based on the residuals
#'   with respect to `t0`.
#' - Step 3: starting from `t0` and scaling by `s` a final location M-estimate
#'   `t` is found iteratively.
#'
#' Each step may be fine-tuned using the parameters described below.
#'
#' @param x A vector of doses. It is expected that each dose is repeated enough
#' times so as to be able to compute a good estimate of scale conditional to
#' the dose during step 2.
#'
#' @param y A vector of responses with the same length than `x`.
#'
#' @param model The string `"fpl"` for the predefined 4PL model or, more
#' generally, a list describing a model with at least two mandatory elements:
#' - `fun`: the model function, which takes doses `x` and parameters `t` as
#'   its two only arguments.
#' - `init`: an initialization function that is able to produce lower and
#'   upper search bounds for step 1. It takes `x` and `y` as arguments, as well
#'   as an optional `extend` argument that controls the extension of the
#'   search region.
#'
#' A third optional element can also be included in the list:
#' - `grad`: the gradient of `fun`, which has the same signature. When `grad`
#'   is present it will be used both for computing a gradient in step 3
#'   (if `qn_gr` is `TRUE`) and for computing the standard errors of the
#'   returned estimates.
#'
#' @param step_1 The loss function used for computing `t0`. It may be a function
#' that takes a vector of residuals as its only argument or one of the following
#' predefined strings:
#' - `"lts"`: `lts_q`-upper-trimmed mean of squared residuals.
#' - `"lms"`: median of squared residuals.
#' - `"ml1"`: mean of absolute residuals. This implements an M-estimate with
#'   `rho(r)` = `|r|`, i.e. L1-regression.
#' - `"sl1"`: median of absolute residuals. This implements an S-estimate with
#'   `rho(r)` = `I(|r| > 1)`.
#' - `"mbi"`: loss for bisquare M-estimate with cutoff point `mbi_k`, scaled by
#'   `1 / mad(y)` so as to make it scale-equivariant.
#' - `"sbi"`: (the default) loss for bisquare S-estimate with cutoff point
#'   `sbi_k`. The root finding search interval will be extended according to
#'   `ms_extend` (see the documentation for `m_scale`).
#'
#' `loss(y - model$fun(x, t))` will be minimized with respect to `t` using a
#' differential evolution routine provided by the `DEoptimR` package. Lower
#' and upper bounds for the search come from `model$init` as documented for
#' the `model` parameter. Other arguments passed to `JDEoptim` may be overridden
#' by passing a `de_args` function. If the optimizer fails to converge to a
#' solution, the entire `drob` function execution is aborted with an error.
#'
#' @param step_2 The function used to compute a scale estimate for each dose.
#' It may be a function taking a vector of residuals (relative to `t0`, as
#' computed in step 1) and returning a scale estimate, or one of the following
#' predefined strings:
#' - `"sl1"`: median of absolute residuals. This is an M-estimate of scale with
#'   `rho(r)` = `I(|r| > 1)`. It is scaled by `sl1_k`.
#' - `"sbi"`: (the default) bisquare M-estimate of scale with cutoff point
#'   `sbi_k`.
#'
#' The default values for `sl1_k` and `sbi_k` make the estimates equal to 1
#' under the standard normal distribution. The function is called once for
#' each different dose, passing it a vector of residuals for such dose.
#'
#' @param step_3 The loss function used to compute `t`. It may be the
#' string `"mbi"` (the default) for a bisquare loss defining an M-estimate
#' with cutoff point `mbi_k` or, more generally, a list containing at the least
#' the following mandatory element:
#' - `rho`: a rho-function taking a vector of scaled residuals.
#'
#' The list may also contain two optional elements:
#' - `psi`: the first derivative of `rho`.
#' - `dpsi`: the second derivative of `rho`.
#'
#' When `psi` and `dpsi` are present (as well as `model$grad`) standard errors
#' for the returned estimates will be computed. Moreover, when `psi` is present
#' (as well as `model$grad` and also `qn_gr` is `TRUE`) a gradient function
#' is composed and passed as the argument for the parameter `gr` of `optim`.
#'
#' This step minimizes `loss((y - model$fun(x, t)) / s)` with respect to `t`.
#' It uses the standard `optim` routine with initial parameter `t0` (also
#' used for parameter scaling, see the `parscale` control parameter of `optim`)
#' and sequentially attempting the following three methods in order:
#' 1. Quasi-Newton L-BFGS-B with the same bounds than in step 1.
#' 2. Quasi-Newton BFGS.
#' 3. Conjugate gradient CG.
#'
#' The result of the first successful optimizer is kept. If all optimizers
#' fail to converge to a solution, the entire `drob` function execution is
#' aborted with an error.
#'
#' Other arguments passed to `optim` may be overridden by passing a `qn_args`
#' function.
#'
#' @param lts_q The proportion of data to be removed by the upper-trimmed mean
#' loss (used when `step_1` is `"lts"`). By default 0.5 to achieve a high
#' breakdown point.
#'
#' @param mbi_k The cutoff point used for bisquare M-estimates (argument
#' `"mbi"` to `step_1` or `step_3`). By default 3.44 in order to achieve about
#' 85% asymptotic efficiency under the normal distribution. Other typical values
#' and their asymptotic efficiencies are: 3.14 (80%), 3.88 (90%), 4.68 (95%).
#'
#' @param sbi_k The cutoff point used for bisquare M-estimates of scale and,
#' consequently, also for bisquare S-estimates (argument `"sbi"` to `step_1`
#' or `step_2`). By default 1.548 in order to achieve a breakdown point of
#' about 0.5.
#'
#' @param sl1_k Scaling factor for the median of absolute residuals (
#' argument `"sl1"` to `step_2`). By default 1.48 in order to get 1 under
#' the standard normal distribution.
#'
#' @param de_args A list that overrides arguments passed to `JDEoptim` in
#' step 1 as in `utils::modifyList(args, de_args)`. By default it is empty.
#'
#' @param qn_args A list that overrides arguments passed to `optim` in
#' step 3 as in `utils::modifyList(args, qn_args)`. By default it is empty.
#'
#' @param qn_gr A flag that indicates if a gradient function is to be built and
#' passed as the argument to the parameter `gr` of the `optim` routine in
#' step 3. Besides `qn_gr` = `TRUE`, for this to be the case both `model$grad`
#' and `loss$psi` must be provided. Otherwise, a finite-difference approximation
#' will be used as per usual. By default `FALSE`.
#'
#' @param ms_extend When computing bisquare M-estimates of scale and,
#' consequently, also when computing bisquare S-estimates, `ms_extend` will
#' be passed to `m_scale` in order to extend the root finding interval (for
#' further details, refer to the documentation of `m_scale`). By default 5.
#'
#' @param bounds A function that takes the `init` list produced by the model
#' and may update `init$lower` and/or `init$upper` based on domain-specific
#' considerations. By default it returns `init` unchanged.
#'
#' @return A list containing the following elements:
#' - `t`: a vector with the final location estimate, produced by step 3.
#' - `t0`: a vector with the initial location estimate, produced by step 1.
#' - `s`: a vector with a scale estimate for each dose, produced by step 2.
#'   `names(s)` will then give the corresponding doses.
#' - `init`: a list with the results of the initialization phase (`model$init`).
#'   It usually includes a non-robust location estimate and its standard error
#'   estimates, besides lower and upper search bounds.
#' - `loss`: the value of the loss function for step 3 evaluated at `t`.
#'
#' If `model$grad`, `loss$psi` and `loss$dpsi` are all provided, the resulting
#' list will also include a further element:
#' - `se`: asymptotic standard error estimates for `t`.
#'
#' @examples
#'
#' set.seed(0)
#' t <- c(10, 10, 40, 20)
#' x <- rep(seq(1, 100, 10), each = 30)
#' n <- length(x)
#' cont <- runif(n) > 0.9
#' e <- rnorm(n, ifelse(cont, 5, 0), (x / 100) * ifelse(cont, 3, 1))
#' y <- fpl$fun(x, t) + e
#' est <- drob(x, y)
#' plot(x, y)
#' lines(x, fpl$fun(x, t), lwd = 2)
#' lines(x, fpl$fun(x, est$init$t), col = 2)
#' lines(x, fpl$fun(x, est$t), col = 3)
#' legend("topleft", legend = c("t", "t_ls", "t_m"), lty = rep(1, 3), col = 1:3)
#'
#' @export
drob <- function( # nolint
  x, y,
  model = "fpl",
  step_1 = "sbi",
  step_2 = "sbi",
  step_3 = "mbi",
  lts_q = 0.5,
  mbi_k = 3.44,
  sbi_k = 1.548,
  sl1_k = 1.48,
  de_args = list(),
  qn_args = list(),
  qn_gr = FALSE,
  ms_extend = 5,
  bounds = identity
) {
  select <- function(arg, ...) if (is.character(arg)) switch(arg, ...) else arg
  mbi <- bisquare(mbi_k)
  sbi <- bisquare(sbi_k)
  model <- select(
    model,
    "fpl" = fpl,
    stop("Invalid model '", model, "'")
  )
  init <- bounds(model$init(x, y))
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
  args <- list(
    fn = function(t) loss(y - model$fun(x, t)), lower = lower, upper = upper,
    fnscale = mloss, tol = 1e-8, maxiter = 500 * length(lower)
  )
  res <- do.call(DEoptimR::JDEoptim, utils::modifyList(args, de_args))
  if (res$convergence == 1) stop("Step 1: optimizer failed")
  t0 <- res$par

  # Step 2 (s, sx)

  scale <- select(
    step_2,
    "sl1" = function(r) median(abs(r)) * sl1_k,
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
    res <- do.call(optim, utils::modifyList(append(args, list(...)), qn_args))
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
    sqrt(((mp / md) * diag(solve(mg))) / length(x))
  }

  list(t = t, t0 = t0, s = s, se = se, init = init, loss = res$value)
}
