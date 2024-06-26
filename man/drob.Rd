% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drob.R
\name{drob}
\alias{drob}
\title{Compute robust estimates of dose-response model parameters}
\usage{
drob(
  x,
  y,
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
)
}
\arguments{
\item{x}{A vector of doses. It is expected that each dose is repeated enough
times so as to be able to compute a good estimate of scale conditional to
the dose during step 2.}

\item{y}{A vector of responses with the same length than \code{x}.}

\item{model}{The string \code{"fpl"} for the predefined 4PL model or, more
generally, a list describing a model with at least two mandatory elements:
\itemize{
\item \code{fun}: the model function, which takes doses \code{x} and parameters \code{t} as
its two only arguments.
\item \code{init}: an initialization function that is able to produce lower and
upper search bounds for step 1. It takes \code{x} and \code{y} as arguments, as well
as an optional \code{extend} argument that controls the extension of the
search region.
}

A third optional element can also be included in the list:
\itemize{
\item \code{grad}: the gradient of \code{fun}, which has the same signature. When \code{grad}
is present it will be used both for computing a gradient in step 3
(if \code{qn_gr} is \code{TRUE}) and for computing the standard errors of the
returned estimates.
}}

\item{step_1}{The loss function used for computing \code{t0}. It may be a function
that takes a vector of residuals as its only argument or one of the following
predefined strings:
\itemize{
\item \code{"lts"}: \code{lts_q}-upper-trimmed mean of squared residuals.
\item \code{"lms"}: median of squared residuals.
\item \code{"ml1"}: mean of absolute residuals. This implements an M-estimate with
\code{rho(r)} = \verb{|r|}, i.e. L1-regression.
\item \code{"sl1"}: median of absolute residuals. This implements an S-estimate with
\code{rho(r)} = \verb{I(|r| > 1)}.
\item \code{"mbi"}: loss for bisquare M-estimate with cutoff point \code{mbi_k}, scaled by
\code{1 / mad(y)} so as to make it scale-equivariant.
\item \code{"sbi"}: (the default) loss for bisquare S-estimate with cutoff point
\code{sbi_k}. The root finding search interval will be extended according to
\code{ms_extend} (see the documentation for \code{m_scale}).
}

\code{loss(y - model$fun(x, t))} will be minimized with respect to \code{t} using a
differential evolution routine provided by the \code{DEoptimR} package. Lower
and upper bounds for the search come from \code{model$init} as documented for
the \code{model} parameter. Other arguments passed to \code{JDEoptim} may be overridden
by passing a \code{de_args} function. If the optimizer fails to converge to a
solution, the entire \code{drob} function execution is aborted with an error.}

\item{step_2}{The function used to compute a scale estimate for each dose.
It may be a function taking a vector of residuals (relative to \code{t0}, as
computed in step 1) and returning a scale estimate, or one of the following
predefined strings:
\itemize{
\item \code{"sl1"}: median of absolute residuals. This is an M-estimate of scale with
\code{rho(r)} = \verb{I(|r| > 1)}. It is scaled by \code{sl1_k}.
\item \code{"sbi"}: (the default) bisquare M-estimate of scale with cutoff point
\code{sbi_k}.
}

The default values for \code{sl1_k} and \code{sbi_k} make the estimates equal to 1
under the standard normal distribution. The function is called once for
each different dose, passing it a vector of residuals for such dose.}

\item{step_3}{The loss function used to compute \code{t}. It may be the
string \code{"mbi"} (the default) for a bisquare loss defining an M-estimate
with cutoff point \code{mbi_k} or, more generally, a list containing at the least
the following mandatory element:
\itemize{
\item \code{rho}: a rho-function taking a vector of scaled residuals.
}

The list may also contain two optional elements:
\itemize{
\item \code{psi}: the first derivative of \code{rho}.
\item \code{dpsi}: the second derivative of \code{rho}.
}

When \code{psi} and \code{dpsi} are present (as well as \code{model$grad}) standard errors
for the returned estimates will be computed. Moreover, when \code{psi} is present
(as well as \code{model$grad} and also \code{qn_gr} is \code{TRUE}) a gradient function
is composed and passed as the argument for the parameter \code{gr} of \code{optim}.

This step minimizes \code{loss((y - model$fun(x, t)) / s)} with respect to \code{t}.
It uses the standard \code{optim} routine with initial parameter \code{t0} (also
used for parameter scaling, see the \code{parscale} control parameter of \code{optim})
and sequentially attempting the following three methods in order:
\enumerate{
\item Quasi-Newton L-BFGS-B with the same bounds than in step 1.
\item Quasi-Newton BFGS.
\item Conjugate gradient CG.
}

The result of the first successful optimizer is kept. If all optimizers
fail to converge to a solution, the entire \code{drob} function execution is
aborted with an error.

Other arguments passed to \code{optim} may be overridden by passing a \code{qn_args}
function.}

\item{lts_q}{The proportion of data to be removed by the upper-trimmed mean
loss (used when \code{step_1} is \code{"lts"}). By default 0.5 to achieve a high
breakdown point.}

\item{mbi_k}{The cutoff point used for bisquare M-estimates (argument
\code{"mbi"} to \code{step_1} or \code{step_3}). By default 3.44 in order to achieve about
85\% asymptotic efficiency under the normal distribution. Other typical values
and their asymptotic efficiencies are: 3.14 (80\%), 3.88 (90\%), 4.68 (95\%).}

\item{sbi_k}{The cutoff point used for bisquare M-estimates of scale and,
consequently, also for bisquare S-estimates (argument \code{"sbi"} to \code{step_1}
or \code{step_2}). By default 1.548 in order to achieve a breakdown point of
about 0.5.}

\item{sl1_k}{Scaling factor for the median of absolute residuals (
argument \code{"sl1"} to \code{step_2}). By default 1.48 in order to get 1 under
the standard normal distribution.}

\item{de_args}{A list that overrides arguments passed to \code{JDEoptim} in
step 1 as in \code{utils::modifyList(args, de_args)}. By default it is empty.}

\item{qn_args}{A list that overrides arguments passed to \code{optim} in
step 3 as in \code{utils::modifyList(args, qn_args)}. By default it is empty.}

\item{qn_gr}{A flag that indicates if a gradient function is to be built and
passed as the argument to the parameter \code{gr} of the \code{optim} routine in
step 3. Besides \code{qn_gr} = \code{TRUE}, for this to be the case both \code{model$grad}
and \code{loss$psi} must be provided. Otherwise, a finite-difference approximation
will be used as per usual. By default \code{FALSE}.}

\item{ms_extend}{When computing bisquare M-estimates of scale and,
consequently, also when computing bisquare S-estimates, \code{ms_extend} will
be passed to \code{m_scale} in order to extend the root finding interval (for
further details, refer to the documentation of \code{m_scale}). By default 5.}

\item{bounds}{A function that takes the \code{init} list produced by the model
and may update \code{init$lower} and/or \code{init$upper} based on domain-specific
considerations. By default it returns \code{init} unchanged.}
}
\value{
A list containing the following elements:
\itemize{
\item \code{t}: a vector with the final location estimate, produced by step 3.
\item \code{t0}: a vector with the initial location estimate, produced by step 1.
\item \code{s}: a vector with a scale estimate for each dose, produced by step 2.
\code{names(s)} will then give the corresponding doses.
\item \code{init}: a list with the results of the initialization phase (\code{model$init}).
It usually includes a non-robust location estimate and its standard error
estimates, besides lower and upper search bounds.
\item \code{loss}: the value of the loss function for step 3 evaluated at \code{t}.
}

If \code{model$grad}, \code{loss$psi} and \code{loss$dpsi} are all provided, the resulting
list will also include a further element:
\itemize{
\item \code{se}: asymptotic standard error estimates for \code{t}.
}
}
\description{
\code{drob} computes an M-estimate of location from data \code{x} and \code{y}, given a
potentially nonlinear and heteroscedastic model. It implements a 3-step
procedure in the spirit of MM-estimation:
\itemize{
\item Step 1: heuristically computes an initial location estimate \code{t0}.
\item Step 2: computes scale estimates \code{s} for each dose, based on the residuals
with respect to \code{t0}.
\item Step 3: starting from \code{t0} and scaling by \code{s} a final location M-estimate
\code{t} is found iteratively.
}

Each step may be fine-tuned using the parameters described below.
}
\examples{

set.seed(0)
t <- c(10, 10, 40, 20)
x <- rep(seq(1, 100, 10), each = 30)
n <- length(x)
cont <- runif(n) > 0.9
e <- rnorm(n, ifelse(cont, 5, 0), (x / 100) * ifelse(cont, 3, 1))
y <- fpl$fun(x, t) + e
est <- drob(x, y)
plot(x, y)
lines(x, fpl$fun(x, t), lwd = 2)
lines(x, fpl$fun(x, est$init$t), col = 2)
lines(x, fpl$fun(x, est$t), col = 3)
legend("topleft", legend = c("t", "t_ls", "t_m"), lty = rep(1, 3), col = 1:3)

}
