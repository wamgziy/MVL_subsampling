# =============================================================================
# Simulation script for the multivariate longitudinal LMM in
#   "Optimal Subsampling for Multivariate Longitudinal Data"
# This file generates the FULL DATA used in the simulation study (Section 4).
# Subsampling itself is performed downstream from the .rds files saved here.
#
# Model (Section 2, eqs. (1)-(2)):
#   For subject i (i = 1, ..., m), with n_i time points and L attributes,
#       Y_i = X_i B + Z_i Gamma_i + E_i,
#   or in vectorized form
#       y^v_i = X^v_i beta^v + Z^v_i gamma^v_i + eps^v_i,
#   with X^v_i = I_L kron X_i, Z^v_i = I_L kron Z_i,
#        gamma^v_i ~ N_{qL}(0, D),
#        eps^v_i  ~ N_{n_i L}(0, Sigma kron I_{n_i}).
#
# Simulation settings (Section 4, "Data generating process"):
#   m  = 500,000 subjects, L = 2 attributes
#   beta^v = (5, 10, 15, 20, 25, 30, 35, 40, 45, 50)'  -> length 10
#       (i.e. seq(1:10)*5; with L = 2 this gives p = d = 5 fixed effects per
#        attribute. NOTE: the paper text writes an 11-entry vector starting
#        with 1, which is a typo -- the code's length-10 vector is the truth.)
#   q = 2 random effects (random intercept + random slope), so D is (qL) x (qL) = 4x4
#   Sigma = diag(3, 2)
#   n_i in {4,5,6} w.p. (1/20, 9/10, 1/20)   (Scenario 1, n ~ 5)
#   n_i in {8,9,10} w.p. (1/20, 1/20, 9/10)  (Scenario 2, n ~ 10)
#   X_i rows: multivariate normal (MVN) OR multivariate t with df = 2 (T2)
# =============================================================================

library(MASS)        # mvrnorm()
library(mvtnorm)     # rmvt()
library(doParallel)  # foreach + %dopar% for parallel subject generation


# -----------------------------------------------------------------------------
# Helper: build a d x d covariance/correlation-like matrix from a single
# parameter phi.
#   - ar = 1: AR(1)-style structure with (i,j) entry = phi^|i-j|
#   - ar = 0: compound-symmetric structure with off-diagonals = phi,
#             diagonal = 1 (the off-diagonals are set by the initial
#             matrix(phi, ...) fill; only the diagonal is overwritten)
# Used below to (a) build the within-row covariance for the covariates X_i
# (with ar = 0, scale 1) and (b) seed the random-effects covariance D
# (with ar = 0, then scaled by 4 and individual entries overwritten).
# NOTE: the name "partial_sigma_square" is a bit misleading -- with ar = 0
# this is just a compound-symmetric correlation matrix.
# -----------------------------------------------------------------------------
partial_sigma_square <- function(phi, d, ar = 1){
  A = matrix(phi, nrow = d, ncol = d)
  diag(A) = 1
  if (ar){
    for (i in 1:d){
      for (j in 1:d){
        A[i, j] = phi^(abs(i-j))
      }
    }
  }
  return (A)
}


# -----------------------------------------------------------------------------
# sim.multiv1: generate the full data for a block of subjects under
#   MVN-distributed covariates (Section 4, "Covariate distribution: MVN").
#
# Arguments:
#   beta_v : vectorized fixed-effects parameter beta^v (length d*L)
#   D      : (q*L) x (q*L) random-effects covariance, matching eq. (1)
#   Sigma  : L x L within-subject error covariance; the per-subject error
#            covariance is Sigma kron I_{n_i} (Section 2, parsimonious form)
#   ni     : integer vector of length m giving n_i for each subject
#   L      : number of response attributes
#
# Layout of the returned matrix (one row per (subject, time, attribute), so
# each subject contributes n_i * L rows):
#   columns 1..(p+L)  -> X1: a "long-format" design useful for some downstream
#                        fitting routines (subject-level X stacked with
#                        attribute-level intercept dummies)
#   column  p+L+1     -> Y : stacked response vector y^v_i
#   next   q*L cols   -> Z : block-diagonal Z^v_i = I_L kron Z_i
#   last   p*L cols   -> X : block-diagonal X^v_i = I_L kron X_i (matches eq. (2))
# -----------------------------------------------------------------------------
sim.multiv1 <- function(beta_v, D, Sigma, ni, L){
  p = length(beta_v) / L                  # d in the paper: # fixed-effects per attribute
  m = length(ni)                          # number of subjects in this block
  q = dim(D)[1] / L                       # # random effects per attribute (random intercept + slope -> q = 2)
  index = rep(1:m, times = ni*L)          # row-to-subject map (length sum(n_i)*L)
  Y = rep(NA, sum(ni*L))
  X  = matrix(NA, nrow = sum(ni*L), ncol = p*L)   # block-diagonal X^v_i stacked
  X1 = matrix(NA, nrow = sum(ni*L), ncol = p+L)   # alternative long-format design
  Z  = matrix(NA, nrow = sum(ni*L), ncol = q*L)   # block-diagonal Z^v_i stacked

  for (i in 1:m){
    n_i = ni[i]

    # Generate the (p-1) non-intercept covariates for subject i.
    # Each covariate column is sampled jointly across the n_i time points from
    # a multivariate normal whose covariance is
    #   CS(phi=0.5, dim = p-1) kron I_{n_i},
    # i.e. the (p-1) covariates are CS-correlated across columns and
    # independent across time points.
    b = matrix(
      mvrnorm(1,
              rep(0, n_i*(p-1)),
              kronecker(partial_sigma_square(0.5, (p-1), ar=0), diag(rep(1, n_i)))),
      nrow = n_i, ncol = (p-1), byrow = FALSE
    )
    c = cbind(rep(1, n_i), b)             # X_i: prepend the intercept column 1_{n_i}
                                          # (Section 2: "first column of X_i is 1_{n_i}")

    # Block-diagonal fixed-effects design X^v_i = I_L kron X_i  (eq. (2))
    X[index == i,] = kronecker(diag(rep(1, L)), c)

    # X1 is an alternative long-format layout (NOT used in the paper's eq. (2)):
    #   - first p columns: X_i replicated L times (1_L kron X_i)
    #   - next  L columns: per-attribute intercept dummies (I_L kron 1_{n_i})
    X1[index == i, (1:p)]         = kronecker(rep(1, L), c)
    X1[index == i, (p+1):(p+L)]   = kronecker(diag(rep(1, L)), rep(1, n_i))

    # Random-effects design Z_i = [1_{n_i} | (1, 2, ..., n_i)']  (Section 4),
    # then Z^v_i = I_L kron Z_i.
    Z[index == i,] = kronecker(diag(rep(1, L)), cbind(rep(1, n_i), 1:n_i))

    # Response y^v_i = X^v_i beta^v + Z^v_i gamma^v_i + eps^v_i  (eq. (2))
    #   gamma^v_i ~ N_{qL}(0, D)
    #   eps^v_i: stacked over the L attributes; here Sigma is diagonal so
    #            we sample each attribute's errors independently with
    #            sd = sqrt(Sigma[l,l]).
    # NOTE: this hard-codes L = 2 (the c(rnorm(...), rnorm(...)) call uses
    # only Sigma[1,1] and Sigma[2,2]); generalizing to other L would require
    # a loop over attributes.
    Y[index == i] = X[index == i,] %*% beta_v +
                    Z[index == i,] %*% mvrnorm(1, rep(0, q*L), D) +
                    c(rnorm(n_i, sd = Sigma[1,1]^0.5),
                      rnorm(n_i, sd = Sigma[2,2]^0.5))
  }

  return (cbind(X1, Y, Z, X))
}


# -----------------------------------------------------------------------------
# sim.multiv1.t: same as sim.multiv1 but with HEAVY-TAILED covariates
# (Section 4, "T2" scenario). The (p-1) non-intercept covariates for subject i
# are drawn from a multivariate t distribution with df = 2 and the same
# CS(phi=0.5) kron I_{n_i} scale matrix as the MVN case.
#
# NOTE: although the function takes a `df = 2` argument, df is HARD-CODED to 2
# inside the rmvt() call below (the argument is ignored). Section 4 always
# uses df = 2, so this matches the paper but is fragile if reused.
# -----------------------------------------------------------------------------
sim.multiv1.t <- function(beta_v, D, Sigma, ni, L, df = 2){
  p = length(beta_v) / L
  m = length(ni)
  q = dim(D)[1] / L
  index = rep(1:m, times = ni*L)
  Y = rep(NA, sum(ni*L))
  X  = matrix(NA, nrow = sum(ni*L), ncol = p*L)
  X1 = matrix(NA, nrow = sum(ni*L), ncol = p+L)
  Z  = matrix(NA, nrow = sum(ni*L), ncol = q*L)

  for (i in 1:m){
    n_i = ni[i]

    # Heavy-tailed covariates: rows from multivariate t_2 with the same
    # CS-kron-I scale matrix used in the MVN version. df = 2 hard-coded here.
    b = matrix(
      rmvt(1,
           kronecker(partial_sigma_square(0.5, (p-1), ar=0), diag(rep(1, n_i))),
           df = 2),
      nrow = n_i, ncol = (p-1), byrow = FALSE
    )
    c = cbind(rep(1, n_i), b)

    X[index == i,]                = kronecker(diag(rep(1, L)), c)
    X1[index == i, (1:p)]         = kronecker(rep(1, L), c)
    X1[index == i, (p+1):(p+L)]   = kronecker(diag(rep(1, L)), rep(1, n_i))
    Z[index == i,]                = kronecker(diag(rep(1, L)), cbind(rep(1, n_i), 1:n_i))

    # Response is generated identically to the MVN case (only X is heavy-tailed;
    # gamma^v_i and eps^v_i remain Gaussian, as in Section 2).
    Y[index == i] = X[index == i,] %*% beta_v +
                    Z[index == i,] %*% mvrnorm(1, rep(0, q*L), D) +
                    c(rnorm(n_i, sd = Sigma[1,1]^0.5),
                      rnorm(n_i, sd = Sigma[2,2]^0.5))
  }

  return (cbind(X1, Y, Z, X))
}


# =============================================================================
# Build the random-effects covariance D (Section 4)
#
# Paper specifies the 4 x 4 matrix
#   D = [[4,   1,   1.5, 2  ],
#        [1,   5,   2,   2.5],
#        [1.5, 2,   3,   2  ],
#        [2,   2.5, 2,   4  ]].
# We define D explicitly so the structure is self-evident from the code.
# =============================================================================
dd <- matrix(c(4,   1,   1.5, 2,
               1,   5,   2,   2.5,
               1.5, 2,   3,   2,
               2,   2.5, 2,   4),
             nrow = 4, ncol = 4, byrow = TRUE)


# =============================================================================
# Scenario: MVN covariates, n_i ~ 5  (Section 4, "Scenario 1")
# =============================================================================
n.subjects <- 500000
nsizes <- c(4, 5, 6)
nii <- sample(nsizes, n.subjects, prob = c(1/20, 9/10, 1/20), replace = TRUE)
table(nii)

# Split the n_i vector into 1000 chunks (rows) for parallel generation.
# With n.subjects = 500,000, each chunk has 500 subjects.
nii1 <- matrix(nii, nrow = 1000, byrow = TRUE)

registerDoParallel(cores = detectCores() - 1)
# Each worker generates one chunk of subjects and returns its rows; foreach
# rbinds them into the full design + response matrix.
# beta^v = seq(1:10)*5 = (5, 10, 15, ..., 50) is the length-10 fixed-effects
# vector (with L = 2 attributes -> p = d = 5 fixed effects per attribute).
# The paper's text gives an 11-entry vector starting with 1; that is a typo.
data.mul.foreach <- foreach(i = 1:1000, .combine = 'rbind',
                            .packages = c('mvtnorm', 'MASS')) %dopar% {
  sim.multiv1(seq(1:10)*5, dd, matrix(c(3, 0, 0, 2), 2, 2), nii1[i,], L = 2)
}
# Prepend a subject-id column. Each subject contributes n_i * L rows
# (L = 2 attributes), hence times = nii * 2.
data.mul.foreach <- cbind(as.factor(rep(1:n.subjects, times = nii*2)), data.mul.foreach)
# Column key for the first 13 columns:
#   id      : subject id
#   int     : intercept column from X1
#   X1..X4  : the 4 simulated non-intercept covariates (p - 1 = 4)
#   D1, D2  : per-attribute intercept dummies (the I_L kron 1 block of X1)
#   Y       : stacked response y^v_i
#   int1, time1, int2, time2 : the 4 columns of Z^v_i = I_L kron [1 | t]
# The remaining columns (not named) are the block-diagonal X^v_i used in eq. (2).
colnames(data.mul.foreach)[1:13] = c("id","int","X1","X2","X3","X4",
                                     "D1","D2","Y",
                                     "int1","time1","int2","time2")

saveRDS(list(data.mul.foreach, nii), "data_n5.rds")


# =============================================================================
# Scenario: MVN covariates, n_i ~ 10  (Section 4, "Scenario 2")
# Same setup as above, but n_i in {8, 9, 10} with most subjects at 10.
# =============================================================================
nsizes <- c(8, 9, 10)
nii <- sample(nsizes, n.subjects, prob = c(1/20, 1/20, 9/10), replace = TRUE)
table(nii)
nii1 <- matrix(nii, nrow = 1000, byrow = TRUE)

registerDoParallel(cores = detectCores() - 1)
data.mul.foreach <- foreach(i = 1:1000, .combine = 'rbind',
                            .packages = c('mvtnorm', 'MASS')) %dopar% {
  sim.multiv1(seq(1:10)*5, dd, matrix(c(3, 0, 0, 2), 2, 2), nii1[i,], L = 2)
}
data.mul.foreach <- cbind(as.factor(rep(1:n.subjects, times = nii*2)), data.mul.foreach)
colnames(data.mul.foreach)[1:13] = c("id","int","X1","X2","X3","X4",
                                     "D1","D2","Y",
                                     "int1","time1","int2","time2")
saveRDS(list(data.mul.foreach, nii), "data_n10.rds")


# =============================================================================
# Scenario: T2 covariates (heavy-tailed)  (Section 4, "Covariate distribution: T2")
# Rebuild D defensively (the previous block leaves dd in scope; redefining
# here keeps each scenario block self-contained).
# =============================================================================
dd <- matrix(c(4,   1,   1.5, 2,
               1,   5,   2,   2.5,
               1.5, 2,   3,   2,
               2,   2.5, 2,   4),
             nrow = 4, ncol = 4, byrow = TRUE)


# -----------------------------------------------------------------------------
# T2 covariates, n_i ~ 5
# -----------------------------------------------------------------------------
n.subjects <- 500000
nsizes <- c(4, 5, 6)
nii <- sample(nsizes, n.subjects, prob = c(1/20, 9/10, 1/20), replace = TRUE)
table(nii)
# 1000 chunks of 500 subjects each.
nii1 <- matrix(nii, nrow = 1000, byrow = TRUE)

registerDoParallel(cores = detectCores() - 1)
data.mul.foreach.t <- foreach(i = 1:1000, .combine = 'rbind',
                              .packages = c('mvtnorm', 'MASS')) %dopar% {
  sim.multiv1.t(seq(1:10)*5, dd, matrix(c(3, 0, 0, 2), 2, 2), nii1[i,], L = 2, 2)
}
data.mul.foreach.t <- cbind(as.factor(rep(1:n.subjects, times = nii*2)), data.mul.foreach.t)
colnames(data.mul.foreach.t)[1:13] = c("id","int","X1","X2","X3","X4",
                                       "D1","D2","Y",
                                       "int1","time1","int2","time2")

saveRDS(list(data.mul.foreach.t, nii), "data_n5_t.rds")


# -----------------------------------------------------------------------------
# T2 covariates, n_i ~ 10
# -----------------------------------------------------------------------------
nsizes <- c(8, 9, 10)
nii <- sample(nsizes, n.subjects, prob = c(1/20, 1/20, 9/10), replace = TRUE)
table(nii)
nii1 <- matrix(nii, nrow = 1000, byrow = TRUE)

registerDoParallel(cores = detectCores() - 1)
data.mul.foreach.t <- foreach(i = 1:1000, .combine = 'rbind',
                              .packages = c('mvtnorm', 'MASS')) %dopar% {
  sim.multiv1.t(seq(1:10)*5, dd, matrix(c(3, 0, 0, 2), 2, 2), nii1[i,], L = 2, 2)
}
data.mul.foreach.t <- cbind(as.factor(rep(1:n.subjects, times = nii*2)), data.mul.foreach.t)
colnames(data.mul.foreach.t)[1:13] = c("id","int","X1","X2","X3","X4",
                                       "D1","D2","Y",
                                       "int1","time1","int2","time2")
saveRDS(list(data.mul.foreach.t, nii), "data_n10_t.rds")
