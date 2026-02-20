#! Generate the parameters for profmatch
.problemparameters_profmatch = function(mom_covs, mom_tols, mom_targets, n_units_aux) {

  if (is.vector(mom_covs)){
    mom_covs = t(as.matrix(mom_covs))
  }

  #! Number of units
  n <- nrow(mom_covs)

  #! Number of dec. vars.
  n_dec_vars <- n

  #! Coeffs. of the obj. fun., cvec
  cvec <- rep.int(1, n)

  #! Constraint matrix, Amat
  row_ind_cur = 0
  #! Mom balance
  rows_mom = NULL
  cols_mom = NULL
  vals_mom = NULL
  n_mom_covs = ncol(mom_covs)
  k = 1
  for (i in seq_len(n_mom_covs)) {
    #! Treated
    rows_mom_plus = rep(row_ind_cur + k, n)
    rows_mom_minus = rep(row_ind_cur + k + 1, n)
    rows_mom = c(rows_mom, rows_mom_plus, rows_mom_minus)
    cols_mom = c(cols_mom, rep(1:n, 2))
    vals_plus = c(mom_covs[, i]-mom_targets[i]-mom_tols[i])
    vals_minus = c(mom_covs[, i]-mom_targets[i]+mom_tols[i])
    vals_mom = c(vals_mom, c(vals_plus, vals_minus))
    k = k + 2
  }
  row_ind_cur = max(rows_mom)
  rows_ind = c(rows_mom)
  cols_ind = c(cols_mom)
  vals = c(vals_mom)
  #! Fixed number of units
  if (!is.null(n_units_aux)) {
    rows_units = rep(row_ind_cur + 1L, n)
    cols_units = 1:n
    vals_units = rep(1, n)
    #! Mom balance
    rows_ind = c(rows_ind, rows_units)
    cols_ind = c(cols_ind, cols_units)
    vals = c(vals, vals_units)
    row_ind_cur = max(rows_units)
  }
  aux = cbind(rows_ind, cols_ind, vals)[order(cols_ind), , drop = FALSE]
  aux = aux[aux[, 3] != 0, , drop = FALSE]
  Amat = simple_triplet_matrix(aux[, 1], aux[, 2], aux[, 3])

  #! Constraint vector, bvec
  bvec = NULL
  #! Mom balance
  bvec_mom = rep(0, length(unique(rows_mom)))
  bvec = c(bvec, bvec_mom)
  #! Fixed number of units
  if (!is.null(n_units_aux)) {
    bvec = c(bvec, n_units_aux)
  }

  #! Sense, sense
  sense = NULL
  #! Mom balance
  sense_covs = rep(c("L", "G"), length(unique(rows_mom))/2)
  sense = c(sense, sense_covs)
  #! Fixed number of units
  if (!is.null(n_units_aux)) {
    sense = c(sense, "E")
  }

  #! Variable types, vtype
  vtype = rep("B", n_dec_vars)

  # Output
  list(n = n,
       n_dec_vars = n_dec_vars,
       cvec = cvec,
       Amat = Amat,
       bvec = bvec,
       sense = sense,
       vtype = vtype)
}
