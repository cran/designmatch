# Solves relaxation problem
.relaxation_n <- function(n_tot, coef, dist_mat, subset_weight, solver, round_cplex, trace) {

  n_dec = (n_tot*(n_tot-1))-sum(1:(n_tot-1))
  #! Nonbipartite matching constraints
  rows_nbm = sort(rep(1:n_tot, n_tot-1))
  temp = matrix(0, nrow = n_tot, ncol = n_tot)
  temp[lower.tri(temp)] = 1:n_dec
  temp = temp + t(temp)
  diag(temp) = NA
  cols_nbm = as.vector(t(temp))
  cols_nbm = cols_nbm[!is.na(cols_nbm)]
  vals_nbm = rep(1, (n_tot-1)*n_tot)

  bvec = rep(1,  length(table(rows_nbm)))
  sense = rep("L", length(table(rows_nbm)))

  aux = cbind(rows_nbm, cols_nbm, vals_nbm)[order(cols_nbm), ]
  Amat = simple_triplet_matrix(i = aux[, 1], j = aux[, 2], v = aux[, 3])

  ub = Inf
  vtype = rep("B", n_dec)
  cvec = coef

  #### SOLVER PART #######

  #! HiGHS
  if (solver == "highs") {
    lhs = rep.int(-Inf, length(sense))
    rhs = rep.int(Inf, length(sense))
    lhs[sense == "G"] = bvec[sense == "G"]
    rhs[sense == "L"] = bvec[sense == "L"]
    lhs[sense == "E"] = bvec[sense == "E"]
    rhs[sense == "E"] = bvec[sense == "E"]

    types <- vtype
    types[types=="B"] <- "I"

    message("  Finding the optimal matches...")
    ptm = proc.time()
    out = highs::highs_solve(L = cvec,
                             lower = 0,
                             upper = ub,
                             A = Amat,
                             lhs = lhs,
                             rhs = rhs,
                             types = types,
                             maximum = TRUE)
    time = (proc.time()-ptm)[3]
    if (out$status %in% c(7, 13)) {
      if (out$status == 7){
        message("  Optimal matches found")
      }
      else {
        message("  Time limit reached!")
      }

      sol = out$primal_solution
      obj = sum((t(dist_mat)[lower.tri(dist_mat)]-(subset_weight*rep(1, n_dec))) * (round(out$primal_solution, 1e-10) == 1))
    }
    else {
      outmessage <- {
        if (out$status == 8) "  Error: problem infeasible!"
        else paste0("  Error: HiGHS solver message: ", out$status_message)
      }
      message(outmessage)

      obj = 0
      sol = NULL
    }
  }

  #! Gurobi
  else if (solver == "gurobi") {
    if (!requireNamespace('gurobi', quietly=TRUE)) {
      stop('Required solver not installed')
    }

    model = list()
    model$obj = cvec
    model$A = Amat
    model$sense = rep(NA, length(sense))
    model$sense[sense=="E"] = '='
    model$sense[sense=="L"] = '<='
    model$sense[sense=="G"] = '>='
    model$rhs = bvec
    model$vtypes = vtype
    model$ub = ub
    model$modelsense = "max"

    params = list(OutputFlag = trace)
    ptm = proc.time()
    out = gurobi::gurobi(model, params)
    time = (proc.time()-ptm)[3]

    if (out$status == "INFEASIBLE") {
      message("  Error: problem infeasible!")
      obj = 0
      sol = NULL
    }

    if (out$status == "OPTIMAL") {
      sol = out$x
      obj = sum((t(dist_mat)[lower.tri(dist_mat)]-(subset_weight*rep(1, n_dec))) * out$x)
    }

  }

  #! CPLEX
  else if (solver == "cplex") {
    if (!requireNamespace('Rcplex', quietly = TRUE)) {
      stop('Required solver not installed')
    }

    ptm = proc.time()
    out = Rcplex::Rcplex(cvec, Amat, bvec, ub = ub, sense = sense, vtype = vtype, n = 1,
                         control = list(trace = trace, round = round_cplex), objsense = "max")
    time = (proc.time()-ptm)[3]

    if (out$status==108) {
      message("  Error: time limit exceeded, no integer solution!")
      obj = 0
      sol = NULL
    }
    else if (is.na(out$obj)) {
      message("  Error: problem infeasible!")
      obj = 0
      sol = NULL
    }
    else {
      sol = out$xopt
      obj = sum((t(dist_mat)[lower.tri(dist_mat)]-(subset_weight*rep(1, n_dec))) * out$xopt)
    }
  }

  #! GLPK
  else if (solver == "glpk") {
    if (!requireNamespace('Rglpk', quietly = TRUE)) {
      stop('Required solver not installed')
    }

    #dir = rep(NA, length(sense))
    #dir[sense=="E"] = '=='
    #dir[sense=="L"] = '<='
    #dir[sense=="G"] = '>='
    #
    #bound = list(lower = list(ind=c(1:length(ub)), val=rep(0,length(ub))),
    #             upper = list(ind=c(1:length(ub)), val=ub))
    #
    #ptm = proc.time()
    #out= Rglpk::Rglpk_solve_LP(cvec, Amat, dir, bvec, bounds = bound, types = vtype, max = TRUE)
    #time = (proc.time()-ptm)[3]
    #
    #if (out$status!=0) {
    #  message("  Error: problem infeasible!")
    #  obj = 0
    #  sol = NULL
    #} else {
    #  sol = out$solution
    #  obj = sum((t(dist_mat)[lower.tri(dist_mat)]-(subset_weight*rep(1, n_dec))) * out$solution)
    #}
    ptm = proc.time()
    Amat = matrix(0, nrow = n_tot, ncol = n_tot)
    res = matrix(0, nrow = n_tot, ncol = n_tot)

    Amat[upper.tri(Amat)] = coef

    max_edge = max(Amat)
    while (max_edge > 0) {
      row = which(Amat == max_edge, arr.ind=TRUE)[1,1]
      col = which(Amat == max_edge, arr.ind=TRUE)[1,2]
      res[row,col] = 1
      Amat[row,] = 0
      Amat[,row] = 0
      Amat[col,] = 0
      Amat[,col] = 0
      max_edge = max(Amat)
    }
    sol = res[upper.tri(res)]
    obj = sum((t(dist_mat)[lower.tri(dist_mat)]-(subset_weight*rep(1, n_dec))) * sol)
    time = (proc.time()-ptm)[3]
  }

  #! Symphony
  else if (solver == "symphony") {
    if (!requireNamespace('Rsymphony', quietly = TRUE)) {
      stop('Required solver not installed')
    }

    dir = rep.int(NA_character_, length(sense))
    dir[sense=="E"] = '=='
    dir[sense=="L"] = '<='
    dir[sense=="G"] = '>='

    bound = list(lower = list(ind=c(1:length(ub)), val=rep(0,length(ub))),
                 upper = list(ind=c(1:length(ub)), val=ub))

    ptm = proc.time()
    out= Rsymphony::Rsymphony_solve_LP(cvec, Amat, dir, bvec, bounds = bound, types = vtype, max = TRUE)
    time = (proc.time()-ptm)[3]

    if (out$status != 0) {
      message("  Error: problem infeasible!")
      obj = 0
      sol = NULL
    }
    else {
      sol = out$solution
      obj = sum((t(dist_mat)[lower.tri(dist_mat)]-(subset_weight*rep(1, n_dec))) * out$solution)
    }
  }

  list(sol = sol, obj = obj, time = time)
}
