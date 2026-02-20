.oneprob_profmatch <- function(level, t_ind, mom, n_units_aux, solver) {
  mom_covs <- mom$covs
  mom_tols <- mom$tols
  mom_targets <- mom$targets

  if (is.null(solver)) {
    solver <- list(t_max = 60 * 15,
                   approximate = 1,
                   name = "highs")
  }

  t_max <- solver$t_max
  approximate <- solver$approximate
  trace <- solver$trace
  round_cplex <- solver$round_cplex
  solver <- solver$name

  .mom_covs <- mom_covs[which(t_ind == level),]

  #! Generate the parameters
  message("  Building the matching problem...")
  prmtrs <- .problemparameters_profmatch(.mom_covs, mom_tols, mom_targets, n_units_aux)

  n = prmtrs$n
  n_dec_vars = prmtrs$n_dec_vars
  cvec = prmtrs$cvec
  Amat = prmtrs$Amat
  bvec = prmtrs$bvec
  sense = prmtrs$sense
  vtype = prmtrs$vtype

  #! Find matches and calculate the elapsed time
  #! HiGHS
  if (solver == "highs") {
    message("  HiGHS optimizer is open...")

    lhs = rep.int(-Inf, length(sense))
    rhs = rep.int(Inf, length(sense))
    lhs[sense == "G"] = bvec[sense == "G"]
    rhs[sense == "L"] = bvec[sense == "L"]
    lhs[sense == "E"] = bvec[sense == "E"]
    rhs[sense == "E"] = bvec[sense == "E"]

    types = vtype
    types[types == "B"] = "I"

    message("  Finding the optimal matches...")
    ptm = proc.time()
    out = highs::highs_solve(L = cvec,
                             lower = 0,
                             upper = 1,
                             A = Amat,
                             lhs = lhs,
                             rhs = rhs,
                             types = types,
                             maximum = TRUE,
                             control = highs::highs_control(time_limit = t_max))
    time = (proc.time()-ptm)[3]

    if (out$status %in% c(7, 13)) {
      if (out$status == 7){
        message("  Optimal matches found")
      }
      else if (out$status == 13){
        message("  Time limit reached!")
      }

      #! Matched units indexes
      id = (1:n_dec_vars)[round(out$primal_solution, 1e-10)==1]

      #! Optimal value of the objective function
      obj_total = out$objective_value
    }
    else {
      outmessage <- {
        if (out$status == 8) "  Error: problem infeasible!"
        else paste0("  Error: HiGHS solver message: ", out$status_message)
      }
      message(outmessage)

      obj_total = NA
      id = NA
      time = NA
    }
  }

  #! Gurobi
  else if (solver == "gurobi") {

    if (!requireNamespace('gurobi', quietly = TRUE)) {
      stop('Required solver not installed')
    }

    message(format("  Gurobi optimizer is open..."))

    model = list()
    model$modelsense = 'max'
    model$obj = cvec
    model$A = Amat
    model$sense = rep(NA, length(sense))
    model$sense[sense=="E"] = '='
    model$sense[sense=="L"] = '<='
    model$sense[sense=="G"] = '>='
    model$rhs = bvec
    model$vtype = vtype

    t_lim = list(TimeLimit = t_max, OutputFlag = trace)

    message("  Finding the optimal matches...")
    ptm = proc.time()
    out = gurobi::gurobi(model, t_lim)
    time = (proc.time()-ptm)[3]

    if (out$status == "INFEASIBLE") {
      message("  Error: problem infeasible!")
      obj_total = NA
      obj_dist_mat = NA
      id = NA
      time = NA
    }
    else if (out$status %in%  c("OPTIMAL", "TIME_LIMIT")) {
      if (out$status == "OPTIMAL") {
        message("  Optimal matches found")
      }
      else {
        message("  Time limit reached, best suboptimal solution given")
      }

      #! Matched units indexes
      id = (1:n_dec_vars)[out$x==1]

      #! Optimal value of the objective function
      obj_total = out$objval
    }
  }

  #! CPLEX
  else if (solver == "cplex") {
    if (!requireNamespace('Rcplex', quietly = TRUE)) {
      stop('Required solver not installed')
    }

    message("  CPLEX optimizer is open...")
    message("  Finding the optimal matches...")

    ptm = proc.time()
    out = Rcplex::Rcplex(objsense = 'max', cvec, Amat, bvec, sense = sense, vtype = vtype, n = 1,
                         control = list(trace = trace, round = round_cplex, tilim = t_max))

    time = (proc.time()-ptm)[3]

    if (out$status == 108) {
      message("  Error: time limit exceeded, no integer solution!")
      obj_total = NA
      obj_dist_mat = NA
      id = NA
      time = NA
    }
    else if (is.na(out$obj)) {
      message("  Error: problem infeasible!")
      obj_total = NA
      obj_dist_mat = NA
      id = NA
      time = NA
    }
    else {
      message("  Optimal matches found")

      #! Matched units indexes
      id = (1:n_dec_vars)[out$xopt==1]

      #! Optimal value of the objective function
      obj_total = out$obj
    }
  }

  #! GLPK
  else if (solver == "glpk") {
    if (!requireNamespace('Rglpk', quietly=TRUE)) {
      stop('Required solver not installed')
    }

    message("  GLPK optimizer is open...")
    dir = rep.int(NA_character_, length(prmtrs$sense))
    dir[prmtrs$sense=="E"] = '=='
    dir[prmtrs$sense=="L"] = '<='
    dir[prmtrs$sense=="G"] = '>='

    message("  Finding the optimal matches...")
    ptm = proc.time()
    out <- Rglpk::Rglpk_solve_LP(cvec, Amat, dir, bvec, types = vtype, max = TRUE)
    time = (proc.time()-ptm)[3]

    if (out$status != 0) {
      message("  Error: problem infeasible!")
      obj_total = NA
      obj_dist_mat = NA
      id = NA
      time = NA
    }
    else {
      message("  Optimal matches found")

      #! Matched units indexes
      id = (1:n_dec_vars)[t_ind==1 & out$solution==1]

      #! Optimal value of the objective function
      obj_total = out$optimum
    }
  }

  #! Symphony
  else {
    if (!requireNamespace('Rsymphony', quietly=TRUE)) {
      stop('Required solver not installed')
    }

    message("  Symphony optimizer is open...")

    dir = rep.int(NA_character_, length(prmtrs$sense))
    dir[prmtrs$sense=="E"] = '=='
    dir[prmtrs$sense=="L"] = '<='
    dir[prmtrs$sense=="G"] = '>='

    message("  Finding the optimal matches...")
    ptm = proc.time()
    out <- Rsymphony::Rsymphony_solve_LP(cvec, Amat, dir, bvec, types = vtype, max = TRUE,
                                         time_limit = t_max)
    time = (proc.time()-ptm)[3]

    if (out$status == 0) {
      message("  Optimal matches found")

      #! Matched units indexes
      id = (1:n_dec_vars)[out$solution==1]

      #! Optimal value of the objective function
      obj_total = out$objval
    }
    else {
      if (out$status == 228) {
        message("  Error: problem exceeded the time limit and no feasible solution is found!")
      }
      else {
        message("  Error: problem infeasible!")
      }

      obj_total = NA
      obj_dist_mat = NA
      id = NA
      time = NA
    }
  }

  #! Output
  list(obj_total = obj_total,
       id = id,
       time = time)
}
