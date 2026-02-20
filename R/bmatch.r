#! bmatch
bmatch <- function(t_ind, dist_mat = NULL, subset_weight = NULL, n_controls = 1, total_groups = NULL,
                   mom = NULL,
                   ks = NULL,
                   exact = NULL,
                   near_exact = NULL,
                   fine = NULL,
                   near_fine = NULL,
                   near = NULL,
                   far = NULL,
                   #use_controls = NULL,
                   solver = NULL) {

  use_controls = NULL

  if (is.null(mom)) {
    mom <- list()
  }

  mom_covs <- mom$covs
  mom_tols <- mom$tols
  mom_targets <- mom$targets

  if (is.null(ks)) {
    ks <- list(n_grid = 10)
  }

  ks_covs <- ks$covs
  ks_n_grid <- ks$n_grid
  ks_tols <- ks$tols

  if (is.null(exact)) {
    exact <- list()
  }

  exact_covs <- exact$covs

  if (is.null(near_exact)) {
    near_exact <- list()
  }

  near_exact_covs <- near_exact$covs
  near_exact_devs <- near_exact$devs

  if (is.null(fine)) {
    fine <- list()
  }

  fine_covs <- fine$covs

  if (is.null(near_fine)) {
    near_fine <- list()
  }

  near_fine_covs <- near_fine$covs
  near_fine_devs <- near_fine$devs

  if (is.null(near)) {
    near <- list()
  }

  near_covs <- near$covs
  near_pairs <- near$pairs
  near_groups <- near$groups

  if (is.null(far)) {
    far <- list()
  }

  far_covs <- far$covs
  far_pairs <- far$pairs
  far_groups <- far$groups

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

  #! CALL ERROR HANDLING

  if (is.null(subset_weight)) {
    subset_weight <- 0
  }

  #! Generate the parameters
  message("  Building the matching problem...")
  prmtrs = .problemparameters(t_ind, dist_mat, subset_weight, n_controls, total_groups,
                              mom_covs, mom_tols, mom_targets,
                              ks_covs, ks_n_grid, ks_tols,
                              exact_covs,
                              near_exact_covs, near_exact_devs,
                              fine_covs,
                              near_fine_covs, near_fine_devs,
                              near_covs, near_pairs, near_groups,
                              far_covs, far_pairs, far_groups,
                              use_controls,
                              approximate)
  n_t <- prmtrs$n_t
  n_c <- prmtrs$n_c

  cvec <- prmtrs$cvec
  Amat <- prmtrs$Amat
  bvec <- prmtrs$bvec
  ub <- prmtrs$ub
  sense <- prmtrs$sense
  vtype <- prmtrs$vtype
  c_index <- prmtrs$c_index

  #! Find matches and calculate the elapsed time
  #! HiGHS
  if (solver == "highs") {
    message("  HiGHS optimizer is open...")

    lhs <- rep.int(-Inf, length(sense))
    rhs <- rep.int(Inf, length(sense))
    lhs[sense == "G"] <- bvec[sense == "G"]
    rhs[sense == "L"] <- bvec[sense == "L"]
    lhs[sense == "E"] <- bvec[sense == "E"]
    rhs[sense == "E"] <- bvec[sense == "E"]

    types <- vtype
    types[types == "B"] <- "I"

    message("  Finding the optimal matches...")
    ptm <- proc.time()
    out <- highs::highs_solve(L = cvec,
                             lower = 0,
                             upper = ub,
                             A = Amat,
                             lhs = lhs,
                             rhs = rhs,
                             types = types,
                             control = highs::highs_control(time_limit = t_max))
    time <- (proc.time()-ptm)[3L]

    if (out$status %in% c(7, 13)) {
      if (out$status == 7) {
        message("  Optimal matches found")
      }
      else {
        message("  Time limit reached!")
      }

      if (approximate == 1) {
        rel <- .relaxation_b(n_t, n_c, out$primal_solution, dist_mat,
                             subset_weight, "highs", round_cplex, trace)
        out$primal_solution <- rel$sol
        out$objval <- rel$obj
        time <- time + rel$time
      }

      #! Matched units indexes
      t_id <- unique(sort(rep(1:n_t, n_c))[round(out$primal_solution[1:(n_t * n_c)], 1e-10) == 1])
      c_id <- (c_index + n_t)[round(out$primal_solution[1:(n_t*n_c)], 1e-10)==1]

      #! Group (or pair) identifier
      group_id_t <- 1:(length(t_id))
      group_id_c <- sort(rep(1:(length(t_id)), n_controls))
      group_id <- c(group_id_t, group_id_c)

      #! Optimal value of the objective function
      obj_total <- out$objval

      if (!is.null(dist_mat)) {
        obj_dist_mat <- sum(c(as.vector(matrix(t(dist_mat), nrow = 1, byrow = TRUE)) * round(out$primal_solution[1:(n_t * n_c)], 1e-10) == 1))
      }
      else {
        obj_dist_mat = NULL
      }
    }
    else {
      outmessage <- {
        if (out$status == 8) "  Error: problem infeasible!"
        else paste0("  Error: HiGHS solver message: ", out$status_message)
      }
      message(outmessage)

      obj_total <- NA
      id <- NA
      time <- NA
    }
  }

  #! Gurobi
  else if (solver == "gurobi") {
    if (!requireNamespace('gurobi', quietly = TRUE)) {
      stop('Required solver not installed')
    }

    message("  Gurobi optimizer is open...")
    model <- list(obj = cvec,
                  A = Amat,
                  sense = rep.int(NA_character_, length(sense)),
                  rhs = bvec,
                  vtypes = vtype,
                  ub = ub)

    model$sense[sense == "E"] <- '='
    model$sense[sense == "L"] <- '<='
    model$sense[sense == "G"] <- '>='

    t_lim <- list(TimeLimit = t_max, OutputFlag = trace)

    message("  Finding the optimal matches...")
    ptm <- proc.time()
    out <- gurobi::gurobi(model, t_lim)
    time <- (proc.time()-ptm)[3]

    if (out$status == "INFEASIBLE") {
      message("  Error: problem infeasible!")
      obj_total = NA
      obj_dist_mat = NA
      t_id = NA
      c_id = NA
      group_id = NA
      time = NA
    }
    else if (out$status %in% c("OPTIMAL", "TIME_LIMIT")) {

      if (out$status == "OPTIMAL") {
        message("  Optimal matches found")
      }
      else {
        message("  Time limit reached, best suboptimal solution given")
      }

      if (approximate == 1) {
        rel = .relaxation_b(n_t, n_c, out$x, dist_mat, subset_weight, "gurobi",
                            round_cplex, trace)
        out$x = rel$sol
        out$objval = rel$obj
        time = time + rel$time
      }

      #! Matched units indexes
      t_id = unique(sort(rep(1:n_t, n_c))[out$x[1:(n_t*n_c)]==1])
      c_id = (c_index+n_t)[out$x[1:(n_t*n_c)]==1]

      #! Group (or pair) identifier
      group_id_t = 1:(length(t_id))
      group_id_c = sort(rep(1:(length(t_id)), n_controls))
      group_id = c(group_id_t, group_id_c)

      #! Optimal value of the objective function
      obj_total = out$objval

      if (!is.null(dist_mat)) {
        obj_dist_mat = sum(c(as.vector(matrix(t(dist_mat), nrow = 1, byrow = TRUE)) * out$x[1:(n_t*n_c)]==1))
      }
      else {
        obj_dist_mat = NULL
      }
    }
  }

  #! CPLEX
  else if (solver == "cplex") {
    if (!requireNamespace('Rcplex', quietly=TRUE)) {
      stop('Required solver not installed')
    }

    message("  CPLEX optimizer is open...")
    message("  Finding the optimal matches...")
    ptm = proc.time()
    out = Rcplex::Rcplex(cvec, Amat, bvec, ub = ub, sense = sense, vtype = vtype, n = 1,
                         control = list(trace = trace, round = round_cplex, tilim = t_max))
    time = (proc.time()-ptm)[3L]

    if (out$status == 108) {
      message("  Error: time limit exceeded, no integer solution!")
      obj_total = NA
      obj_dist_mat = NA
      t_id = NA
      c_id = NA
      group_id = NA
      time = NA
    }
    else if (is.na(out$obj)) {
      message("  Error: problem infeasible!")
      obj_total = NA
      obj_dist_mat = NA
      t_id = NA
      c_id = NA
      group_id = NA
      time = NA
    }
    else {
      message("  Optimal matches found")

      if (approximate == 1) {
        rel = .relaxation_b(n_t, n_c, out$xopt, dist_mat, subset_weight, "cplex",
                            round_cplex, trace)
        out$xopt = rel$sol
        out$obj = rel$obj
        time = time + rel$time
      }

      #! Matched controls indexes
      t_id = unique(sort(rep(1:n_t, n_c))[out$xopt[1:(n_t*n_c)]==1])
      c_id = (c_index+n_t)[out$xopt[1:(n_t*n_c)]==1]

      #! Group (or pair) identifier
      group_id_t = 1:(length(t_id))
      group_id_c = sort(rep(1:(length(t_id)), n_controls))
      group_id = c(group_id_t, group_id_c)

      #! Optimal value of the objective function
      obj_total = out$obj

      if (!is.null(dist_mat)) {
        obj_dist_mat = sum(c(as.vector(matrix(t(dist_mat), nrow = 1, byrow = TRUE)) * out$xopt[1:(n_t*n_c)]==1))
      }
      else {
        obj_dist_mat = NULL
      }
    }
  }

  #! GLPK
  else if (solver == "glpk") {
    if (!requireNamespace('Rglpk', quietly = TRUE)) {
      stop('Required solver not installed')
    }

    message("  GLPK optimizer is open...")
    dir = rep.int(NA_character_, length(prmtrs$sense))
    dir[prmtrs$sense=="E"] = '=='
    dir[prmtrs$sense=="L"] = '<='
    dir[prmtrs$sense=="G"] = '>='

    bound = list(lower = list(ind=c(1:length(ub)), val=rep(0,length(ub))),
                 upper = list(ind=c(1:length(ub)), val=ub))

    message("  Finding the optimal matches...")
    ptm = proc.time()
    out= Rglpk::Rglpk_solve_LP(cvec, Amat, dir, bvec, bounds = bound, types = vtype, max = FALSE)
    time = (proc.time()-ptm)[3]

    if (out$status != 0) {
      message("  Error: problem infeasible!")
      obj_total = NA
      obj_dist_mat = NA
      t_id = NA
      c_id = NA
      group_id = NA
      time = NA
    }
    else {
      message("  Optimal matches found")

      if (approximate == 1) {
        rel = .relaxation_b(n_t, n_c, out$solution, dist_mat, subset_weight, "glpk", round_cplex, trace)
        out$solution = rel$sol
        out$optimum = rel$obj
        time = time + rel$time
      }


      #! Matched controls indexes
      t_id = unique(sort(rep(1:n_t, n_c))[out$solution[1:(n_t*n_c)]==1])
      c_id = (c_index+n_t)[out$solution[1:(n_t*n_c)]==1]

      #! Group (or pair) identifier
      group_id_t = 1:(length(t_id))
      group_id_c = sort(rep(1:(length(t_id)), n_controls))
      group_id = c(group_id_t, group_id_c)

      #! Optimal value of the objective function
      obj_total = out$optimum

      if (!is.null(dist_mat)) {
        obj_dist_mat = sum(c(as.vector(matrix(t(dist_mat), nrow = 1, byrow = TRUE)) * out$solution[1:(n_t*n_c)]==1))
      } else {
        obj_dist_mat = NULL
      }

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

    bound = list(lower = list(ind=c(1:length(ub)), val=rep(0,length(ub))),
                 upper = list(ind=c(1:length(ub)), val=ub))

    message("  Finding the optimal matches...")
    ptm = proc.time()
    out= Rsymphony::Rsymphony_solve_LP(cvec, Amat, dir, bvec, bounds = bound, types = vtype, max = FALSE, time_limit = t_max)
    time = (proc.time()-ptm)[3]

    if (out$status == 0) {
      message("  Optimal matches found")

      if (approximate == 1) {
        rel = .relaxation_b(n_t, n_c, out$solution, dist_mat, subset_weight, "symphony", round_cplex, trace)
        out$solution = rel$sol
        out$objval = rel$obj
        time = time + rel$time
      }

      #! Matched controls indexes
      t_id = unique(sort(rep(1:n_t, n_c))[out$solution[1:(n_t*n_c)]==1])
      c_id = (c_index+n_t)[out$solution[1:(n_t*n_c)]==1]

      #! Group (or pair) identifier
      group_id_t = 1:(length(t_id))
      group_id_c = sort(rep(1:(length(t_id)), n_controls))
      group_id = c(group_id_t, group_id_c)

      #! Optimal value of the objective function
      obj_total = out$objval

      if (!is.null(dist_mat)) {
        obj_dist_mat = sum(c(as.vector(matrix(t(dist_mat), nrow = 1, byrow = TRUE)) * out$solution[1:(n_t*n_c)]==1))
      }
      else {
        obj_dist_mat = NULL
      }
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
      t_id = NA
      c_id = NA
      group_id = NA
      time = NA
    }
  }

  #! Output
  list(obj_total = obj_total,
       obj_dist_mat = obj_dist_mat,
       t_id = t_id,
       c_id = c_id,
       group_id = group_id,
       time = time,
       status = out$status)
}
