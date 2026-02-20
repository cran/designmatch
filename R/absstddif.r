absstddif <- function(X_mat, t_ind, std_dif) {
     n_vrs <- ncol(X_mat)
     n_obs <- nrow(X_mat)

     vapply(seq_len(n_vrs), function(j) {
         yes_before <- unlist(X_mat[t_ind == 1, j])
         no_before <- unlist(X_mat[t_ind == 0, j])
         pooled_sd <- sqrt((var(yes_before, na.rm = TRUE) + var(no_before, na.rm = TRUE))/2)
         pooled_sd * std_dif
     }, numeric(1L))
}
