#' Get optimal span value
#'
#' Get optimal span value for QC-RLSC normalization
#'
#' @param qc_data data.frame of only the QC samples from one block/batch, where
#'   the first column gives the run order and the second column gives the
#'   abundance value
#' @param lam vector of (numeric) possible polynomial degrees to test in the
#'   loess fitting; defaults to 1,
#' @return list containing elements for: Lambda, Alpha, MSE
#' @details This internal function is called by the \code{get_params} function
#' @author Lisa Bramer, Kelly Stratton
#' @rdname optimal_span
#' @name optimal_span
optimal_span <- function(qc_data, lam = c(1, 2)) {
  ### initial checks ###

  # check that qc_data is data.frame with 2 columns #
  if (class(qc_data) != "data.frame") {
    stop("qc_data must be a data.frame")
  }
  if (ncol(qc_data) != 2) {
    stop("qc_data must have 2 columns")
  }

  # check that lam is numeric vector #
  if (!is.numeric(lam)) {
    stop("lam must be a numeric vector")
  }
  if (!is.vector(lam)) {
    stop("lam must be a numeric vector")
  }

  ### end of initial checks ###

  qc_data = qc_data[!is.na(qc_data[, 2]), ] # remove QCs that have NA values for abundance

  # determine number of QC samples #
  n_qc = nrow(qc_data) - 1

  # initialize results list #
  all_res = list()

  # iterate over polynomial degree options #
  for (i in lam) {
    # set current lambda value #
    cur_lam = i

    # create a sequence of possible alpha values to be considered #
    alp_vals = seq(from = (cur_lam + 1),
                   to = n_qc,
                   by = 1) / n_qc

    res = matrix(NA, nrow = (n_qc + 1), ncol = length(alp_vals))

    for (j in 1:length(alp_vals)) {
      cur_alp = alp_vals[j]
      pred = NULL

      for (k in 1:(n_qc + 1)) {
        # LOO-CV
        data = qc_data[-(k), ]
        model = loess(
          data[, 2] ~ data[, 1],
          span = cur_alp,
          degree = cur_lam,
          control = loess.control(surface = "direct")
        )
        pred[k] = predict(model, newdata = qc_data[k, 1])
      }
      res[, j] = pred
    }
    all_res[[i]] = res
  }
  temp = lapply(all_res, function(x) {
    apply((x - qc_data[, 2]) ^ 2, 2, mean)
  })

  res = list()

  for (m in 1:length(temp)) {
    temp2 = temp[[m]]
    alp = seq(from = lam[m] + 1,
              to = n_qc,
              by = 1) / n_qc
    res[[m]] = data.frame(Lambda = rep(lam[m]),
                          Alpha = alp,
                          MSE = temp2)
  }

  final = do.call(rbind, res)
  final2 = final[which.min(final$MSE),]
  return(list(
    Lambda = final2$Lambda,
    Alpha = final2$Alpha,
    MSE = final2$MSE,
    all_res = final
  ))
}
