#' Proximal operator of the order weighted l1-norm
#'
#' Computes the proximal operator of the ordered weighted l1-norm
#'
#' @param v The input vector
#'
#' @return The projection of \code{v} onto the set of nondecreasing vectors,
#'   obtained by solving an isotonic regression.
#'
#' @export
#' @examples prox.isotonic(c(1,3,-2,4,5))


prox.owl <- function(v, t=0, opts=list()) {
  # v <- x
  
  w <- opts$weights
  v_abs <- abs(v)

  sorting <- sort(v_abs, decreasing = TRUE, index.return= TRUE)
  ix <- sorting$ix

  v_abs <- v_abs[ix]
  
  # print(ajsjgahshg)
  v_abs <- Iso::pava(v_abs - w, decreasing = TRUE)
  v_abs[v_abs < 0] <- 0

  # undo the sorting
  inv.ix <- phonTools::zeros(ix)
  inv.ix[ix] <- seq(length(v))
  v_abs <- v_abs[inv.ix]

  return (sign(v) * v_abs)
}

