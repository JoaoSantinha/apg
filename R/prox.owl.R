#' Proximal operator of the order weighted l1-norm
#'
#' Computes the proximal operator of the ordered weighted l1-norm
#'
#' @param x The input vector
#'
#' @return The projection of \code{x} onto the set of nondecreasing vectors,
#'   obtained by solving an isotonic regression.
#'
#' @export
#' @examples prox.isotonic(c(1,3,-2,4,5))


prox.owl <- function(x, t, opts=list()) {
  v <- x
  
  w <- opts$weights
  
  v.abs <- abs(v)
  sorting <- sort(v.abs, decreasing = TRUE, index.return= TRUE)
  ix <- sorting$ix
  v.abs <- v.abs[ix]
  v.abs <- pava(v.abs - w, decreasing = FALSE)
  v.abs[v.abs < 0] <- 0 
  
  # undo the sorting
  inv.ix <- phonTools::zeros(ix)
  inv.ix[ix] <- seq(length(v))
  v_abs <- v_abs[inv.ix]
  
  return (sign(v) * v_abs)
}
