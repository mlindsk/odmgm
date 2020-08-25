## ---------------------------------------------------------
##                NON-EXPORTED HELPERS
## ---------------------------------------------------------

## MAPS
.map_chr     <- function(x, fun, ...) vapply(X = x, FUN = fun, FUN.VALUE = character(1), ...)
.map_int     <- function(x, fun, ...) vapply(X = x, FUN = fun, FUN.VALUE = integer(1), ...)
.map_dbl     <- function(x, fun, ...) vapply(X = x, FUN = fun, FUN.VALUE = numeric(1), ...)
.map_lgl     <- function(x, fun, ...) vapply(X = x, FUN = fun, FUN.VALUE = logical(1), ...)


## STRINGS
es_to_vs     <- function(e) strsplit(e, "\\|")
vs_to_es     <- function(e) lapply(e, paste0, collapse = "|")
rev_es       <- function(e) .map_chr(es_to_vs(e), function(x) paste0(rev(x), collapse = "|"))
sort_        <- function(x) paste0(sort(x), collapse = "|")
.split_chars <- function(x) unlist(strsplit(x, ""))

## MISC
`%[int%`     <- function(a, x) if (neq_empt_int(x)) return(a[x]) else return(a)
`%[chr%`     <- function(a, x) if (neq_empt_chr(x)) return(a[x]) else return(0)
# `%[[chr_df%` <- function(a, x) if (neq_empt_chr(x)) return(a[, x]) else return(0)
push         <- function(l, el, name = NULL) c(l, structure(list(el), names = name))
## '%ni%'       <- Negate('%in%')

## ---------------------------------------------------------
##                  EXPORTED HELPERS
## ---------------------------------------------------------

#' To Single Chars
#'
#' Convert all values in a data frame or matrix of characters to a single character representation
#'
#' @param x Data frame or matrix of characters
#' @examples
#'
#' ## Example 1 - pure discrete
#' d1 <- data.frame(x = c("hhh", "f"), y = c("f", "hhh"), stringsAsFactors = FALSE)
#' to_single_chars(d1)
#'
#' ## Example 2 - mixed variables
#' library(dplyr)
#' d2 <- cbind(d1, data.frame(z = runif(2), w = runif(2)))
#' 
#' d2_disc <- d2 %>%
#'         select_if(is.character) %>%
#'         to_single_chars()
#'
#' d2_cont <- d2 %>% select_if(is.numeric)
#'
#' cbind(d2_disc, d2_cont)
#' @export
to_single_chars <- function(x) {
  ## Implicitly assumes that no columns has more than length(letters) = 26 unique levels
  apply(x, 2, function(z) {
    f <- as.factor(z)
    levels(f) <- letters[1:length(levels(f))]
    as.character(f)
  })
}
