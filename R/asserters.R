neq_null           <- function(x) !is.null(x)
neq_empt_chr       <- function(x) !identical(x, character(0))
neq_empt_num       <- function(x) !identical(x, numeric(0))
neq_empt_int       <- function(x) !identical(x, integer(0))
neq_empt_lst       <- function(x) !identical(x, list())
is_matrix_chr      <- function(x) is.matrix(x) && typeof(x) == "character"
is_data_frame_chr  <- function(x) is.data.frame(x) && all(.map_lgl(x, is.character))
identical_colnames <- function(x, y) identical(colnames(x), colnames(y))

only_single_chars  <- function(A) {
  for (i in seq_along(nrow(A))) {
    for (j in seq_along(ncol(A)))
      if (nchar(A[i,j]) != 1L) return(FALSE)
  }
  return(TRUE)
}
