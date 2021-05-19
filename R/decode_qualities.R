#' Decode FastQ quality strings into Phred scores
#'
#' @param qualities a string scalar containing FastQ encoded quality scores
#' @param offset the ASCII offset (33 or 64) of the encoded phred scores
#'
#' @return a vector of Phred scores
#' @export
#'
#' @examples
#' decode_qualities("???#;ABAAAH")
#'
#' decode_qualities("WZZZVX[]", offset=64)
decode_qualities <- function(qualities, offset=33) {

  assertthat::assert_that(assertthat::is.scalar(qualities))
  assertthat::assert_that(assertthat::is.string(qualities))
  assertthat::assert_that(assertthat::is.scalar(offset))
  assertthat::assert_that(assertthat::is.number(offset))

  if (!(offset==33 | offset==64)) {
    stop("Offset can only be 33 or 64")
  }

  as.integer(charToRaw(qualities)) - offset -> phred_scores

  if (any(phred_scores < 1)) {
    base::stop("Negative phred scores produced - check offset")
  }

  return(phred_scores)
}
