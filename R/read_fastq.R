read_fastq <- function(file) {
  assertthat::assert_that(assertthat::is.readable(file))
  assertthat::assert_that(assertthat::has_extension(file,"fq"))

  base::scan(file, character()) -> file.lines
  file.lines[c(T,F,F,F)] -> ids
  file.lines[c(F,T,F,F)] -> sequences
  file.lines[c(F,F,F,T)] -> qualities

  if (!all(base::startsWith(ids,"@"))) {
    stop("Some ID lines didn't start with @")
  }

  stringr::str_sub(ids,2) -> ids

  if (!all(base::nchar(sequences)==base::nchar(qualities))) {
    base::stop("Some sequences were a different length to the qualities")
  }

  if (any(base::duplicated(ids))) {
    base::stop("Some IDs are duplicated")
  }

  tibble::tibble(
    ID = ids,
    Bases=sequences,
    Qualities=qualities,
    GC=gc_content(sequences)
  ) %>%
    return()

}

gc_content <- function(seq) {

  assertthat::assert_that(is.character(seq))

  seq <- base::toupper(seq)

  if (any(stringr::str_detect(seq,"[^GATC]"))) {
    base::warning("Non GATC characters found in sequences")
  }

  stringr::str_replace_all(seq,"[^GC]","") -> just_gc

  return(100*(base::nchar(just_gc)/base::nchar(seq)))

}


