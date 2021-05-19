library(fastqR)
library(testthat)

test_that("Reading FastQ files", {

  expect_output(
    str(
      suppressWarnings(
        read_fastq(system.file("good.fq", package = "fastqR"))
      )),
    "tibble"
  )
  expect_warning(read_fastq(system.file("good.fq", package = "fastqR")))
  expect_error(read_fastq(system.file("broken_format.fq", package = "fastqR")))
  expect_error(read_fastq(system.file("duplicate_ids.fq", package = "fastqR")))
  expect_error(read_fastq(system.file("mismatched_lengths.fq", package = "fastqR")))
  expect_error(read_fastq(system.file("wrong_extension.fa", package = "fastqR")))
})


test_that("Check GC content", {
  expect_equal(gc_content("GGATCTTAGG"), 50)
  expect_warning(gc_content("GGATCNTTAGG"))
})


test_that("Check Quality Decoding", {
  expect_equal(decode_qualities("B"),33)
  expect_equal(decode_qualities("B", offset=64),2)
  expect_equal(decode_qualities("AB"),c(32,33))
  expect_error(decode_qualities(c("A","B")))
  expect_equal(decode_qualities("@"),31)
  expect_error(decode_qualities("@",offset=64))
  expect_error(decode_qualities("A",offset="sanger"))
  expect_error(decode_qualities("A",offset=63))

})
