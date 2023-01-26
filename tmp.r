library(seqinr)
library(dplyr)
library(stringr)

dengue <- read.fasta('dengue.fasta')
dengue <- dengue[[1]]

myrho <- function(sequence='', astr=''){
  if (is.null(sequence) | is.null(astr))
  {
    print("invalid input")
    return (0)
  }
  srch_tbl <- seqinr::count(sequence, wordsize = str_length(astr))
  monomer_tbl <- seqinr::count(sequence, wordsize = 1)
  srch_mono <- str_split(tolower(astr), '')
  srch_mono <- srch_mono[[1]]
  
  denom <- 1
  for (i in 1:length(srch_mono)) {
    denom <- denom * monomer_tbl[srch_mono[i]]/length(sequence)
  }
  numer <- srch_tbl[tolower(astr)]/length(sequence)
  return(numer/denom)
}

res <- myrho(dengue, 'CAG')
print(res)