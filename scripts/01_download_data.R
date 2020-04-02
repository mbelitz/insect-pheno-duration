# Script to download data

source("scripts/00_pkg_functions.R")

data_fun <- function(binomial) {
  dat <- get_records(binomial)
  
  write.csv(dat, paste("downloaded_data/", binomial, ".csv", 
                       sep = "", row.names(FALSE)))
  
  print(paste("Downloaded", binomial, sep = " "))
}

#Download Odonata data

odes <- list(
"Ischnura verticalis",
"Sympetrum corruptum",
"Argia vivida",
"Pachydiplax longipennis",
"Enallagma civile",
"Hetaerina americana",
"Libellula lydia",
"Calopteryx maculata",
"Anax junius")

mclapply(odes, data_fun, mc.cores = 9)

# Download Bee data

bees <- list(
  "Bombus impatiens",
  "Bombus pensylvanicus",
  "Bombus fervidus",
  "Bombus vagans",
  "Bombus griseocollis",
  "Xylocopa virginica",
  "Bombus melanopygus",
  "Bombus terricola",
  "Bombus flavifrons"
)

mclapply(bees, data_fun, mc.cores = 10)
