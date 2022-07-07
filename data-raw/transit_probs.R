
# default transition probabilities ====
tp0 <- readRDS("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/code/package_functions/transitProbs_27subtypes_1350cells_Luo2017&Liu2021.rds")

usethis::use_data(tp0, overwrite = TRUE, internal = TRUE)
