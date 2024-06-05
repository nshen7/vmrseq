pkgname <- "vmrseq"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('vmrseq')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("cell_1")
### * cell_1

flush(stderr()); flush(stdout())

### Name: cell_1
### Title: cell_1
### Aliases: cell_1
### Keywords: datasets

### ** Examples

data(cell_1)
cell_1




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
