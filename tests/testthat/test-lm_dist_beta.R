
### A real data example of chloroplast accumulation responses to a blue microbeam ###

# Load packages
library(cellssm)

## Bayes
# Load data
data("cell1", "cell2", "cell3", "cell4", "chloroplast_mvtime")
cell_list <- list(cell1, cell2, cell3, cell4)

# Specify the path of the output directory of [ssm_individual] or [ssm_KFAS]
# below, the path of system files is specified to show an example
ssm_file <- stringr::str_split(system.file("extdata", "individual_model.stan",
                                           package = "cellssm"), "/")[[1]]
ssm_path <- paste(ssm_file[-length(ssm_file)], collapse = "/")

# Linear regression
glist <- lm_dist_beta(cell_list = cell_list, mvtime = chloroplast_mvtime,
                      ssm_path = ssm_path,
                      ssm_method = "Bayes", res_name = "chloroplast",
                      ex_name = "microbeam", unit1 = "micrometer", unit2 = "min")

# Test
test_that("length and class of the output are correct", {
  expect_equal(length(glist), 5)
  expect_equal(class(glist[[1]])[2], "ggplot")
})


