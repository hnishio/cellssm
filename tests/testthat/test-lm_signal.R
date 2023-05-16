
### A real data example of chloroplast accumulation responses to a blue microbeam ###

# Load packages
library(cellssm)

# Load data
data("cell1", "cell2", "cell3", "cell4", "chloroplast_mvtime")
cell_list <- list(cell1, cell2, cell3, cell4)

# Linear regression
glist <- lm_signal(cell_list = cell_list, mvtime = chloroplast_mvtime,
                   ex_name = "microbeam", unit1 = "micrometer", unit2 = "min")

# Test
test_that("length and class of the output are correct", {
  expect_equal(length(glist), 6)
  expect_equal(class(glist[[1]])[2], "ggplot")
})



### A simulated data example of Paramecium escape responses from a laser heating ###

# Load data
data("Paramecium", "Paramecium_mvtime")
cell_list <- list(Paramecium)

# Linear regression
glist <- lm_signal(cell_list = cell_list, mvtime = Paramecium_mvtime,
                   graph_title = "Experiment", ex_name = "heat",
                   unit1 = "millimeter", unit2 = "sec")

# Test
test_that("length and class of the output are correct", {
  expect_equal(length(glist), 2)
  expect_equal(class(glist[[1]])[2], "ggplot")
})

