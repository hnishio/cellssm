
### A real data example of chloroplast accumulation responses to a blue microbeam ###

# Load packages
library(cellssm)

# Load real data of chloroplast movements
data("cell1", "cell2", "cell3", "cell4")
cell_list <- list(cell1, cell2, cell3, cell4)

# Plotting
glist <- dist_vis(cell_list = cell_list,
                  res_name = "chloroplast", ex_name = "microbeam",
                  unit1 = "micrometer", unit2 = "min")

# Test
test_that("length and class of the output are correct", {
  expect_equal(length(glist), length(cell_list))
  expect_equal(class(glist[[1]])[2], "ggplot")
})



### A simulated data example of Paramecium escape responses from a laser heating ###

# Load simulated data of Paramecium movement
data("Paramecium")
cell_list <- list(Paramecium)

# Plotting
glist <- dist_vis(cell_list = cell_list,
                  res_name = "Paramecium", ex_name = "heat",
                  unit1 = "millimeter", unit2 = "sec")

# Test
test_that("length and class of the output are correct", {
  expect_equal(length(glist), length(cell_list))
  expect_equal(class(glist[[1]])[2], "ggplot")
})

