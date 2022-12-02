
### A real data example of chloroplast accumulation responses to a blue microbeam ###

# Load packages
library(cellssm)

# Load data of chloroplast movements
data("cell1", "cell2", "cell3", "cell4", "visual")
cell_list <- list(cell1, cell2, cell3, cell4)

# Predict movement

# When you want to compare the statistical and visual estimation of the start time
nomodel(cell_list = cell_list, visual = visual, out = "04_nomodel",
        res_name = "chloroplast", ex_name = "microbeam",
        unit1 = "micrometer", unit2 = "min")

# Test
test_that("length of the output is correct", {
  expect_equal(length(list.files("04_nomodel")), 51)
})

unlink("04_nomodel", recursive = T)



### A simulated data example of Paramecium escape responses from a laser heating ###

# Load data
data("Paramecium")
cell_list <- list(Paramecium)

# Predict movement
nomodel(cell_list = cell_list, out = "14_nomodel",
        ex_sign = "positive", fold = 1, df_name = "experiment",
        res_name = "Paramecium", ex_name = "heat",
        unit1 = "millimeter", unit2 = "sec")

# Test
test_that("length of the output is correct", {
  expect_equal(length(list.files("14_nomodel")), 11)
})

unlink("14_nomodel", recursive = T)

