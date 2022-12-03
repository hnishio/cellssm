
### A real data example of chloroplast accumulation responses to a blue microbeam ###

# Load packages
library(cellssm)

# Load data of chloroplast movements
data("cell1", "visual")
cell_list <- list(cell1[,1:5])

# Execution of state-space modeling

# Test
test_that("length of the output is correct", {
  skip_if_not_installed("KFAS")
  ssm_KFAS(cell_list = cell_list, visual = visual[1:3,], out = "03_ssm_KFAS",
           res_name = "chloroplast", ex_name = "microbeam",
           unit1 = "micrometer", unit2 = "min")
  expect_equal(length(list.files("03_ssm_KFAS/csv")), 5)
  expect_equal(length(list.files("03_ssm_KFAS/pdf")), 3)
})

unlink("03_ssm_KFAS", recursive = T)



### A simulated data example of Paramecium escape responses from a laser heating ###

# Load data
data("Paramecium")
cell_list <- list(Paramecium[,1:5])


# Execution of state-space modeling

# Test
  test_that("length of the output is correct", {
    skip_if_not_installed("KFAS")
    ssm_KFAS(cell_list = cell_list, out = "13_ssm_KFAS",
             ex_sign = "positive", df_name = "experiment",
             res_name = "Paramecium", ex_name = "heat",
             unit1 = "millimeter", unit2 = "sec")
    expect_equal(length(list.files("13_ssm_KFAS/csv")), 5)
    expect_equal(length(list.files("13_ssm_KFAS/pdf")), 3)
  })

unlink("13_ssm_KFAS", recursive = T)

