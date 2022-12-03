
### A real data example of chloroplast accumulation responses to a blue microbeam ###

# Load packages
library(cellssm)

# Load data of chloroplast movements
data("cell1", "cell2", "cell3", "cell4", "visual")
cell_list <- list(cell1, cell2, cell3, cell4)

# Execution of state-space modeling

# When you want to compare the statistical and visual estimation of the start time
ssm_KFAS(cell_list = cell_list, visual = visual, out = "03_ssm_KFAS",
         res_name = "chloroplast", ex_name = "microbeam",
         unit1 = "micrometer", unit2 = "min")

# Test
if(requireNamespace("KFAS", quietly = TRUE)){
  test_that("length of the output is correct", {
    expect_equal(length(list.files("03_ssm_KFAS/csv")), 52)
    expect_equal(length(list.files("03_ssm_KFAS/pdf")), 50)
  })
}

unlink("03_ssm_KFAS", recursive = T)



### A simulated data example of Paramecium escape responses from a laser heating ###

# Load data
data("Paramecium")
cell_list <- list(Paramecium)


# Execution of state-space modeling
ssm_KFAS(cell_list = cell_list, out = "13_ssm_KFAS",
         ex_sign = "positive", df_name = "experiment",
         res_name = "Paramecium", ex_name = "heat",
         unit1 = "millimeter", unit2 = "sec")

# Test
if(requireNamespace("KFAS", quietly = TRUE)){
  test_that("length of the output is correct", {
    expect_equal(length(list.files("13_ssm_KFAS/csv")), 12)
    expect_equal(length(list.files("13_ssm_KFAS/pdf")), 10)
  })
}

unlink("13_ssm_KFAS", recursive = T)

