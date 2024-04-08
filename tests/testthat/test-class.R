test_that("Test sur la classe RankedRSMatrix", {

library(Matrix)

x <- c(1.2, 2.9, 2.0) ; xr <- c(1.2, 2.0, 2.9) 
j <- c(1L, 0L, 2L) ; jr <- c(1L, 2L, 0L)
p <- c(0L, 1L, 3L)
dim <- c(2L,3L)
  
mat <- new(RankedRSMatrix, xr, jr, p, dim) # constructeur sur des données rangées.
expect_equal(mat$xr, xr)
expect_equal(mat$jr, jr)
expect_equal(mat$p, p)
expect_equal(mat$dim, dim)
expect_equal(mat$nrow(), dim[1])
expect_equal(mat$ncol(), dim[2])
expect_equal(mat$ncol(), dim[2])
expect_equal(mat$at(0L,0L), NA_real_)
expect_equal(mat$at(0L,1L), x[1])

# ATTENTION le j suit la convention c++ de commencer à compter de 0. Valide pour RankedRSMatrix, mais pas pour dgRMatrix
m2 <- sparseMatrix(x = x, j = j + 1L, p = p, dims = dim, repr = "R")

expect_equal(mat$unrank(), m2)

m3 <- new(RankedRSMatrix, m2) # constructeur sur des données NON rangées.

expect_equal(m3, mat)

# tri <- mat2triplet(m2) |> as.data.frame()
# sparseMatrix(i = tri$i, j = tri$j, x = tri$x, repr = "R")




})

