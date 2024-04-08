# Tests de base.

test_that("Test de base pour meaps continu", {

M1 <- matrix(1:16 +.1, nrow = 4 )
emp1 <- rep(10., 4)
act1 <- rep(11., 4)
f1 <- (sum(act1) / sum(emp1)) - 1
f1 <- rep(f1,4)
shuf1 <- matrix(c(1:4), nrow = 1)


.f_distrib <- function(emp, f, act, attract = 1L) {
  
  accessibiliy <- cumsum(emp * attract)
  absorption <- -log(f)/max(accessibiliy)
  passants <- act * exp(-absorption * accessibiliy)
  
  dplyr::lag(passants, default = act) - passants
}

ligne1 <- .f_distrib(emp1, f1, act1[1])
ligne2 <- .f_distrib(emp1-ligne1, f1, act1[2])
ref1.1 <- one_distrib_continu(act1[1], f1[1], rep(1.,4), M1[1, ], emp1)

expect_equal(ref1.1, ligne1)

res1 <- meaps_continu(M1, emp1, act1, f1, shuf1, nthreads = 1)

expect_equal(ligne1, res1[1,])
expect_equal(ligne2, res1[2,])

})     


test_that("Comparaison meaps_continu avec meaps_multishuf", {
  d <- genere_data(n=4, k=4, nshuf=10, densite="uniforme", tx_nas = 0, dist_not_equal = TRUE)
 
  
  meaps_c <- meaps_continu(
    dist=d$dist,
    emplois=d$emplois,
    actifs=d$actifs, 
    attraction = "constant",
    f=d$fuite,
    shuf=d$shuf,
    progress = FALSE, nthreads = 1L) |> as.matrix()
  
  meaps_ref <- meaps_multishuf(
    rkdist=d$rkdist,
    emplois=d$emplois,
    actifs=d$actifs, 
    modds=d$modds, 
    f=d$fuite,
    shuf=d$shuf,
    progress = FALSE)
  
  
  expect_true(all(meaps_c >= 0))
  expect_equal(sum(meaps_c), sum(d$emplois))
  expect_equal(rowSums(meaps_c), d$actifs*(1-d$fuite))
  expect_equal(colSums(meaps_c), d$emplois)
  expect_equal(meaps_c, meaps_ref, tolerance = 0.01)
})

