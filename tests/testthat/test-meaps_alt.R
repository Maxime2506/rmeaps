test_that("meaps alt", {
  d <- genere_data(n=10, k=10, nshuf=16, densite="uniforme")
  
  meapss <- meaps_multishuf(
    rkdist=d$rkdist,
    emplois=d$emplois,
    actifs=d$actifs, 
    modds=d$modds, 
    f=d$fuite,
    shuf=d$shuf,
    progress = FALSE)
  meapst <- meaps_alt(
    rkdist=d$rkdist,
    emplois=d$emplois,
    actifs=d$actifs, 
    modds=d$modds, 
    f=d$fuite,
    shuf=d$shuf,
    progress = FALSE)
  expect_snapshot(meapst)
  expect_equal(sum(meapst), sum(d$emplois))
  expect_equal(rowSums(meapst), d$actifs*(1-d$fuite))
  expect_equal(colSums(meapst), d$emplois)
  expect_equal(sum(!abs(meapss-meapst)<0.0001) < 0.001*length(meapss), TRUE)
})
