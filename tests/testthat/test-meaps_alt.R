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
  #expect_snapshot(meapst)
  expect_equal(sum(meapst), sum(d$emplois))
  expect_equal(rowSums(meapst), d$actifs*(1-d$fuite))
  expect_equal(colSums(meapst), d$emplois)
  expect_equal(sum(!abs(meapss-meapst)<0.0001) < 0.001*length(meapss), TRUE)
})

# test avec des nas

test_that("meaps alt avec NA", {
  d <- genere_data(n=8, k=8, nshuf=1, densite="uniforme")
  d$dist[d$dist>1.0] <- NA
  d$rkdist <- matrixStats::rowRanks(d$dist, ties = "first")
  meaps <- meaps_oneshuf(
    rkdist=d$rkdist,
    emplois=d$emplois,
    actifs=d$actifs, 
    modds=d$modds, 
    f=d$fuite,
    shuf=d$shuf)
  meapst <- meaps_alt(
    rkdist=d$rkdist,
    emplois=d$emplois,
    actifs=d$actifs, 
    modds=d$modds, 
    f=d$fuite,
    shuf=d$shuf,
    progress=FALSE)
  expect_equal(meapst, meaps)
  expect_equal(all(abs(rowSums(meapst)-d$actifs*(1-d$fuite))<0.0001), TRUE)
  expect_equal(all(abs(colSums(meapst)-d$emplois)<0.0001), TRUE)
})

# test des odds grands

test_that("meaps alt avec NA", {
  d <- genere_data(n=4, k=4, nshuf=1, densite="uniforme")
  d$rkdist <- matrixStats::rowRanks(d$dist, ties = "first")
  d$modds[1,1] <- 1000
  meapst <- meaps_alt(
    rkdist=d$rkdist,
    emplois=d$emplois,
    actifs=d$actifs, 
    modds=d$modds, 
    f=d$fuite,
    shuf=d$shuf,
    progress=FALSE)
  expect_equal(rowSums(meapst), d$actifs*(1-d$fuite))
  expect_equal(colSums(meapst), d$emplois)
})

# test des odds petits

test_that("meaps alt avec NA", {
  d <- genere_data(n=8, k=8, nshuf=1, densite="uniforme")
  d$rkdist <- matrixStats::rowRanks(d$dist, ties = "first")
  d$modds[1,1] <- 0.0001
  meapst <- meaps_alt(
    rkdist=d$rkdist,
    emplois=d$emplois,
    actifs=d$actifs, 
    modds=d$modds, 
    f=d$fuite,
    shuf=d$shuf,
    progress=FALSE)
  expect_equal(rowSums(meapst), d$actifs*(1-d$fuite))
  expect_equal(colSums(meapst), d$emplois)
})