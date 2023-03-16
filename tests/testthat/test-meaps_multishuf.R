test_that("meaps multishuf", {
  d <- genere_data(n=10, k=10, nshuf=16, densite="uniforme")
  
  meapst <- meaps_multishuf(
    rkdist=d$rkdist,
    emplois=d$emplois,
    actifs=d$actifs, 
    modds=d$modds, 
    f=d$fuite,
    shuf=d$shuf,
    progress = FALSE)
  
  expect_equal(sum(meapst), sum(d$emplois))
  expect_equal(rowSums(meapst), d$actifs*(1-d$fuite))
  expect_equal(colSums(meapst), d$emplois)
})

# test avec des nas. en raison de la saturation peut s'écarter nettement pour les fins de shufs.

test_that("meaps multishuf avec NAs", {
  d <- genere_data(n=8, k=8, nshuf=50, densite="uniforme")
  d$dist[d$dist>1.0] <- NA
  d$rkdist <- matrixStats::rowRanks(d$dist, ties = "first")
  
  meapst <- meaps_multishuf(
    rkdist=d$rkdist,
    emplois=d$emplois,
    actifs=d$actifs, 
    modds=d$modds, 
    f=d$fuite,
    shuf=d$shuf,
    progress=FALSE)
  
  expect_equal(all(abs(rowSums(meapst)-d$actifs*(1-d$fuite)) / d$actifs*(1-d$fuite) < 0.05), TRUE)
  expect_equal(all(abs(colSums(meapst)-d$emplois) / d$emplois <0.05), TRUE)
})

# test des grands odds 

test_that("meaps multishuf avec NAs et grand OR", {
  d <- genere_data(n=4, k=4, nshuf=1, densite="uniforme")
  d$rkdist <- matrixStats::rowRanks(d$dist, ties = "first")
  d$modds[1,1] <- 1000
  
  meapst <- meaps_multishuf(
    rkdist=d$rkdist,
    emplois=d$emplois,
    actifs=d$actifs, 
    modds=d$modds, 
    f=d$fuite,
    shuf=d$shuf,
    progress=FALSE)
  
  expect_equal(all(abs(rowSums(meapst)-d$actifs*(1-d$fuite))<0.0001), TRUE)
  expect_equal(all(abs(colSums(meapst)-d$emplois)<0.0001), TRUE)
})

# test des odds petits

test_that("meaps multishuf avec NAs et petit OR", {
  d <- genere_data(n=8, k=8, nshuf=1, densite="uniforme")
  d$rkdist <- matrixStats::rowRanks(d$dist, ties = "first")
  d$modds[1,1] <- 0.0001
  
  meapst <- meaps_multishuf(
    rkdist=d$rkdist,
    emplois=d$emplois,
    actifs=d$actifs, 
    modds=d$modds, 
    f=d$fuite,
    shuf=d$shuf,
    progress=FALSE)
  
  expect_equal(all(abs(rowSums(meapst)-d$actifs*(1-d$fuite))<0.0001), TRUE)
  expect_equal(all(abs(colSums(meapst)-d$emplois)<0.0001), TRUE)
})