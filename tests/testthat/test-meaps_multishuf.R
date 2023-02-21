test_that("meaps_oneshuf et _multishuf sur un échantillon moyen", {
  
  d <- genere_data(n=10, k=10, densite="uniforme")
  
  meaps1 <- meaps_oneshuf(
    rkdist=d$rkdist,
    emplois=d$emplois,
    actifs=d$actifs, 
    modds=d$modds, 
    f=d$fuite,
    shuf=1:length(d$actifs))
  
  meapss <- meaps_multishuf(
    rkdist=d$rkdist,
    emplois=d$emplois,
    actifs=d$actifs, 
    modds=d$modds, 
    f=d$fuite,
    shuf=matrix(1:length(d$actifs), nrow=8, ncol=length(d$actifs), byrow = TRUE),
    progress = FALSE)
  
  expect_equal(sum(meapss), sum(d$emplois))
  expect_equal(rowSums(meapss), d$actifs*(1-d$fuite))
  expect_equal(colSums(meapss), d$emplois)
  expect_equal(meaps1, meapss)
})

test_that("meaps_oneshuf et _multishuf sur un échantillon moyen", {
  
  d <- genere_data(n=10, k=10, nshuf=128, densite="uniforme")
  
  meapss <- meaps_multishuf(
    rkdist=d$rkdist,
    emplois=d$emplois,
    actifs=d$actifs, 
    modds=d$modds, 
    f=d$fuite,
    shuf=d$shuf,
    nthreads = 1,
    progress = FALSE)

  meapsst <- meaps_multishuf(
    rkdist=d$rkdist,
    emplois=d$emplois,
    actifs=d$actifs, 
    modds=d$modds, 
    f=d$fuite,
    shuf=d$shuf,
    nthreads = 4L,
    progress = FALSE)
    
  expect_equal(sum(meapss), sum(d$emplois))
  expect_equal(rowSums(meapss), d$actifs*(1-d$fuite))
  expect_equal(colSums(meapss), d$emplois)
  expect_equal(meapss, meapsst)
})

