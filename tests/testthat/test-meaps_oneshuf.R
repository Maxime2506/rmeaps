test_that("meaps_oneshuf sur un petit échantillon", {
  
  d <- genere_data(n=2, k=2, densite="uniforme")
  
  meaps1 <- meaps_oneshuf(
    rkdist=d$rkdist,
    emplois=d$emplois,
    actifs=d$actifs, 
    modds=d$modds, 
    f=d$fuite,
    shuf=1:length(d$actifs))
  
  expect_equal(sum(meaps1), sum(d$emplois))
  expect_equal(rowSums(meaps1), d$actifs*(1-d$fuite))
  expect_equal(colSums(meaps1), d$emplois)

  })

test_that("meaps_oneshuf sur un grand échantillon", {
  
  d <- genere_data(n=50, k=25, densite="uniforme")
  
  meaps1 <- meaps_oneshuf(
    rkdist=d$rkdist,
    emplois=d$emplois,
    actifs=d$actifs, 
    modds=d$modds, 
    f=d$fuite,
    shuf=1:length(d$actifs))
  
  expect_equal(abs(sum(meaps1)-sum(d$emplois))<0.0001, TRUE)
  expect_equal(all(abs(rowSums(meaps1)-d$actifs*(1-d$fuite))<0.0001), TRUE)
  expect_equal(all(abs(colSums(meaps1)-d$emplois)<0.0001), TRUE)
  
})

test_that("meaps_oneshuf sur un petit échantillon, timing", {
  
  d <- genere_data(n=20, k=20, densite="uniforme")
  
  bench <- bench::mark(one = meaps_oneshuf(
    rkdist=d$rkdist,
    emplois=d$emplois,
    actifs=d$actifs, 
    modds=d$modds, 
    f=d$fuite,
    shuf=1:length(d$actifs)))
  # on teste qu'on fait plus de 5 itérations par seconde (en monoprocesseur)
  expect_equal(bench[["itr/sec"]]>5, TRUE)
})
