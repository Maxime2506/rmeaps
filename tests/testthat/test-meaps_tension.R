test_that("meaps_multishuf==_tension sur un échantillon moyen", {
  
  d <- genere_data(n=20, k=20, densite="uniforme")
  
  meapss <- meaps_multishuf(
    rkdist=d$rkdist,
    emplois=d$emplois,
    actifs=d$actifs, 
    modds=d$modds, 
    f=d$fuite,
    shuf=matrix(1:length(d$actifs), nrow=8, ncol=length(d$actifs), byrow = TRUE),
    progress = FALSE)
  meapst <- meaps_tension(
    rkdist=d$rkdist,
    emplois=d$emplois,
    actifs=d$actifs, 
    modds=d$modds, 
    f=d$fuite,
    shuf=matrix(1:length(d$actifs), nrow=8, ncol=length(d$actifs), byrow = TRUE),
    progress = FALSE)
  
  expect_equal(meapss, meapst$flux)
})
