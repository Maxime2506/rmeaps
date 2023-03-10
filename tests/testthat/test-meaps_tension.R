test_that("meaps_multishuf==_tension sur un échantillon moyen", {
  
  d <- genere_data(n=10, k=10, densite="uniforme", ties = "first")
  shuf_u <- rbind(1:length(d$actifs),
                  length(d$actifs):1)
  meapss <- meaps_multishuf(
    rkdist=d$rkdist,
    emplois=d$emplois,
    actifs=d$actifs, 
    modds=d$modds, 
    f=d$fuite,
    shuf=shuf_u,
    progress = FALSE)
  
  meapst <- meaps_tension(
    rkdist=d$rkdist,
    emplois=d$emplois,
    actifs=d$actifs, 
    modds=d$modds, 
    f=d$fuite,
    shuf=shuf_u,
    progress = FALSE)
  
  expect_equal(meapst$flux, meapss)
  expect_snapshot(meapst$tension)
})