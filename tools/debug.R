devtools::clean_dll()
devtools::load_all()

library(tidyverse)
library(Matrix)

loc <- data.frame(  x = rep(-1:1*50, 3), 
                    y = c(rep(-1, 3), rep(0, 3), rep(1,3)),
                    actifs = rep(c(15, 150, 30), 3))

ze <- data.frame(  x = rep(c(-.5, 0, .5, 1)*50, 3), 
                    y = c(rep(-1, 4), rep(0, 4), rep(1,4)),
                    emplois = rep(c(10, 100, 20), 4))

cas <- full_join(
    loc, ze, by = c("x", "y")
)

points_cas <- cas |> 
    select(x,y) |> 
    as.matrix() |> 
    sf::st_multipoint(dim = "XY") |> 
    sf::st_sfc() |> 
    sf::st_cast(to = "POINT")

dist_cas <- sf::st_distance(points_cas)
dist_cas[dist_cas > 60] <- NA_real_
dist_cas[dist_cas==0] <- .5

dist_cas <- as(as(as(dist_cas, "dMatrix"), "generalMatrix"), "TsparseMatrix")
dist_cas <- data.frame(fromidINS = dist_cas@i, toidINS = dist_cas@j, x = dist_cas@x)


shuf_cas <- matrix(1:9, nrow=1)
cas |> select(emplois, actifs) |> summarise(across(.cols = c(emplois, actifs), .fns = \(x) sum(x, na.rm = TRUE)))
cas <- cas |> 
    mutate(fuite = 0.125,
    gfrom = case_when(
    x < 0 ~ 0L,
    x == 0 ~ 1L,
    x > 0 ~ 2L),
    gto = gfrom)

#dist1 <- as(as(as(dist_cas, "dMatrix"), "generalMatrix"), "TsparseMatrix")
#dist2 <- data.frame(fromidINS = dist1@i, toidINS = dist1@j, x = dist1@x)
#is_triplet(dist2)


res_cont <- meaps_continu(dist = dist_cas, emplois = ze$emplois, actifs = loc$emplois, f = rep(.125, 9), shuf =shuf_cas, 
                        attraction = "constant")





dist0 <- matrix(c(50, 1000, 25000, 26000, 50000,
                  200, 50, 23000, 25000, 49000,
                  20000, 22000, 50, 200, 26000,
                  21000, 21500, 200, 50, 24000,
                  49000, 49000, 25600, 24000, 50), byrow = TRUE, nrow = 5)

emplois0 <- c(100, 1500, 2000, 500, 600)
habitant0 <- c(500, 1600, 1800, 400, 700) 

sum(habitant0)
sum(emplois0)
fuite <- sum(habitant0)/sum(emplois0) - 1
f <- rep(fuite, 5L)
f

gfrom <- c(1L, 1L, 2L, 2L, 3L)
gto <- gfrom

shuf <- matrix(1:5, nrow = 1)

dist1 <- as(as(as(dist0, "dMatrix"), "generalMatrix"), "TsparseMatrix")
dist2 <- data.frame(fromidINS = dist1@i, toidINS = dist1@j, x = dist1@x)
is_triplet(dist2)
dist3 <- prep_meaps_dist(dist2)

res_agg <- meaps_optim (dist_prep = dist3$RankedMat, emplois = emplois0, actifs = habitant0, fuite = f, shuf =shuf, 
                        groups_from = gfrom, groups_to =gto,
                        attraction = "constant")
sum(res_agg)

res_cont <- meaps_continu(dist = dist2, emplois = emplois0, actifs = habitant0, f = f, shuf =shuf, 
                        attraction = "constant")

sparseMatrix(i = res_cont$fromidINS + 1L, j = res_cont$toidINS + 1L, x = res_cont$flux)


