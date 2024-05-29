library(tidyverse)
library(sf)
library(Matrix)
devtools::load_all()


jr = c(0:4,4:0)
p = c(0,5,10)
xr= c(1, 2, 2, 3, 4, 2, 2, 3, 4, 4)
actifs = c(5.5, 11)
emplois = c(3,3,3,3,3)

fuites = c(.1, .1)
shuf = matrix(c(0:1, 1:0), byrow = TRUE, nrow = 2)

newmultishuf(jr, p, xr, emplois, actifs, fuites, parametres = 0, shuf, nthreads = 1)






