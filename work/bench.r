 ## MEAPS avec file d'attente et hétérogénité absorption/fuite

library(tidyverse)
library(conflicted)
source("work/radiation functions.r")
library(tictoc)
library(matrixStats)
library(Rcpp)
library(rmeaps)
conflict_prefer("filter","dplyr")
Rcpp::sourceCpp("work/meaps.cpp")
Rcpp::sourceCpp("work/meaps_rcpp.cpp")

n <- 5000
k <- 4500
bins <- 1.2/0.05
# Rcpp::sourceCpp("radiation/meaps.cpp")
# Rcpp::sourceCpp("radiation/meaps_rcpp.cpp")
# Rcpp::sourceCpp("radiation/meaps_omp.cpp")

hab <- cbind(pos_cnorm(n=n, sigma = 0.2), f=0.1, g=1)
emp <- cbind(pos_cnorm(n=k, sigma = 0.05), p=1, g=1)
s <- make_tibs(emp, hab)
shufs <- do.call(rbind,map(1:50, ~sample.int(n,n)))

(meaps(s$rk,f = s$f, p = s$p, shuf = 1:n, exact=FALSE)$emps) |> colSums() |> summary()
(meaps_cpp(s$rk,f = s$f, p = s$p, shuf = 1:n)$emps) |> colSums() |> summary()
(meaps_cpp_emp(s$rk,f = s$f, p = s$p, shuf = 1:n)) |> colSums() |> summary()
(meaps_rcpp(s$rk,f = s$f, p = s$p, shuf = 1:n)) |> colSums() |> summary()
(meaps_single(rkdist=s$rk, emplois=rep(1, k), actifs = rep(1, n),  f=rep(.1,n), modds=matrix(1, nrow=n, ncol=k), shuf = 1:n ))|> colSums() |> summary()
tic();mm <- rmeaps_bstp(s, shufs, workers=15); toc();
# cpp11::cpp_source("radiation/meaps_cpp11.cpp")
m  <-  matrix(1, nrow=n, ncol=k)
emplois=rep(1, k)
actifs = rep(1, n)
f=rep(.1,n)
shuf = 1:n
bench::mark(r = meaps(s$rk,f = s$f, p = s$p, shuf = 1:n),
            rmeaps = meaps_single(rkdist=s$rk, emplois=emplois, actifs =actifs,  f=f, modds=m, shuf = shuf ),
            # cppe = meaps_cpp_emp(s$rk,f = s$f, p = s$p, shuf = 1:n),
            rcpp = meaps_cpp(s$rk,f = s$f, p = s$p, shuf = 1:n),
            check=FALSE)
m1 <- meaps(s$rk,f = s$f, p = s$p, shuf = 1:n)
m2 <- meaps_cpp(s$rk,f = s$f, p = s$p, shuf = 1:n)$emps
m3 <- meaps_rcpp(s$rk,f = s$f, p = s$p,shuf = 1:n)

## test cpp11
Rcpp::sourceCpp("work/meaps.cpp")
# cpp11::cpp_source("radiation/meaps_cpp11.cpp")

meaps_cpp(rbind(c(1L,2L,3L), c(2L,1L,3L),c(2L,1L,3L)), f=c(0.1, 0.1, 0.1), p=c(1,1,1), shuf = 1:3)
meaps(rbind(c(1L,2L,3L), c(2L,1L,3L),c(2L,1L,3L)), f=c(0.1, 0.1, 0.1), p=c(1,1,1), shuf = 1:3)
meaps_cpp(rbind(c(1L,2L,3L), c(2L,1L,3L),c(2L,1L,3L),c(2L,1L,3L)), f=c(0.1, 0.1, 0.1, 0.1), p=c(1,1,1), shuf = 1:4)
meaps(rbind(c(1L,2L,3L), c(2L,1L,3L),c(2L,1L,3L),c(2L,1L,3L)), f=c(0.1, 0.1, 0.1, 0.1), p=c(1,1,1), shuf = 1:4)
meaps_cpp_corevec(rbind(c(1L,2L,3L), c(2L,1L,3L),c(2L,1L,3L),c(2L,1L,3L)), f=c(0.1, 0.1, 0.1, 0.1), p=c(1,1,1), shuf = 1:4)
meaps_cpp(rbind(c(1L,2L,3L), c(2L,1L,3L),c(2L,1L,3L),c(2L,1L,3L)), f=c(0.1, 0.1, 0.1, 0.1), p=c(1,1,1), shuf = 1:4)
meaps(rbind(c(1L,2L,3L), c(2L,1L,3L),c(2L,1L,3L),c(2L,1L,3L)), f=c(0.1, 0.1, 0.1, 0.1), p=c(1,1,1), shuf = 1:4)

bench::mark()

meaps(rbind(c(1L,2L,3L), c(2L,1L,3L), c(2L,1L,3L)), f=c(0.1, 0.1, 0.1), p=c(1,1,1), shuf = c(1L,2L, 3L))
meaps(rbind(c(1L,2L,3L), c(2L,1L,3L), c(2L,1L,3L), c(2L,1L,3L), c(2L,1L,3L)), f=c(0.1, 0.1, 0.1, 0.1, 0.1), p=c(1,1,1), shuf = 1L:5L)

meaps(rkdist,f = f, p = p, shuf = 1:k)


## ergodicité
library(furrr)
plan("multisession", workers = 8)

coco <- reduce(map(1:100, ~meaps_cpp(rkdist,f = f, p = p, shuf = sample.int(k,k))$occup),
               `+`)/100

fdata <- matrixStats::colSds(coco)/  matrixStats::colMeans2(coco)
