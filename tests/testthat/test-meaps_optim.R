

M1 <- matrix(1:16 +.1, nrow = 4 )
emp1 <- rep(10., 4)
act1 <- rep(11., 4)
f1 <- (sum(act1) / sum(emp1)) - 1
f1 <- rep(f1,4)
shuf1 <- matrix(c(1:4), nrow = 1)
gfrom1 <- c(1L,1L,2L,2L)
gto1 <- c(1L,1L, 2L, 2L)


meaps_optim(M1, emp1, act1, f = f1, shuf = shuf1,  groups_from = gfrom1, groups_to =gto1)

meaps_optim(M1, emp1, act1, f = f1, shuf = shuf1,  groups_from = gfrom1, groups_to =gto1,
            attraction = "marche", alpha = 5, beta = 10)

meaps_optim(M1, emp1, act1, f = f1, shuf = shuf1,  groups_from = gfrom1, groups_to =gto1,
            attraction = "logistique", alpha = 5, beta = 10)
