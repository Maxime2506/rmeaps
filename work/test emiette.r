library(rmeaps)

actif <- c(p1 = 100, p2 = 1000, p3 = 1, p4 = 1, p5 = 1000, p6 = 1000, p7 = 200, p8 = 80, p9 = 100, p10 = 2)
shuf <- emiette(actif, 2, seuil = 500)

actif2 <- actif[-1] |> sort()
reordonne_shuf(shuf, actif2)
