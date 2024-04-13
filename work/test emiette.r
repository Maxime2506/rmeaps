library(rmeaps)

actif <- c(p1 = 100, p2 = 1000, p3 = 1, p4 = 1)
shuf <- emiette(actif, 2, seuil = 500)

actif2 <- actif |> sort()
reordonne_shuf(shuf, actif2)
