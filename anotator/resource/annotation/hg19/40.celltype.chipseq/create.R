


tab = read.csv("1-s2.0-S2211124715007640-mmc2.csv",header =T, as.is = T)

sel = (tab$mcv == 40)

tab = (tab[sel,c(1:9)])

str(tab)

min(tab[[3]] - tab[[2]])
max(tab[[3]] - tab[[2]])
hist(tab[[3]] - tab[[2]])


stop("")
write.table(file="constitutive.CTCF.txt", tab[sel,c(1:9)], sep = "\t", row.names=F)
