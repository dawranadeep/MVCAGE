memory.size(max = T)
setwd("/Users/rdaw/Documents/MVCAGE/ocean")
library(gstat)
library(sp)





#  First, separate the data in a dataframe
load("data/oceanExample.RData")
df = cbind(ocean@data$CHA, ocean@data$X2)


#  Standardize and complete cases
cc = ocean@coords
mmm = apply(cc, 2, min)
cc = t(t(cc) - mmm)
MMM = apply(cc, 2, max)
cc = t(t(cc) /MMM)

o1 = ocean
idx = which(complete.cases(df))
ddd = data.frame(df[idx, ])
names(ddd) = c("CHA", "ROMS")
o1@data = ddd
o1@coords = cc[idx, ]
proj4string(o1) <- CRS("+init=epsg:28992")


write.csv(ocean@coords[idx, ] , "output/cc2.csv", row.names =    F)



#   Fit the variogram
v.CHA <- variogram(CHA ~ 1, o1, cutoff = 0.35)
v.ROMS <- variogram(ROMS ~ 1, o1, cutoff = 0.35)

f = function(x) attr(m.fit <<- fit.variogram(v.CHA, vgm(,"Mat",nugget=NA,kappa=x)),"SSErr")
opt = optimize(f, c(0.01, 1))
v.CHA.fit = fit.variogram(v.CHA, vgm(,"Mat",nugget=NA,kappa=opt$minimum), fit.kappa = F )
v.CHA.fit
plot(v.CHA, v.CHA.fit)




(v.ROMS.fit <- fit.variogram(v.ROMS, vgm(, "Mat", nugget = NA, kappa=1), fit.kappa = F))
#(v.CHA.fit <- fit.variogram(v.CHA, vgm(, "Mat", nugget = 0.1, range = 0.1, kappa=0.7), fit.kappa = F ))

plot(v.CHA, v.CHA.fit)
plot(v.ROMS, v.ROMS.fit)

# Cross variogram
(g <- gstat(NULL, id = "ROMS", form = ROMS ~ 1, data=o1))
(g <- gstat(g, id = "CHA", form = CHA ~ 1, data=o1))
v.cross <- variogram(g, cutoff = 0.35)
plot(v.cross, cutoff = 0.35)


(g <- gstat(g, id = "ROMS", model = v.ROMS.fit, fill.all=T))
(g <- fit.lmc(v.cross, g, fit.method=6, correct.diagonal=1.0))
plot(v.cross, model=g$model)




v11 <- g$model$ROMS
v22 <- g$model$CHA
v12 <- g$model$ROMS.CHA


#DD = spDists(o1@coords)
#c2 = matrix(runif(20000), ncol = 2)
DD = spDists(o1@coords)


### Estimated Covariance matrix
covMat11 <- variogramLine(v11, dist_vector = DD, covariance = TRUE)
covMat22 <- variogramLine(v22, dist_vector = DD, covariance = TRUE)
covMat12 <- variogramLine(v12, dist_vector = DD, covariance = TRUE)

C11 <- covMat11 - v11$psill[1] * diag(nrow(covMat11))
C22 <- covMat22 - v22$psill[1] * diag(nrow(covMat11))
C12 <- covMat12 - v12$psill[1] * diag(nrow(covMat11))



### Save the covariances
write.csv(C11, "output/C11.csv", row.names =   F, col.names = F)
write.csv(C12, "output/C12.csv", row.names =   F, col.names = F)
write.csv(C22, "output/C22.csv", row.names =   F, col.names = F)
write.csv(o1@coords , "output/cc.csv", row.names =    F, col.names = F)
write.csv(ddd, "output/ocean.csv", row.names = F, col.names = F)


#cc = ocean@coords
knotcc = t(t(oceanKnots) - mmm)
knotcc = t(t(knotcc) / MMM)
write.csv(knotcc, "output/oceanKnots.csv", row.names=F, col.names=F)


  