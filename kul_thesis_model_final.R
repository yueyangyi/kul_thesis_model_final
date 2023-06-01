####################################
### KU Leuven Master's Thesis
####################################
### required packages
remove(list=ls())
library(raster)
library(rgdal)
library(INLA)
library(splancs)
library(fields)
library(dplyr)
library(doBy)
library(viridisLite)
library(ape)
library(gstat)
library(readxl)

####################################
### data cleaning
####################################
# import data free to use
ntl_utc_15_30as_grid <- raster("~/netherlands data/harmntl_utc_2015_30arcsec_grid.tif")
slope_utc_3as_grid <- raster("~/netherlands data/slope_utc_3arcsec_grid.tif")
pop_utc_15_30as_grid <- raster("~/netherlands data/ghspop_utc_2015_30arcsec_grid.tif")

# prepare for original pop mapping at 30 arcsec
slope_utc_resampled_grid <- resample(slope_utc_3as_grid, pop_utc_15_30as_grid, method = 'ngb')
ntl_utc_15_resampled_grid <- resample(ntl_utc_15_30as_grid, pop_utc_15_30as_grid, method = 'ngb')

stacked_utc_30as <- stack(pop_utc_15_30as_grid, ntl_utc_15_resampled_grid, slope_utc_resampled_grid)
stacked_utc_30as_df <- as.data.frame(stacked_utc_30as, xy = TRUE)
colnames(stacked_utc_30as_df) <- c("x", "y", "pop_utc_15_30as", "ntl_utc_15_30as", "slope_utc_30as") 
stacked_utc_30as_df[is.na(stacked_utc_30as_df)]<-0 

utc_outline <- readOGR(dsn="~/netherlands data",layer="prov_utc_2015_outline")
crs.geo <- CRS("+proj=longlat +datum=WGS84 +no_defs")
utc_outline  <- spTransform(utc_outline, crs.geo)
utc.ov.30as <- over(SpatialPoints(stacked_utc_30as_df[,1:2], proj4string=crs.geo), utc_outline)

stacked_utc_30as_df <- cbind(stacked_utc_30as_df, utc.ov.30as)
stacked_utc_30as_na <- stacked_utc_30as_df[is.na(stacked_utc_30as_df$statnaam),]# |stacked_utc_30as_df$pop_utc_15_30as<1
stacked_utc_30as_fit <- stacked_utc_30as_df[!is.na(stacked_utc_30as_df$statnaam),]# &stacked_utc_30as_df$pop_utc_15_30as>=1

utc.x.res <- dim(pop_utc_15_30as_grid)[2] 
utc.y.res <- dim(pop_utc_15_30as_grid)[1] 

utc.seq.x.grid <- seq(from=pop_utc_15_30as_grid@extent[1],to=pop_utc_15_30as_grid@extent[2],length=utc.x.res)
utc.seq.y.grid <- seq(from=pop_utc_15_30as_grid@extent[3],to=pop_utc_15_30as_grid@extent[4],length=utc.y.res)

utc_30as_na_sub <- stacked_utc_30as_na[,c(1,2,6)]

stack_utc_30as_real <- stacked_utc_30as_fit[,c(1,2,3)]
stack_utc_30as_real <- as.data.frame(mapply(c,stack_utc_30as_real, utc_30as_na_sub))
stack_utc_30as_real_order <- orderBy(~y + x, data=stack_utc_30as_real)
pop_utc_2015_30arcsec_grid.mat <- log(matrix(as.matrix(as.numeric(stack_utc_30as_real_order$pop_utc_15_30as)+1), 
                                             nrow = utc.x.res))

# set colour
mapping.color.c <- function(n = 201) {
  return(viridis(n))
}

utc_borderline <- readOGR(dsn="~/netherlands data",layer="prov_utc_2015_borderline")

### visualise spatial model and make comparison 
par(mfrow=c(1,1), mar = c(5, 5, 2, 2))
# plot original pop map on sqrt scale
image.plot(utc.seq.x.grid,utc.seq.y.grid,pop_utc_2015_30arcsec_grid.mat,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid,utc.seq.y.grid,pop_utc_2015_30arcsec_grid.mat, add=T, lwd=2, labcex=1)
lines(utc_borderline, lwd=4)

### construct finer grid
utc.seq.x.grid.15as <- seq(from=pop_utc_15_30as_grid@extent[1],to=pop_utc_15_30as_grid@extent[2],length=utc.x.res*2)
utc.seq.y.grid.15as <- seq(from=pop_utc_15_30as_grid@extent[3],to=pop_utc_15_30as_grid@extent[4],length=utc.y.res*2)
utc.pred.grid.15as <- as.matrix(expand.grid(x=utc.seq.x.grid.15as,y=utc.seq.y.grid.15as))

# resample ntl data to 3 arcsec resolution
# empty predict grids as basis
utc.empty.raster.15as<-rasterFromXYZ(cbind(utc.pred.grid.15as,as.vector(rep(NA, nrow(utc.pred.grid.15as)))))

# smoothly disaggregate ntl and resample to predict grid
ntl_utc_resampled_15as <- disaggregate(ntl_utc_15_30as_grid, fact=c(2, 2), method='bilinear')
ntl_utc_resampled_15as <- resample(ntl_utc_resampled_15as, utc.empty.raster.15as, method = 'ngb')

bsa_utc_15_3as_grid <- raster("~/netherlands data/ghsbuilts_utc_2015_100m_grid.tif")

bsa_utc_resampled_15as <- resample(bsa_utc_15_3as_grid, utc.empty.raster.15as, method = 'ngb')
slope_utc_resampled_15as <- resample(slope_utc_3as_grid, utc.empty.raster.15as, method = 'ngb')

pop_utc_resampled_30as <- disaggregate(pop_utc_15_30as_grid, fact=c(2, 2))
pop_utc_resampled_30as <- resample(pop_utc_resampled_30as, utc.empty.raster.15as, method = 'ngb')

pop_utc_15_9as_grid <- raster("~/netherlands data/ghspop_utc_2015_9arcsec_grid.tif")
pop_utc_resampled_9as <- raster::aggregate(pop_utc_15_9as_grid, fact=c(2, 2), fun=sum)
pop_utc_resampled_9as <- resample(pop_utc_resampled_9as, utc.empty.raster.15as, method = 'ngb')

stacked_utc_pred_15as <- stack(ntl_utc_resampled_15as, bsa_utc_resampled_15as, slope_utc_resampled_15as, 
                               pop_utc_resampled_30as, pop_utc_resampled_9as)
stacked_utc_pred_15as_df <- as.data.frame(stacked_utc_pred_15as, xy = TRUE)
colnames(stacked_utc_pred_15as_df) <- c("x", "y", "ntl_utc_15_15as", "bsa_utc_15_15as","slope_utc_15as", 
                                        "pop_utc_15_30as", "pop_utc_15_15as_fake")
stacked_utc_pred_15as_df[is.na(stacked_utc_pred_15as_df)]<-0 

utc.ov.15as <- over(SpatialPoints(stacked_utc_pred_15as_df[,1:2], proj4string=crs.geo), utc_outline)

stacked_utc_pred_15as_df <- cbind(stacked_utc_pred_15as_df, utc.ov.15as)
stacked_utc_pred_15as_na <- stacked_utc_pred_15as_df[is.na(stacked_utc_pred_15as_df$statnaam),] # |stacked_utc_pred_15as_df$pop_utc_15_30as<1
stacked_utc_pred_15as_fit <- stacked_utc_pred_15as_df[!is.na(stacked_utc_pred_15as_df$statnaam),] # &stacked_utc_pred_15as_df$pop_utc_15_30as>=1

utc_cbsstat <- readOGR(dsn="~/netherlands data",layer="cbsstat_utc_2015_500m")
utc_ghssmod <- readOGR(dsn="~/netherlands data",layer="ghssmod_utc_2015_1km")
utc_muniout <- readOGR(dsn="~/netherlands data",layer="muni_utc_2015_outline")

uct.cbsstat.15as <- over(SpatialPoints(stacked_utc_pred_15as_fit[,1:2], proj4string=crs.geo), utc_cbsstat)
uct.ghssmod.15as <- over(SpatialPoints(stacked_utc_pred_15as_fit[,1:2], proj4string=crs.geo), utc_ghssmod)
uct.muniout.15as <- over(SpatialPoints(stacked_utc_pred_15as_fit[,1:2], proj4string=crs.geo), utc_muniout)

stacked_utc_pred_15as_fit <- cbind(stacked_utc_pred_15as_fit, uct.cbsstat.15as, uct.ghssmod.15as, uct.muniout.15as)
colnames(stacked_utc_pred_15as_fit)[8:13] <- c("prov_utc_15", "pop_utc_15_15as_real", "hhs_utc_15_15as", 
                                               "psd_utc_15_15as", "stt_utc_15_30as", "muni_utc_15")
stacked_utc_pred_15as_fit[is.na(stacked_utc_pred_15as_fit)]<-0 
stacked_utc_pred_15as_fit[stacked_utc_pred_15as_fit==-99997]<-0

####################################
### exploratory data analysis 1: assess assumptions
####################################
# eda for data at 30 arcsec
utc.loc_est_30as <- cbind(stacked_utc_30as_fit$x, stacked_utc_30as_fit$y)
utc.30as.sp <- SpatialPointsDataFrame(SpatialPoints(utc.loc_est_30as), 
                                      data = data.frame(counts = log(stacked_utc_30as_fit$pop_utc_15_30as+1)))
utc.30as.g <- gstat(id = "counts", formula = counts ~ 1, data = utc.30as.sp)
utc.30as.evgm <- variogram(utc.30as.g, cutoff = 0.2, width = 0.003)
utc.30as.revgm <- variogram(utc.30as.g, cutoff = 0.2, width = 0.003, cressie = TRUE)
utc.30as.aevgm <- variogram(utc.30as.g, cutoff = 0.2, width = 0.003, alpha = c(0, 45, 90, 135))

# eda for fake data at 15 arcsec
utc.loc_pred_15as <- cbind(stacked_utc_pred_15as_fit$x, stacked_utc_pred_15as_fit$y)
utc.15as.sp.f <- 
  SpatialPointsDataFrame(SpatialPoints(utc.loc_pred_15as), 
                         data = data.frame(counts = log(stacked_utc_pred_15as_fit$pop_utc_15_15as_fake+1)))
utc.15as.g.f <- gstat(id = "counts", formula = counts ~ 1, data = utc.15as.sp.f)
utc.15as.evgm.f <- variogram(utc.15as.g.f, cutoff = 0.2, width = 0.003)
utc.15as.revgm.f <- variogram(utc.15as.g.f, cutoff = 0.2, width = 0.003, cressie = TRUE)
utc.15as.aevgm.f <- variogram(utc.15as.g.f, cutoff = 0.2, width = 0.003, alpha = c(0, 45, 90, 135))

# eda for real data at 15 arcsec
utc.15as.sp.r <- 
  SpatialPointsDataFrame(SpatialPoints(utc.loc_pred_15as), 
                         data = data.frame(counts = log(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as_real)+1)))
utc.15as.g.r <- gstat(id = "counts", formula = counts ~ 1, data = utc.15as.sp.r)
utc.15as.evgm.r <- variogram(utc.15as.g.r, cutoff = 0.2, width = 0.003)
utc.15as.revgm.r <- variogram(utc.15as.g.r, cutoff = 0.2, width = 0.003, cressie = TRUE)
utc.15as.aevgm.r <- variogram(utc.15as.g.r, cutoff = 0.2, width = 0.003, alpha = c(0, 45, 90, 135))

# plot semivariogram for utc
# check assumption of second-order stationary
par(mfrow=c(1,1), mar = c(5, 5, 2, 2))
plot(cbind(utc.30as.evgm$dist, utc.30as.evgm$gamma), 
     xlim=c(0,0.2),ylim=c(0,7.5),col="blue",t="o",ylab=expression(hat(gamma)), xlab="Distance (degree)", lty=1, pch=20)
lines(cbind(utc.30as.revgm$dist, utc.30as.revgm$gamma),t="o",xlab="",ylab="",col="blue", lty=3, pch=20)
lines(cbind(utc.15as.evgm.f$dist, utc.15as.evgm.f$gamma),t="o",xlab="",ylab="",col="green", lty=1, pch=20)
lines(cbind(utc.15as.revgm.f$dist, utc.15as.revgm.f$gamma),t="o",xlab="",ylab="",col="green", lty=3, pch=20)
lines(cbind(utc.15as.evgm.r$dist, utc.15as.evgm.r$gamma),t="o",xlab="",ylab="",col="red", lty=1, pch=20)
lines(cbind(utc.15as.revgm.r$dist, utc.15as.revgm.r$gamma),t="o",xlab="",ylab="",col="red", lty=3, pch=20)
legend(x = "bottomright", legend=c("GHS-POP at 30 arcsec, classic", "GHS-POP at 30 arcsec, robust",
                                   "GHS-POP at 15 arcsec, classic", "GHS-POP at 15 arcsec, robust", 
                                   "CBS-POP at 15 arcsec, classic", "CBS-POP at 15 arcsec, robust"), 
       col=c("blue", "blue", "green", "green", "red", "red"), lty=c(1,3,1,3,1,3), pch=c(20,20,20,20,20,20), 
       lwd=c(1,1,1,1,1,1), cex=0.8)

# check assumption of isotropy
utc.30as.aevgm.0 <- subset(utc.30as.aevgm, dir.hor==0, select=c(dist,gamma))
utc.30as.aevgm.45 <- subset(utc.30as.aevgm, dir.hor==45, select=c(dist,gamma))
utc.30as.aevgm.90 <- subset(utc.30as.aevgm, dir.hor==90, select=c(dist,gamma))
utc.30as.aevgm.135 <- subset(utc.30as.aevgm, dir.hor==135, select=c(dist,gamma))

utc.15as.aevgm.0.f <- subset(utc.15as.aevgm.f, dir.hor==0, select=c(dist,gamma))
utc.15as.aevgm.45.f <- subset(utc.15as.aevgm.f, dir.hor==45, select=c(dist,gamma))
utc.15as.aevgm.90.f <- subset(utc.15as.aevgm.f, dir.hor==90, select=c(dist,gamma))
utc.15as.aevgm.135.f <- subset(utc.15as.aevgm.f, dir.hor==135, select=c(dist,gamma))

utc.15as.aevgm.0.r <- subset(utc.15as.aevgm.r, dir.hor==0, select=c(dist,gamma))
utc.15as.aevgm.45.r <- subset(utc.15as.aevgm.r, dir.hor==45, select=c(dist,gamma))
utc.15as.aevgm.90.r <- subset(utc.15as.aevgm.r, dir.hor==90, select=c(dist,gamma))
utc.15as.aevgm.135.r <- subset(utc.15as.aevgm.r, dir.hor==135, select=c(dist,gamma))

par(mfrow=c(1,1), mar = c(5, 5, 2, 2))
plot(cbind(utc.30as.aevgm.0$dist, utc.30as.aevgm.0$gamma), 
     xlim=c(0,0.2),ylim=c(0,8),col="blue",t="o",ylab=expression(hat(gamma)), xlab="Distance (degree)", lty=1, pch=15)
lines(cbind(utc.30as.aevgm.45$dist, utc.30as.aevgm.45$gamma),t="o",xlab="",ylab="",col="blue", lty=1, pch=16)
lines(cbind(utc.30as.aevgm.90$dist, utc.30as.aevgm.90$gamma),t="o",xlab="",ylab="",col="blue", lty=1, pch=17)
lines(cbind(utc.30as.aevgm.135$dist, utc.30as.aevgm.135$gamma),t="o",xlab="",ylab="",col="blue", lty=1, pch=18)
lines(cbind(utc.15as.aevgm.0.f$dist, utc.15as.aevgm.0.f$gamma),t="o",xlab="",ylab="",col="green", lty=1, pch=15)
lines(cbind(utc.15as.aevgm.45.f$dist, utc.15as.aevgm.45.f$gamma),t="o",xlab="",ylab="",col="green", lty=1, pch=16)
lines(cbind(utc.15as.aevgm.90.f$dist, utc.15as.aevgm.90.f$gamma),t="o",xlab="",ylab="",col="green", lty=1, pch=17)
lines(cbind(utc.15as.aevgm.135.f$dist, utc.15as.aevgm.135.f$gamma),t="o",xlab="",ylab="",col="green", lty=1, pch=18)
lines(cbind(utc.15as.aevgm.0.r$dist, utc.15as.aevgm.0.r$gamma),t="o",xlab="",ylab="",col="red", lty=1, pch=15)
lines(cbind(utc.15as.aevgm.45.r$dist, utc.15as.aevgm.45.r$gamma),t="o",xlab="",ylab="",col="red", lty=1, pch=16)
lines(cbind(utc.15as.aevgm.90.r$dist, utc.15as.aevgm.90.r$gamma),t="o",xlab="",ylab="",col="red", lty=1, pch=17)
lines(cbind(utc.15as.aevgm.135.r$dist, utc.15as.aevgm.135.r$gamma),t="o",xlab="",ylab="",col="red", lty=1, pch=18)
legend(x = "bottomright", legend=c("GHS-POP at 30 arcsec, 0°", "GHS-POP at 30 arcsec, 45°", 
                                   "GHS-POP at 30 arcsec, 90°","GHS-POP at 30 arcsec, 135°",
                                   "GHS-POP at 15 arcsec, 0°", "GHS-POP at 15 arcsec, 45°", 
                                   "GHS-POP at 15 arcsec, 90°","GHS-POP at 15 arcsec, 135°",
                                   "CBS-POP at 15 arcsec, 0°", "CBS-POP at 15 arcsec, 45°", 
                                   "CBS-POP at 15 arcsec, 90°","CBS-POP at 15 arcsec, 135°"), 
       col=c("blue", "blue", "blue", "blue", "green", "green", "green", "green", "red", "red", "red", "red"), 
       lty=c(1,1,1,1,1,1,1,1,1,1,1,1), pch=c(15,16,17,18,15,16,17,18,15,16,17,18), 
       lwd=c(1,1,1,1,1,1,1,1,1,1,1,1), cex=0.8)


####################################
### exploratory data analysis 2: assess ability of approximating real data
####################################
# full data
res.15as.full <- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as_real) - stacked_utc_pred_15as_fit$pop_utc_15_15as_fake
RMSE.15as.full <- sqrt(mean(res.15as.full^2))
RMSE.15as.full
correl.15as.full <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as_real), stacked_utc_pred_15as_fit$pop_utc_15_15as_fake) 
correl.15as.full 

# unweighted sampling
stacked_utc_pred_15as_fit[,ncol(stacked_utc_pred_15as_fit)+1] <- 1:dim(stacked_utc_pred_15as_fit)[1]
colnames(stacked_utc_pred_15as_fit)[ncol(stacked_utc_pred_15as_fit)] <- "ind4samp"

# unweighted 1% data
r_1 <- 0.01
set.seed(666)
utc_rnum_u1 <- sample(1:dim(stacked_utc_pred_15as_fit)[1], r_1*dim(stacked_utc_pred_15as_fit)[1], replace=F)
stacked_utc_sampled_u1 <- stacked_utc_pred_15as_fit[stacked_utc_pred_15as_fit$ind4samp %in% utc_rnum_u1,]

res.15as.u1 <- as.numeric(stacked_utc_sampled_u1$pop_utc_15_15as_real) - stacked_utc_sampled_u1$pop_utc_15_15as_fake
RMSE.15as.u1 <- sqrt(mean(res.15as.u1^2))
RMSE.15as.u1
correl.15as.u1 <- cor(as.numeric(stacked_utc_sampled_u1$pop_utc_15_15as_real), stacked_utc_sampled_u1$pop_utc_15_15as_fake) 
correl.15as.u1 

# unweighted 2% data
r_2 <- 0.02
set.seed(666)
utc_rnum_u2 <- sample(1:dim(stacked_utc_pred_15as_fit)[1], r_2*dim(stacked_utc_pred_15as_fit)[1], replace=F)
stacked_utc_sampled_u2 <- stacked_utc_pred_15as_fit[stacked_utc_pred_15as_fit$ind4samp %in% utc_rnum_u2,]

res.15as.u2 <- as.numeric(stacked_utc_sampled_u2$pop_utc_15_15as_real) - stacked_utc_sampled_u2$pop_utc_15_15as_fake
RMSE.15as.u2 <- sqrt(mean(res.15as.u2^2))
RMSE.15as.u2
correl.15as.u2 <- cor(as.numeric(stacked_utc_sampled_u2$pop_utc_15_15as_real), stacked_utc_sampled_u2$pop_utc_15_15as_fake) 
correl.15as.u2 

# unweighted 5% data
r_3 <- 0.05
set.seed(666)
utc_rnum_u3 <- sample(1:dim(stacked_utc_pred_15as_fit)[1], r_3*dim(stacked_utc_pred_15as_fit)[1], replace=F)
stacked_utc_sampled_u3 <- stacked_utc_pred_15as_fit[stacked_utc_pred_15as_fit$ind4samp %in% utc_rnum_u3,]

res.15as.u3 <- as.numeric(stacked_utc_sampled_u3$pop_utc_15_15as_real) - stacked_utc_sampled_u3$pop_utc_15_15as_fake
RMSE.15as.u3 <- sqrt(mean(res.15as.u3^2))
RMSE.15as.u3
correl.15as.u3 <- cor(as.numeric(stacked_utc_sampled_u3$pop_utc_15_15as_real), stacked_utc_sampled_u3$pop_utc_15_15as_fake) 
correl.15as.u3 

# unweighted 10% data
r_4 <- 0.1
set.seed(666)
utc_rnum_u4 <- sample(1:dim(stacked_utc_pred_15as_fit)[1], r_4*dim(stacked_utc_pred_15as_fit)[1], replace=F)
stacked_utc_sampled_u4 <- stacked_utc_pred_15as_fit[stacked_utc_pred_15as_fit$ind4samp %in% utc_rnum_u4,]

res.15as.u4 <- as.numeric(stacked_utc_sampled_u4$pop_utc_15_15as_real) - stacked_utc_sampled_u4$pop_utc_15_15as_fake
RMSE.15as.u4 <- sqrt(mean(res.15as.u4^2))
RMSE.15as.u4
correl.15as.u4 <- cor(as.numeric(stacked_utc_sampled_u4$pop_utc_15_15as_real), stacked_utc_sampled_u4$pop_utc_15_15as_fake) 
correl.15as.u4 

# unweighted 20% data
r_5 <- 0.2
set.seed(666)
utc_rnum_u5 <- sample(1:dim(stacked_utc_pred_15as_fit)[1], r_5*dim(stacked_utc_pred_15as_fit)[1], replace=F)
stacked_utc_sampled_u5 <- stacked_utc_pred_15as_fit[stacked_utc_pred_15as_fit$ind4samp %in% utc_rnum_u5,]

res.15as.u5 <- as.numeric(stacked_utc_sampled_u5$pop_utc_15_15as_real) - stacked_utc_sampled_u5$pop_utc_15_15as_fake
RMSE.15as.u5 <- sqrt(mean(res.15as.u5^2))
RMSE.15as.u5
correl.15as.u5 <- cor(as.numeric(stacked_utc_sampled_u5$pop_utc_15_15as_real), stacked_utc_sampled_u5$pop_utc_15_15as_fake) 
correl.15as.u5 

# unweighted 50% data
r_6 <- 0.5
set.seed(666)
utc_rnum_u6 <- sample(1:dim(stacked_utc_pred_15as_fit)[1], r_6*dim(stacked_utc_pred_15as_fit)[1], replace=F)
stacked_utc_sampled_u6 <- stacked_utc_pred_15as_fit[stacked_utc_pred_15as_fit$ind4samp %in% utc_rnum_u6,]

res.15as.u6 <- as.numeric(stacked_utc_sampled_u6$pop_utc_15_15as_real) - stacked_utc_sampled_u6$pop_utc_15_15as_fake
RMSE.15as.u6 <- sqrt(mean(res.15as.u6^2))
RMSE.15as.u6
correl.15as.u6 <- cor(as.numeric(stacked_utc_sampled_u6$pop_utc_15_15as_real), stacked_utc_sampled_u6$pop_utc_15_15as_fake) 
correl.15as.u6 

# weighted sampling with weights coming from ghs-pop data at 30 arcsec
# weighted 1% data
set.seed(999)
utc_rnum_w1 <- sample(1:dim(stacked_utc_pred_15as_fit)[1], r_1*dim(stacked_utc_pred_15as_fit)[1], replace=F, 
                      prob=stacked_utc_pred_15as_fit$pop_utc_15_30as)
stacked_utc_sampled_w1 <- stacked_utc_pred_15as_fit[stacked_utc_pred_15as_fit$ind4samp %in% utc_rnum_w1,]

res.15as.w1 <- as.numeric(stacked_utc_sampled_w1$pop_utc_15_15as_real) - stacked_utc_sampled_w1$pop_utc_15_15as_fake
RMSE.15as.w1 <- sqrt(mean(res.15as.w1^2))
RMSE.15as.w1
correl.15as.w1 <- cor(as.numeric(stacked_utc_sampled_w1$pop_utc_15_15as_real), stacked_utc_sampled_w1$pop_utc_15_15as_fake) 
correl.15as.w1 

# weighted 2% data
set.seed(999)
utc_rnum_w2 <- sample(1:dim(stacked_utc_pred_15as_fit)[1], r_2*dim(stacked_utc_pred_15as_fit)[1], replace=F, 
                      prob=stacked_utc_pred_15as_fit$pop_utc_15_30as)
stacked_utc_sampled_w2 <- stacked_utc_pred_15as_fit[stacked_utc_pred_15as_fit$ind4samp %in% utc_rnum_w2,]

res.15as.w2 <- as.numeric(stacked_utc_sampled_w2$pop_utc_15_15as_real) - stacked_utc_sampled_w2$pop_utc_15_15as_fake
RMSE.15as.w2 <- sqrt(mean(res.15as.w2^2))
RMSE.15as.w2
correl.15as.w2 <- cor(as.numeric(stacked_utc_sampled_w2$pop_utc_15_15as_real), stacked_utc_sampled_w2$pop_utc_15_15as_fake) 
correl.15as.w2 

# weighted 5% data
set.seed(999)
utc_rnum_w3 <- sample(1:dim(stacked_utc_pred_15as_fit)[1], r_3*dim(stacked_utc_pred_15as_fit)[1], replace=F, 
                      prob=stacked_utc_pred_15as_fit$pop_utc_15_30as)
stacked_utc_sampled_w3 <- stacked_utc_pred_15as_fit[stacked_utc_pred_15as_fit$ind4samp %in% utc_rnum_w3,]

res.15as.w3 <- as.numeric(stacked_utc_sampled_w3$pop_utc_15_15as_real) - stacked_utc_sampled_w3$pop_utc_15_15as_fake
RMSE.15as.w3 <- sqrt(mean(res.15as.w3^2))
RMSE.15as.w3
correl.15as.w3 <- cor(as.numeric(stacked_utc_sampled_w3$pop_utc_15_15as_real), stacked_utc_sampled_w3$pop_utc_15_15as_fake) 
correl.15as.w3 

# weighted 10% data
set.seed(999)
utc_rnum_w4 <- sample(1:dim(stacked_utc_pred_15as_fit)[1], r_4*dim(stacked_utc_pred_15as_fit)[1], replace=F, 
                      prob=stacked_utc_pred_15as_fit$pop_utc_15_30as)
stacked_utc_sampled_w4 <- stacked_utc_pred_15as_fit[stacked_utc_pred_15as_fit$ind4samp %in% utc_rnum_w4,]

res.15as.w4 <- as.numeric(stacked_utc_sampled_w4$pop_utc_15_15as_real) - stacked_utc_sampled_w4$pop_utc_15_15as_fake
RMSE.15as.w4 <- sqrt(mean(res.15as.w4^2))
RMSE.15as.w4
correl.15as.w4 <- cor(as.numeric(stacked_utc_sampled_w4$pop_utc_15_15as_real), stacked_utc_sampled_w4$pop_utc_15_15as_fake) 
correl.15as.w4 

# weighted 20% data
set.seed(999)
utc_rnum_w5 <- sample(1:dim(stacked_utc_pred_15as_fit)[1], r_5*dim(stacked_utc_pred_15as_fit)[1], replace=F, 
                      prob=stacked_utc_pred_15as_fit$pop_utc_15_30as)
stacked_utc_sampled_w5 <- stacked_utc_pred_15as_fit[stacked_utc_pred_15as_fit$ind4samp %in% utc_rnum_w5,]

res.15as.w5 <- as.numeric(stacked_utc_sampled_w5$pop_utc_15_15as_real) - stacked_utc_sampled_w5$pop_utc_15_15as_fake
RMSE.15as.w5 <- sqrt(mean(res.15as.w5^2))
RMSE.15as.w5
correl.15as.w5 <- cor(as.numeric(stacked_utc_sampled_w5$pop_utc_15_15as_real), stacked_utc_sampled_w5$pop_utc_15_15as_fake) 
correl.15as.w5 

# weighted 50% data
set.seed(999)
utc_rnum_w6 <- sample(1:dim(stacked_utc_pred_15as_fit)[1], r_6*dim(stacked_utc_pred_15as_fit)[1], replace=F, 
                      prob=stacked_utc_pred_15as_fit$pop_utc_15_30as)
stacked_utc_sampled_w6 <- stacked_utc_pred_15as_fit[stacked_utc_pred_15as_fit$ind4samp %in% utc_rnum_w6,]

res.15as.w6 <- as.numeric(stacked_utc_sampled_w6$pop_utc_15_15as_real) - stacked_utc_sampled_w6$pop_utc_15_15as_fake
RMSE.15as.w6 <- sqrt(mean(res.15as.w6^2))
RMSE.15as.w6
correl.15as.w6 <- cor(as.numeric(stacked_utc_sampled_w6$pop_utc_15_15as_real), stacked_utc_sampled_w6$pop_utc_15_15as_fake) 
correl.15as.w6 

### import border
utc_border <- readOGR(dsn="~/netherlands data",layer="prov_utc_2015_border")
utc_border_df <- as.data.frame(utc_border, xy = TRUE)
utc_border_df <- subset(utc_border_df, select = c(coords.x1, coords.x2))

# plot where to sample
col_u <- rgb(red=0.2, green=0.2, blue=1.0, alpha=0.2)
col_w <- rgb(red=1.0, green=0.2, blue=0.2, alpha=0.2)
par(mfrow=c(1,1), mar = c(1, 1, 1, 1))
plot(utc_borderline,lwd=4,col="black")
# points(utc.loc_pred_15as, pch=19, col="black", cex=0.05)
points(cbind(stacked_utc_sampled_u6$x, stacked_utc_sampled_u6$y), pch=19, col=col_u, cex=0.25)
points(cbind(stacked_utc_sampled_u5$x, stacked_utc_sampled_u5$y), pch=19, col=col_u, cex=0.5)
points(cbind(stacked_utc_sampled_u4$x, stacked_utc_sampled_u4$y), pch=19, col=col_u, cex=1)
points(cbind(stacked_utc_sampled_u3$x, stacked_utc_sampled_u3$y), pch=19, col=col_u, cex=2)
points(cbind(stacked_utc_sampled_u2$x, stacked_utc_sampled_u2$y), pch=19, col=col_u, cex=4)
points(cbind(stacked_utc_sampled_u1$x, stacked_utc_sampled_u1$y), pch=19, col=col_u, cex=8)
points(cbind(stacked_utc_sampled_w6$x, stacked_utc_sampled_w6$y), pch=19, col=col_w, cex=0.25)
points(cbind(stacked_utc_sampled_w5$x, stacked_utc_sampled_w5$y), pch=19, col=col_w, cex=0.5)
points(cbind(stacked_utc_sampled_w4$x, stacked_utc_sampled_w4$y), pch=19, col=col_w, cex=1)
points(cbind(stacked_utc_sampled_w3$x, stacked_utc_sampled_w3$y), pch=19, col=col_w, cex=2)
points(cbind(stacked_utc_sampled_w2$x, stacked_utc_sampled_w2$y), pch=19, col=col_w, cex=4)
points(cbind(stacked_utc_sampled_w1$x, stacked_utc_sampled_w1$y), pch=19, col=col_w, cex=8)
legend(x = "topright", legend=c("unweighted sample","weighted sample"), 
       col=c(col_u,col_w),pch=c(19,19),
       cex=0.8)

####################################
### set up inputs and non-informative priors
####################################
# set up non-informative pirors
### make mesh
utc.bound <- inla.nonconvex.hull(utc.loc_est_30as, concave = 0.1, resolution = c(20, 20))
utc.mesh <- inla.mesh.2d(boundary=utc.bound, max.edge=c(0.012, 0.03), cutoff = 0.01) 

utc.spde <- inla.spde2.matern(mesh=utc.mesh, alpha=2, 
                              B.tau=matrix(c(0, 1, 0),nrow=1,ncol=3),
                              B.kappa=matrix(c(0, 0, 1),nrow=1,ncol=3), 
                              prior.variance.nominal = 1,
                              theta.prior.prec=0.1)

utc.s.index <- inla.spde.make.index(name="utc.spatial.field", n.spde=utc.spde$n.spde) 

utc.fixed.vague <- list(mean.intercept = 0, 
                        prec.intercept = 0, 
                        mean = 0,
                        prec = 0.001)

utc.prec.vague <- list(initial=4, prior="loggamma", 
                       param=c(1, 0.00005), fixed=FALSE)

utc.A.est.30as <- inla.spde.make.A(mesh = utc.mesh, loc= utc.loc_est_30as)

utc.A.pred.15as <- inla.spde.make.A(mesh=utc.mesh, loc=utc.loc_pred_15as)

par(mfrow=c(1,1), mar = c(1, 1, 1, 1))
plot(utc.mesh,main="",asp=1)
points(utc.loc_pred_15as, pch=19, col=4, cex=0.05)
points(utc.loc_est_30as, pch=19, col=2, cex=0.3) 
lines(utc_borderline,lwd=4,col=1)
legend(x = "topright", legend=c("location for estimation","location for prediction"), 
       col=c(2,4),pch=c(19,19),
       cex=0.8)

# set up inputs
stacked_utc_pred_15as_na <- stacked_utc_pred_15as_na[,-7]

### CAUTION!!! for fake data case only
stacked_utc_pred_15as_fit <- stacked_utc_pred_15as_fit[,-9]
colnames(stacked_utc_pred_15as_fit)[7] <- "pop_utc_15_15as"

stacked_utc_sampled_u1 <- stacked_utc_sampled_u1[,-9]
colnames(stacked_utc_sampled_u1)[7] <- "pop_utc_15_15as"
stacked_utc_sampled_u2 <- stacked_utc_sampled_u2[,-9]
colnames(stacked_utc_sampled_u2)[7] <- "pop_utc_15_15as"
stacked_utc_sampled_u3 <- stacked_utc_sampled_u3[,-9]
colnames(stacked_utc_sampled_u3)[7] <- "pop_utc_15_15as"
stacked_utc_sampled_u4 <- stacked_utc_sampled_u4[,-9]
colnames(stacked_utc_sampled_u4)[7] <- "pop_utc_15_15as"
stacked_utc_sampled_u5 <- stacked_utc_sampled_u5[,-9]
colnames(stacked_utc_sampled_u5)[7] <- "pop_utc_15_15as"
stacked_utc_sampled_u6 <- stacked_utc_sampled_u6[,-9]
colnames(stacked_utc_sampled_u6)[7] <- "pop_utc_15_15as"

stacked_utc_sampled_w1 <- stacked_utc_sampled_w1[,-9]
colnames(stacked_utc_sampled_w1)[7] <- "pop_utc_15_15as"
stacked_utc_sampled_w2 <- stacked_utc_sampled_w2[,-9]
colnames(stacked_utc_sampled_w2)[7] <- "pop_utc_15_15as"
stacked_utc_sampled_w3 <- stacked_utc_sampled_w3[,-9]
colnames(stacked_utc_sampled_w3)[7] <- "pop_utc_15_15as"
stacked_utc_sampled_w4 <- stacked_utc_sampled_w4[,-9]
colnames(stacked_utc_sampled_w4)[7] <- "pop_utc_15_15as"
stacked_utc_sampled_w5 <- stacked_utc_sampled_w5[,-9]
colnames(stacked_utc_sampled_w5)[7] <- "pop_utc_15_15as"
stacked_utc_sampled_w6 <- stacked_utc_sampled_w6[,-9]
colnames(stacked_utc_sampled_w6)[7] <- "pop_utc_15_15as"

### CAUTION!!! for real data case only
stacked_utc_pred_15as_fit <- stacked_utc_pred_15as_fit[,-7]
colnames(stacked_utc_pred_15as_fit)[8] <- "pop_utc_15_15as"

stacked_utc_sampled_u1 <- stacked_utc_sampled_u1[,-7]
colnames(stacked_utc_sampled_u1)[8] <- "pop_utc_15_15as"
stacked_utc_sampled_u2 <- stacked_utc_sampled_u2[,-7]
colnames(stacked_utc_sampled_u2)[8] <- "pop_utc_15_15as"
stacked_utc_sampled_u3 <- stacked_utc_sampled_u3[,-7]
colnames(stacked_utc_sampled_u3)[8] <- "pop_utc_15_15as"
stacked_utc_sampled_u4 <- stacked_utc_sampled_u4[,-7]
colnames(stacked_utc_sampled_u4)[8] <- "pop_utc_15_15as"
stacked_utc_sampled_u5 <- stacked_utc_sampled_u5[,-7]
colnames(stacked_utc_sampled_u5)[8] <- "pop_utc_15_15as"
stacked_utc_sampled_u6 <- stacked_utc_sampled_u6[,-7]
colnames(stacked_utc_sampled_u6)[8] <- "pop_utc_15_15as"

stacked_utc_sampled_w1 <- stacked_utc_sampled_w1[,-7]
colnames(stacked_utc_sampled_w1)[8] <- "pop_utc_15_15as"
stacked_utc_sampled_w2 <- stacked_utc_sampled_w2[,-7]
colnames(stacked_utc_sampled_w2)[8] <- "pop_utc_15_15as"
stacked_utc_sampled_w3 <- stacked_utc_sampled_w3[,-7]
colnames(stacked_utc_sampled_w3)[8] <- "pop_utc_15_15as"
stacked_utc_sampled_w4 <- stacked_utc_sampled_w4[,-7]
colnames(stacked_utc_sampled_w4)[8] <- "pop_utc_15_15as"
stacked_utc_sampled_w5 <- stacked_utc_sampled_w5[,-7]
colnames(stacked_utc_sampled_w5)[8] <- "pop_utc_15_15as"
stacked_utc_sampled_w6 <- stacked_utc_sampled_w6[,-7]
colnames(stacked_utc_sampled_w6)[8] <- "pop_utc_15_15as"


####################################
# worst case scenario
# full model w/o transformation + vague prior
# est at 30 arcsec full data + pred at 3 arcsec
####################################
utc.stack.est.30as <- inla.stack(data = list(pop=as.numeric(stacked_utc_30as_fit$pop_utc_15_30as)),
                                 A = list(utc.A.est.30as, 1),
                                 effects = list(c(utc.s.index, list(Intercept=1)), 
                                                list(ntl=as.numeric(stacked_utc_30as_fit$ntl_utc_15_30as), 
                                                     slope=as.numeric(stacked_utc_30as_fit$slope_utc_30as)
                                                )),
                                 tag="est")

utc.stack.pred.15as <- inla.stack(data = list(pop=NA),
                                  A = list(utc.A.pred.15as, 1),
                                  effects = list(c(utc.s.index, list(Intercept=1)),
                                                 list(ntl=as.numeric(stacked_utc_pred_15as_fit$ntl_utc_15_15as), 
                                                      slope=as.numeric(stacked_utc_pred_15as_fit$slope_utc_15as)                                                     
                                                 )),
                                  tag="pred.response") 

utc.stack.latent <- inla.stack(data = list(xi=NA),
                               A = list(utc.A.pred.15as),
                               effects = list(utc.s.index),
                               tag="pred.latent")

utc.join.stack.pred.15as <- inla.stack(utc.stack.est.30as, utc.stack.latent, utc.stack.pred.15as)

utc.formula <- pop ~ -1 + Intercept + ntl + slope + f(utc.spatial.field, model=utc.spde)

utc.output.pred.15as <- inla(utc.formula, 
                             data=inla.stack.data(utc.join.stack.pred.15as, spde=utc.spde), family="gaussian",
                             control.fixed=utc.fixed.vague, 
                             control.family=list(hyper=list(prec=utc.prec.vague)), 
                             control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as), compute=TRUE),
                             control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as)

utc.output.pred.15as$dic$dic
utc.output.pred.15as$waic$waic

# compute negative sum of log cpo
slcpo <- function(output, na.rm=TRUE){
  -sum(log(output$cpo$cpo), na.rm = na.rm)
}
slcpo(utc.output.pred.15as)

### extract parameter estimates
utc.fixed.out <- round(utc.output.pred.15as$summary.fixed[,1:5],3) 
utc.fixed1_marg <- utc.output.pred.15as$marginals.fixed[[1]]
utc.fixed2_marg <- utc.output.pred.15as$marginals.fixed[[2]]
utc.fixed3_marg <- utc.output.pred.15as$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as$marginals.hyperpar[[1]])
utc.sigma2e_m1 <- inla.emarginal(function(x) x, utc.sigma2e_marg)
utc.sigma2e_m2 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg)
utc.sigma2e_stdev <- sqrt(utc.sigma2e_m2 - utc.sigma2e_m1^2)
utc.sigma2e_quantiles <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg)

utc.spde.result <- inla.spde2.result(utc.output.pred.15as, name="utc.spatial.field", utc.spde)

# variance of spatial effects
utc.var.nom.marg <- utc.spde.result$marginals.variance.nominal[[1]]
utc.var.nom.m1 <- inla.emarginal(function(x) x, utc.var.nom.marg)
utc.var.nom.m2 <- inla.emarginal(function(x) x^2, utc.var.nom.marg)
utc.var.nom.stdev <- sqrt(utc.var.nom.m2 - utc.var.nom.m1^2)
utc.var.nom.quantiles <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg)

utc.range.nom.marg <- utc.spde.result$marginals.range.nominal[[1]]
utc.range.nom.m1 <- inla.emarginal(function(x) x, utc.range.nom.marg)
utc.range.nom.m2 <- inla.emarginal(function(x) x^2, utc.range.nom.marg)
utc.range.nom.stdev <- sqrt(utc.range.nom.m2 - utc.range.nom.m1^2)
utc.range.nom.quantiles <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg)

### extract predicted values
utc.index.pred.response <- inla.stack.index(utc.join.stack.pred.15as, "pred.response")$data 
utc.index.pred.latent <- inla.stack.index(utc.join.stack.pred.15as, "pred.latent")$data 

utc.pred.response.mean <- utc.output.pred.15as$summary.linear.predictor[utc.index.pred.response, "mean"]
utc.pred.response.sd <- utc.output.pred.15as$summary.linear.predictor[utc.index.pred.response, "sd"]

utc.pred.latent.mean <- utc.output.pred.15as$summary.linear.predictor[utc.index.pred.latent, "mean"]
utc.pred.latent.sd <- utc.output.pred.15as$summary.linear.predictor[utc.index.pred.latent, "sd"]

utc_pred_na_15as_sub <- stacked_utc_pred_15as_na[,c(1,2,7)]

utc_pred_resp_mean <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.mean)
utc_pred_resp_mean_df <- as.data.frame(mapply(c,utc_pred_resp_mean, utc_pred_na_15as_sub))
utc_pred_resp_mean_order <- orderBy(~y + x, data=utc_pred_resp_mean_df)
utc.pred.response.mean.grid <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order$utc.pred.response.mean)), nrow=utc.x.res*2)

utc_pred_resp_sd <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd)
utc_pred_resp_sd_df <- as.data.frame(mapply(c,utc_pred_resp_sd, utc_pred_na_15as_sub))
utc_pred_resp_sd_order <- orderBy(~y + x, data=utc_pred_resp_sd_df)
utc.pred.response.sd.grid <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order$utc.pred.response.sd)), nrow=utc.x.res*2)

utc_pred_lat_mean <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean)
utc_pred_lat_mean_df <- as.data.frame(mapply(c,utc_pred_lat_mean, utc_pred_na_15as_sub))
utc_pred_lat_mean_order <- orderBy(~y + x, data=utc_pred_lat_mean_df)
utc.pred.latent.mean.grid <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order$utc.pred.latent.mean)), nrow=utc.x.res*2)

utc_pred_lat_sd <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd)
utc_pred_lat_sd_df <- as.data.frame(mapply(c,utc_pred_lat_sd, utc_pred_na_15as_sub))
utc_pred_lat_sd_order <- orderBy(~y + x, data=utc_pred_lat_sd_df)
utc.pred.latent.sd.grid <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order$utc.pred.latent.sd)), nrow=utc.x.res*2)


# compare predictions with real data with a 3-arcsec resolution
stack_utc_pred_15as_real <- stacked_utc_pred_15as_fit[,c(1,2,8)] # 7 for fake 8 for real
stack_utc_pred_15as_real_df <- as.data.frame(mapply(c,stack_utc_pred_15as_real, utc_pred_na_15as_sub))
stack_utc_pred_15as_real_order <- orderBy(~y + x, data=stack_utc_pred_15as_real_df)
pop_utc_2015_15arcsec_grid.mat <- matrix(as.matrix(as.numeric(stack_utc_pred_15as_real_order$pop_utc_15_15as)), nrow = utc.x.res*2)

par(mfrow=c(1,1), mar = c(5, 5, 2, 2))
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,pop_utc_2015_15arcsec_grid.mat,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,2000)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,pop_utc_2015_15arcsec_grid.mat, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,utc.pred.response.mean.grid/4,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,2000)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,utc.pred.response.mean.grid/4, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,utc.pred.latent.mean.grid/4,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-1.5,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,utc.pred.latent.mean.grid/4, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid-utc.pred.latent.mean.grid)/4,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0.2,2.5)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid-utc.pred.latent.mean.grid)/4, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

res.15as<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean/4
RMSE.15as <- sqrt(mean(res.15as^2))
RMSE.15as 
correl.15as <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean/4) 
correl.15as 

####################################
# worst case scenario
# full model w/ sqrt transformation + vague prior
# est at 30 arcsec full data + pred at 3 arcsec
####################################
utc.stack.est.30as.sqrt <- inla.stack(data = list(sqrtpop=sqrt(as.numeric(stacked_utc_30as_fit$pop_utc_15_30as))),
                                      A = list(utc.A.est.30as, 1),
                                      effects = list(c(utc.s.index, list(Intercept=1)), 
                                                     list(ntl=as.numeric(stacked_utc_30as_fit$ntl_utc_15_30as), 
                                                          slope=as.numeric(stacked_utc_30as_fit$slope_utc_30as)
                                                     )),
                                      tag="est")

utc.stack.pred.15as.sqrt <- inla.stack(data = list(sqrtpop=NA),
                                       A = list(utc.A.pred.15as, 1),
                                       effects = list(c(utc.s.index, list(Intercept=1)),
                                                      list(ntl=as.numeric(stacked_utc_pred_15as_fit$ntl_utc_15_15as), 
                                                           slope=as.numeric(stacked_utc_pred_15as_fit$slope_utc_15as)                                                     
                                                      )),
                                       tag="pred.response") 

utc.join.stack.pred.15as.sqrt <- inla.stack(utc.stack.est.30as.sqrt, utc.stack.latent, utc.stack.pred.15as.sqrt)

utc.formula.sqrt <- sqrtpop ~ -1 + Intercept + ntl + slope + f(utc.spatial.field, model=utc.spde)

utc.output.pred.15as.sqrt <- inla(utc.formula.sqrt, 
                                  data=inla.stack.data(utc.join.stack.pred.15as.sqrt, spde=utc.spde), family="gaussian",
                                  control.fixed=utc.fixed.vague, 
                                  control.family=list(hyper=list(prec=utc.prec.vague)), 
                                  control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.sqrt), compute=TRUE),
                                  control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.sqrt)

utc.output.pred.15as.sqrt$dic$dic
utc.output.pred.15as.sqrt$waic$waic
slcpo(utc.output.pred.15as.sqrt)

### extract parameter estimates
utc.fixed.out.sqrt <- round(utc.output.pred.15as.sqrt$summary.fixed[,1:5],3) 
utc.fixed1_marg.sqrt  <- utc.output.pred.15as.sqrt$marginals.fixed[[1]]
utc.fixed2_marg.sqrt  <- utc.output.pred.15as.sqrt$marginals.fixed[[2]]
utc.fixed3_marg.sqrt  <- utc.output.pred.15as.sqrt$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.sqrt <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.sqrt$marginals.hyperpar[[1]])
utc.sigma2e_m1.sqrt <- inla.emarginal(function(x) x, utc.sigma2e_marg.sqrt)
utc.sigma2e_m2.sqrt <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.sqrt)
utc.sigma2e_stdev.sqrt <- sqrt(utc.sigma2e_m2.sqrt - utc.sigma2e_m1.sqrt^2)
utc.sigma2e_quantiles.sqrt <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.sqrt)

utc.spde.result.sqrt <- inla.spde2.result(utc.output.pred.15as.sqrt, name="utc.spatial.field", utc.spde)

# variance of spatial effects
utc.var.nom.marg.sqrt <- utc.spde.result.sqrt$marginals.variance.nominal[[1]]
utc.var.nom.m1.sqrt <- inla.emarginal(function(x) x, utc.var.nom.marg.sqrt)
utc.var.nom.m2.sqrt <- inla.emarginal(function(x) x^2, utc.var.nom.marg.sqrt)
utc.var.nom.stdev.sqrt <- sqrt(utc.var.nom.m2.sqrt - utc.var.nom.m1.sqrt^2)
utc.var.nom.quantiles.sqrt <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.sqrt)

utc.range.nom.marg.sqrt <- utc.spde.result.sqrt$marginals.range.nominal[[1]]
utc.range.nom.m1.sqrt <- inla.emarginal(function(x) x, utc.range.nom.marg.sqrt)
utc.range.nom.m2.sqrt <- inla.emarginal(function(x) x^2, utc.range.nom.marg.sqrt)
utc.range.nom.stdev.sqrt <- sqrt(utc.range.nom.m2.sqrt - utc.range.nom.m1.sqrt^2)
utc.range.nom.quantiles.sqrt <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.sqrt)

### extract predicted values
utc.index.pred.response.sqrt <- inla.stack.index(utc.join.stack.pred.15as.sqrt, "pred.response")$data 
utc.index.pred.latent.sqrt <- inla.stack.index(utc.join.stack.pred.15as.sqrt, "pred.latent")$data 

utc.pred.response.mean.sqrt <- utc.output.pred.15as.sqrt$summary.linear.predictor[utc.index.pred.response.sqrt, "mean"]
utc.pred.response.sd.sqrt <- utc.output.pred.15as.sqrt$summary.linear.predictor[utc.index.pred.response.sqrt, "sd"]

utc.pred.latent.mean.sqrt <- utc.output.pred.15as.sqrt$summary.linear.predictor[utc.index.pred.latent.sqrt, "mean"]
utc.pred.latent.sd.sqrt <- utc.output.pred.15as.sqrt$summary.linear.predictor[utc.index.pred.latent.sqrt, "sd"]

utc_pred_resp_mean.sqrt <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.mean.sqrt)
utc_pred_resp_mean_df.sqrt <- as.data.frame(mapply(c,utc_pred_resp_mean.sqrt, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.sqrt <- orderBy(~y + x, data=utc_pred_resp_mean_df.sqrt)
utc.pred.response.mean.grid.sqrt <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.sqrt$utc.pred.response.mean.sqrt)), nrow=utc.x.res*2)

utc_pred_resp_sd.sqrt <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.sqrt)
utc_pred_resp_sd_df.sqrt <- as.data.frame(mapply(c,utc_pred_resp_sd.sqrt, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.sqrt <- orderBy(~y + x, data=utc_pred_resp_sd_df.sqrt)
utc.pred.response.sd.grid.sqrt <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.sqrt$utc.pred.response.sd.sqrt)), nrow=utc.x.res*2)

utc_pred_lat_mean.sqrt <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.sqrt)
utc_pred_lat_mean_df.sqrt <- as.data.frame(mapply(c,utc_pred_lat_mean.sqrt, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.sqrt <- orderBy(~y + x, data=utc_pred_lat_mean_df.sqrt)
utc.pred.latent.mean.grid.sqrt <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.sqrt$utc.pred.latent.mean.sqrt)), nrow=utc.x.res*2)

utc_pred_lat_sd.sqrt <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.sqrt)
utc_pred_lat_sd_df.sqrt <- as.data.frame(mapply(c,utc_pred_lat_sd.sqrt, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.sqrt <- orderBy(~y + x, data=utc_pred_lat_sd_df.sqrt)
utc.pred.latent.sd.grid.sqrt <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.sqrt$utc.pred.latent.sd.sqrt)), nrow=utc.x.res*2)


# compare predictions with real data with a 3-arcsec resolution
pop_utc_2015_15arcsec_grid.mat.sqrt <- sqrt(matrix(as.matrix(as.numeric(stack_utc_pred_15as_real_order$pop_utc_15_15as)), nrow = utc.x.res*2))

par(mfrow=c(1,1), mar = c(5, 5, 2, 2))
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,pop_utc_2015_15arcsec_grid.mat.sqrt,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,80)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,pop_utc_2015_15arcsec_grid.mat.sqrt, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,utc.pred.response.mean.grid.sqrt/sqrt(4),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,80)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,utc.pred.response.mean.grid.sqrt/sqrt(4), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,utc.pred.latent.mean.grid.sqrt/sqrt(4),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-1.5,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,utc.pred.latent.mean.grid.sqrt/sqrt(4), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.sqrt-utc.pred.latent.mean.grid.sqrt)/sqrt(4),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0.2,2.5)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.sqrt-utc.pred.latent.mean.grid.sqrt)/sqrt(4), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.sqrt,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.sqrt, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.sqrt,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.sqrt, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

res.15as.sqrt<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - (utc.pred.response.mean.sqrt)^2/4
RMSE.15as.sqrt <- sqrt(mean(res.15as.sqrt^2))
RMSE.15as.sqrt 
correl.15as.sqrt <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), (utc.pred.response.mean.sqrt)^2/4) 
correl.15as.sqrt 

####################################
# worst case scenario
# full model w/ log transformation + vague prior
# est at 30 arcsec full data + pred at 3 arcsec
####################################
utc.stack.est.30as.log <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_30as_fit$pop_utc_15_30as)+1)),
                                     A = list(utc.A.est.30as, 1),
                                     effects = list(c(utc.s.index, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_30as_fit$ntl_utc_15_30as), 
                                                         slope=as.numeric(stacked_utc_30as_fit$slope_utc_30as)
                                                    )),
                                     tag="est")

utc.stack.pred.15as.log <- inla.stack(data = list(logpop=NA),
                                      A = list(utc.A.pred.15as, 1),
                                      effects = list(c(utc.s.index, list(Intercept=1)),
                                                     list(ntl=as.numeric(stacked_utc_pred_15as_fit$ntl_utc_15_15as), 
                                                          slope=as.numeric(stacked_utc_pred_15as_fit$slope_utc_15as)                                                     
                                                     )),
                                      tag="pred.response") 

utc.join.stack.pred.15as.log <- inla.stack(utc.stack.est.30as.log, utc.stack.latent, utc.stack.pred.15as.log)

utc.formula.log <- logpop ~ -1 + Intercept + ntl + slope + f(utc.spatial.field, model=utc.spde)

utc.output.pred.15as.log <- inla(utc.formula.log, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.log, spde=utc.spde), family="gaussian",
                                 control.fixed=utc.fixed.vague, 
                                 control.family=list(hyper=list(prec=utc.prec.vague)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.log), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.log)

utc.output.pred.15as.log$dic$dic
utc.output.pred.15as.log$waic$waic
slcpo(utc.output.pred.15as.log)

### extract parameter estimates
utc.fixed.out.log <- round(utc.output.pred.15as.log$summary.fixed[,1:5],3) 
utc.fixed1_marg.log  <- utc.output.pred.15as.log$marginals.fixed[[1]]
utc.fixed2_marg.log  <- utc.output.pred.15as.log$marginals.fixed[[2]]
utc.fixed3_marg.log  <- utc.output.pred.15as.log$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.log <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.log$marginals.hyperpar[[1]])
utc.sigma2e_m1.log <- inla.emarginal(function(x) x, utc.sigma2e_marg.log)
utc.sigma2e_m2.log <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.log)
utc.sigma2e_stdev.log <- sqrt(utc.sigma2e_m2.log - utc.sigma2e_m1.log^2)
utc.sigma2e_quantiles.log <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.log)

utc.spde.result.log <- inla.spde2.result(utc.output.pred.15as.log, name="utc.spatial.field", utc.spde)

# variance of spatial effects
utc.var.nom.marg.log <- utc.spde.result.log$marginals.variance.nominal[[1]]
utc.var.nom.m1.log <- inla.emarginal(function(x) x, utc.var.nom.marg.log)
utc.var.nom.m2.log <- inla.emarginal(function(x) x^2, utc.var.nom.marg.log)
utc.var.nom.stdev.log <- sqrt(utc.var.nom.m2.log - utc.var.nom.m1.log^2)
utc.var.nom.quantiles.log <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.log)

utc.range.nom.marg.log <- utc.spde.result.log$marginals.range.nominal[[1]]
utc.range.nom.m1.log <- inla.emarginal(function(x) x, utc.range.nom.marg.log)
utc.range.nom.m2.log <- inla.emarginal(function(x) x^2, utc.range.nom.marg.log)
utc.range.nom.stdev.log <- sqrt(utc.range.nom.m2.log - utc.range.nom.m1.log^2)
utc.range.nom.quantiles.log <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.log)

### extract predicted values
utc.index.pred.response.log <- inla.stack.index(utc.join.stack.pred.15as.log, "pred.response")$data 
utc.index.pred.latent.log <- inla.stack.index(utc.join.stack.pred.15as.log, "pred.latent")$data 

utc.pred.response.mean.log <- utc.output.pred.15as.log$summary.linear.predictor[utc.index.pred.response.log, "mean"]
utc.pred.response.sd.log <- utc.output.pred.15as.log$summary.linear.predictor[utc.index.pred.response.log, "sd"]

utc.pred.latent.mean.log <- utc.output.pred.15as.log$summary.linear.predictor[utc.index.pred.latent.log, "mean"]
utc.pred.latent.sd.log <- utc.output.pred.15as.log$summary.linear.predictor[utc.index.pred.latent.log, "sd"]

utc_pred_resp_mean.log <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.log.mod/4+1))
utc_pred_resp_mean_df.log <- as.data.frame(mapply(c,utc_pred_resp_mean.log, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.log <- orderBy(~y + x, data=utc_pred_resp_mean_df.log)
utc.pred.response.mean.grid.log <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.log$`log(utc.pred.response.mean.log.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.log <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.log)
utc_pred_resp_sd_df.log <- as.data.frame(mapply(c,utc_pred_resp_sd.log, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.log <- orderBy(~y + x, data=utc_pred_resp_sd_df.log)
utc.pred.response.sd.grid.log <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.log$utc.pred.response.sd.log)), nrow=utc.x.res*2)

utc_pred_lat_mean.log <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.log)
utc_pred_lat_mean_df.log <- as.data.frame(mapply(c,utc_pred_lat_mean.log, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.log <- orderBy(~y + x, data=utc_pred_lat_mean_df.log)
utc.pred.latent.mean.grid.log <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.log$utc.pred.latent.mean.log)), nrow=utc.x.res*2)

utc_pred_lat_sd.log <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.log)
utc_pred_lat_sd_df.log <- as.data.frame(mapply(c,utc_pred_lat_sd.log, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.log <- orderBy(~y + x, data=utc_pred_lat_sd_df.log)
utc.pred.latent.sd.grid.log <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.log$utc.pred.latent.sd.log)), nrow=utc.x.res*2)

# compare predictions with real data with a 3-arcsec resolution
pop_utc_2015_15arcsec_grid.mat.log <- log(matrix(as.matrix(as.numeric(stack_utc_pred_15as_real_order$pop_utc_15_15as)+1), nrow = utc.x.res*2))

par(mfrow=c(1,1), mar = c(5, 5, 2, 2))
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,pop_utc_2015_15arcsec_grid.mat.log,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,pop_utc_2015_15arcsec_grid.mat.log, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.log),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.log), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.log-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-1.5,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.log-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.log-utc.pred.latent.mean.grid.log)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0.2,2.5)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.log-utc.pred.latent.mean.grid.log)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.log,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.log, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.log,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.log, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.log.mod <- exp(utc.pred.response.mean.log)-1
utc.pred.response.mean.log.mod[utc.pred.response.mean.log.mod<0]<-0
res.15as.log<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.log.mod/4
RMSE.15as.log <- sqrt(mean(res.15as.log^2))
RMSE.15as.log 
correl.15as.log <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.log.mod/4) 
correl.15as.log 

####################################
# worst case scenario
# model w/ log transformation w/o covariates + vague prior
# est at 30 arcsec full data + pred at 3 arcsec
####################################
utc.stack.est.30as.noboth <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_30as_fit$pop_utc_15_30as)+1)),
                                        A = list(utc.A.est.30as),
                                        effects = list(c(utc.s.index, list(Intercept=1))),
                                        tag="est")

utc.stack.pred.15as.noboth <- inla.stack(data = list(logpop=NA),
                                         A = list(utc.A.pred.15as),
                                         effects = list(c(utc.s.index, list(Intercept=1))),
                                         tag="pred.response") 

utc.join.stack.pred.15as.noboth <- inla.stack(utc.stack.est.30as.noboth, utc.stack.latent, utc.stack.pred.15as.noboth)

utc.formula.noboth <- logpop ~ -1 + Intercept + f(utc.spatial.field, model=utc.spde)

utc.output.pred.15as.noboth <- inla(utc.formula.noboth, 
                                    data=inla.stack.data(utc.join.stack.pred.15as.noboth, spde=utc.spde), family="gaussian",
                                    control.fixed=utc.fixed.vague, 
                                    control.family=list(hyper=list(prec=utc.prec.vague)), 
                                    control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.noboth), compute=TRUE),
                                    control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.noboth)

utc.output.pred.15as.noboth$dic$dic
utc.output.pred.15as.noboth$waic$waic
slcpo(utc.output.pred.15as.noboth)

### extract parameter estimates
utc.fixed.out.noboth <- round(utc.output.pred.15as.noboth$summary.fixed[,1:5],3) 
utc.fixed1_marg.noboth <- utc.output.pred.15as.noboth$marginals.fixed[[1]]

# variance of unstructured residual
utc.sigma2e_marg.noboth <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.noboth$marginals.hyperpar[[1]])
utc.sigma2e_m1.noboth <- inla.emarginal(function(x) x, utc.sigma2e_marg.noboth)
utc.sigma2e_m2.noboth <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.noboth)
utc.sigma2e_stdev.noboth <- sqrt(utc.sigma2e_m2.noboth - utc.sigma2e_m1.noboth^2)
utc.sigma2e_quantiles.noboth <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.noboth)

utc.spde.result.noboth <- inla.spde2.result(utc.output.pred.15as.noboth, name="utc.spatial.field", utc.spde)

# variance of spatial effects
utc.var.nom.marg.noboth <- utc.spde.result.noboth$marginals.variance.nominal[[1]]
utc.var.nom.m1.noboth <- inla.emarginal(function(x) x, utc.var.nom.marg.noboth)
utc.var.nom.m2.noboth <- inla.emarginal(function(x) x^2, utc.var.nom.marg.noboth)
utc.var.nom.stdev.noboth <- sqrt(utc.var.nom.m2.noboth - utc.var.nom.m1.noboth^2)
utc.var.nom.quantiles.noboth <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.noboth)

utc.range.nom.marg.noboth <- utc.spde.result.noboth$marginals.range.nominal[[1]]
utc.range.nom.m1.noboth <- inla.emarginal(function(x) x, utc.range.nom.marg.noboth)
utc.range.nom.m2.noboth <- inla.emarginal(function(x) x^2, utc.range.nom.marg.noboth)
utc.range.nom.stdev.noboth <- sqrt(utc.range.nom.m2.noboth - utc.range.nom.m1.noboth^2)
utc.range.nom.quantiles.noboth <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.noboth)

### extract predicted values
utc.index.pred.response.noboth <- inla.stack.index(utc.join.stack.pred.15as.noboth, "pred.response")$data 
utc.index.pred.latent.noboth <- inla.stack.index(utc.join.stack.pred.15as.noboth, "pred.latent")$data 

utc.pred.response.mean.noboth <- utc.output.pred.15as.noboth$summary.linear.predictor[utc.index.pred.response.noboth, "mean"]
utc.pred.response.sd.noboth <- utc.output.pred.15as.noboth$summary.linear.predictor[utc.index.pred.response.noboth, "sd"]

utc.pred.latent.mean.noboth <- utc.output.pred.15as.noboth$summary.linear.predictor[utc.index.pred.latent.noboth, "mean"]
utc.pred.latent.sd.noboth <- utc.output.pred.15as.noboth$summary.linear.predictor[utc.index.pred.latent.noboth, "sd"]

utc_pred_resp_mean.noboth <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.mean.noboth)
utc_pred_resp_mean_df.noboth <- as.data.frame(mapply(c,utc_pred_resp_mean.noboth, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.noboth <- orderBy(~y + x, data=utc_pred_resp_mean_df.noboth)
utc.pred.response.mean.grid.noboth <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.noboth$utc.pred.response.mean.noboth)), nrow=utc.x.res*2)

utc_pred_resp_sd.noboth <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.noboth)
utc_pred_resp_sd_df.noboth <- as.data.frame(mapply(c,utc_pred_resp_sd.noboth, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.noboth <- orderBy(~y + x, data=utc_pred_resp_sd_df.noboth)
utc.pred.response.sd.grid.noboth <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.noboth$utc.pred.response.sd.noboth)), nrow=utc.x.res*2)

utc_pred_lat_mean.noboth <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.noboth)
utc_pred_lat_mean_df.noboth <- as.data.frame(mapply(c,utc_pred_lat_mean.noboth, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.noboth <- orderBy(~y + x, data=utc_pred_lat_mean_df.noboth)
utc.pred.latent.mean.grid.noboth <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.noboth$utc.pred.latent.mean.noboth)), nrow=utc.x.res*2)

utc_pred_lat_sd.noboth <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.noboth)
utc_pred_lat_sd_df.noboth <- as.data.frame(mapply(c,utc_pred_lat_sd.noboth, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.noboth <- orderBy(~y + x, data=utc_pred_lat_sd_df.noboth)
utc.pred.latent.sd.grid.noboth <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.noboth$utc.pred.latent.sd.noboth)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.noboth-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-1.1,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.noboth-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.noboth-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-1.5,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.noboth-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.noboth-utc.pred.latent.mean.grid.noboth)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0.2,2.5)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.noboth-utc.pred.latent.mean.grid.noboth)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.noboth,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.noboth, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.noboth,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.noboth, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.noboth.mod <- exp(utc.pred.response.mean.noboth)-1
utc.pred.response.mean.noboth.mod[utc.pred.response.mean.noboth.mod<0]<-0
res.15as.noboth<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.noboth.mod/4
RMSE.15as.noboth <- sqrt(mean(res.15as.noboth^2))
RMSE.15as.noboth 
correl.15as.noboth <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.noboth.mod/4) 
correl.15as.noboth 

####################################
# worst case scenario + best prior
# full model w/ log transformation + normal prior est at 3 arcsec full data
# est at 30 arcsec full data + pred at 3 arcsec
####################################
# all data from utc at 15as as informed prior
# calculate prior distribution of parameters from supportive pop 
utc.stack.est.15as.log <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.pred.15as, 1),
                                     effects = list(c(utc.s.index, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_pred_15as_fit$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_pred_15as_fit$slope_utc_15as)
                                                    )),
                                     tag="est")


utc.output.est.15as.log <- inla(utc.formula.log, 
                                data=inla.stack.data(utc.stack.est.15as.log, spde=utc.spde), family="gaussian", 
                                control.fixed=utc.fixed.vague, 
                                control.family=list(hyper=list(prec=utc.prec.vague)), 
                                control.predictor=list(A=inla.stack.A(utc.stack.est.15as.log), compute=TRUE),
                                control.compute=list(cpo=TRUE, dic=TRUE), verbose=TRUE)

summary(utc.output.est.15as.log)

### extract parameter estimates for model w/ ntl and slope and w/ log transformation
utc.15as.fixed.out <- round(utc.output.est.15as.log$summary.fixed[,1:5],3) 

# variance of unstructured residual
utc.15as.sigma2e_marg <- inla.tmarginal(function(x) 1/x, utc.output.est.15as.log$marginals.hyperpar[[1]])
utc.15as.sigma2e_m1 <- inla.emarginal(function(x) x, utc.15as.sigma2e_marg)
utc.15as.sigma2e_m2 <- inla.emarginal(function(x) x^2, utc.15as.sigma2e_marg)
utc.15as.sigma2e_stdev <- sqrt(utc.15as.sigma2e_m2 - utc.15as.sigma2e_m1^2)
utc.15as.sigma2e_quantiles <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.15as.sigma2e_marg)

utc.15as.spde.result <- inla.spde2.result(utc.output.est.15as.log, name="utc.spatial.field", spde=utc.spde)

# variance of spatial effects
utc.15as.var.nom.marg <- utc.15as.spde.result$marginals.variance.nominal[[1]]
utc.15as.var.nom.m1 <- inla.emarginal(function(x) x, utc.15as.var.nom.marg)
utc.15as.var.nom.m2 <- inla.emarginal(function(x) x^2, utc.15as.var.nom.marg)
utc.15as.var.nom.stdev <- sqrt(utc.15as.var.nom.m2 - utc.15as.var.nom.m1^2)
utc.15as.var.nom.quantiles <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.15as.var.nom.marg)

utc.15as.range.nom.marg <- utc.15as.spde.result$marginals.range.nominal[[1]]
utc.15as.range.nom.m1 <- inla.emarginal(function(x) x, utc.15as.range.nom.marg)
utc.15as.range.nom.m2 <- inla.emarginal(function(x) x^2, utc.15as.range.nom.marg)
utc.15as.range.nom.stdev <- sqrt(utc.15as.range.nom.m2 - utc.15as.range.nom.m1^2)
utc.15as.range.nom.quantiles <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.15as.range.nom.marg)

# set up priors
utc.15as.spde.normal <- inla.spde2.matern(mesh=utc.mesh, alpha=2, 
                                          B.tau=matrix(c(0, 1, 0),nrow=1,ncol=3),
                                          B.kappa=matrix(c(0, 0, 1),nrow=1,ncol=3), 
                                          theta.prior.mean=c(utc.15as.spde.result$summary.theta[[2]][1], 
                                                             utc.15as.spde.result$summary.theta[[2]][2]),
                                          theta.prior.prec=c(1/(utc.15as.spde.result$summary.theta[[3]][1]^2), 
                                                             1/(utc.15as.spde.result$summary.theta[[3]][2]^2)))

utc.15as.s.index.normal <- inla.spde.make.index(name="utc.15as.spatial.field.normal", n.spde=utc.15as.spde.normal$n.spde) 

utc.15as.spde.pc <- inla.spde2.pcmatern(utc.mesh, alpha = 2, 
                                        prior.range = c(utc.15as.range.nom.quantiles[1], 0.05),
                                        prior.sigma = c(sqrt(utc.15as.var.nom.quantiles[3]), 0.05))

utc.15as.s.index.pc <- inla.spde.make.index(name = "utc.15as.spatial.field.pc", n.spde = utc.15as.spde.pc$n.spde)

utc.15as.prec.normal <- list(initial=utc.output.est.15as.log$internal.summary.hyperpar[[1]][1], prior="normal", 
                             param=c(utc.output.est.15as.log$internal.summary.hyperpar[[1]][1], 
                                     1/(utc.output.est.15as.log$internal.summary.hyperpar[[2]][1]^2)), fixed=FALSE)

utc.15as.prec.pc <- list(prec=list(prior="pc.prec", param=c(sqrt(utc.15as.sigma2e_quantiles[3]), 0.05)))

utc.15as.fixed.normal <- list(mean.intercept = 0, 
                              prec.intercept = 0, 
                              mean = list(ntl=utc.output.est.15as.log$summary.fixed[[1]][2], 
                                          slope=utc.output.est.15as.log$summary.fixed[[1]][3]),
                              prec = list(ntl=1/(utc.output.est.15as.log$summary.fixed[[2]][2]^2),
                                          slope=1/(utc.output.est.15as.log$summary.fixed[[2]][3]^2)))

utc.stack.est.30as.wsbp.n <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_30as_fit$pop_utc_15_30as)+1)),
                                        A = list(utc.A.est.30as, 1),
                                        effects = list(c(utc.15as.s.index.normal, list(Intercept=1)), 
                                                       list(ntl=as.numeric(stacked_utc_30as_fit$ntl_utc_15_30as), 
                                                            slope=as.numeric(stacked_utc_30as_fit$slope_utc_30as)
                                                       )),
                                        tag="est")

utc.stack.pred.15as.bp.n <- inla.stack(data = list(logpop=NA),
                                       A = list(utc.A.pred.15as, 1),
                                       effects = list(c(utc.15as.s.index.normal, list(Intercept=1)),
                                                      list(ntl=as.numeric(stacked_utc_pred_15as_fit$ntl_utc_15_15as), 
                                                           slope=as.numeric(stacked_utc_pred_15as_fit$slope_utc_15as)                                                     
                                                      )),
                                       tag="pred.response") 

utc.stack.latent.bp.n <- inla.stack(data = list(xi=NA),
                                    A = list(utc.A.pred.15as),
                                    effects = list(utc.15as.s.index.normal),
                                    tag="pred.latent") 

utc.join.stack.pred.15as.wsbp.n <- inla.stack(utc.stack.est.30as.wsbp.n, utc.stack.latent.bp.n, utc.stack.pred.15as.bp.n)

utc.formula.log.15as.normal <- logpop ~ -1 + Intercept + ntl + slope + f(utc.15as.spatial.field.normal, model=utc.15as.spde.normal)

utc.output.pred.15as.wsbp.n <- inla(utc.formula.log.15as.normal, 
                                    data=inla.stack.data(utc.join.stack.pred.15as.wsbp.n, spde=utc.15as.spde.normal), family="gaussian",
                                    control.fixed=utc.15as.fixed.normal, 
                                    control.family=list(hyper=list(prec=utc.15as.prec.normal)), 
                                    control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.wsbp.n), compute=TRUE),
                                    control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.wsbp.n)

utc.output.pred.15as.wsbp.n$dic$dic
utc.output.pred.15as.wsbp.n$waic$waic
slcpo(utc.output.pred.15as.wsbp.n)

### extract parameter estimates
utc.fixed.out.wsbp.n <- round(utc.output.pred.15as.wsbp.n$summary.fixed[,1:5],3) 
utc.fixed1_marg.wsbp.n  <- utc.output.pred.15as.wsbp.n$marginals.fixed[[1]]
utc.fixed2_marg.wsbp.n  <- utc.output.pred.15as.wsbp.n$marginals.fixed[[2]]
utc.fixed3_marg.wsbp.n  <- utc.output.pred.15as.wsbp.n$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.wsbp.n <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.wsbp.n$marginals.hyperpar[[1]])
utc.sigma2e_m1.wsbp.n <- inla.emarginal(function(x) x, utc.sigma2e_marg.wsbp.n)
utc.sigma2e_m2.wsbp.n <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.wsbp.n)
utc.sigma2e_stdev.wsbp.n <- sqrt(utc.sigma2e_m2.wsbp.n - utc.sigma2e_m1.wsbp.n^2)
utc.sigma2e_quantiles.wsbp.n <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.wsbp.n)

utc.spde.result.wsbp.n <- inla.spde2.result(utc.output.pred.15as.wsbp.n, name="utc.15as.spatial.field.normal", utc.15as.spde.normal)

# variance of spatial effects
utc.var.nom.marg.wsbp.n <- utc.spde.result.wsbp.n$marginals.variance.nominal[[1]]
utc.var.nom.m1.wsbp.n <- inla.emarginal(function(x) x, utc.var.nom.marg.wsbp.n)
utc.var.nom.m2.wsbp.n <- inla.emarginal(function(x) x^2, utc.var.nom.marg.wsbp.n)
utc.var.nom.stdev.wsbp.n <- sqrt(utc.var.nom.m2.wsbp.n - utc.var.nom.m1.wsbp.n^2)
utc.var.nom.quantiles.wsbp.n <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.wsbp.n)

utc.range.nom.marg.wsbp.n <- utc.spde.result.wsbp.n$marginals.range.nominal[[1]]
utc.range.nom.m1.wsbp.n <- inla.emarginal(function(x) x, utc.range.nom.marg.wsbp.n)
utc.range.nom.m2.wsbp.n <- inla.emarginal(function(x) x^2, utc.range.nom.marg.wsbp.n)
utc.range.nom.stdev.wsbp.n <- sqrt(utc.range.nom.m2.wsbp.n - utc.range.nom.m1.wsbp.n^2)
utc.range.nom.quantiles.wsbp.n <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.wsbp.n)

### extract predicted values
utc.index.pred.response.wsbp.n <- inla.stack.index(utc.join.stack.pred.15as.wsbp.n, "pred.response")$data 
utc.index.pred.latent.wsbp.n <- inla.stack.index(utc.join.stack.pred.15as.wsbp.n, "pred.latent")$data 

utc.pred.response.mean.wsbp.n <- utc.output.pred.15as.wsbp.n$summary.linear.predictor[utc.index.pred.response.wsbp.n, "mean"]
utc.pred.response.sd.wsbp.n <- utc.output.pred.15as.wsbp.n$summary.linear.predictor[utc.index.pred.response.wsbp.n, "sd"]

utc.pred.latent.mean.wsbp.n <- utc.output.pred.15as.wsbp.n$summary.linear.predictor[utc.index.pred.latent.wsbp.n, "mean"]
utc.pred.latent.sd.wsbp.n <- utc.output.pred.15as.wsbp.n$summary.linear.predictor[utc.index.pred.latent.wsbp.n, "sd"]

utc_pred_resp_mean.wsbp.n <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.mean.wsbp.n)
utc_pred_resp_mean_df.wsbp.n <- as.data.frame(mapply(c,utc_pred_resp_mean.wsbp.n, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.wsbp.n <- orderBy(~y + x, data=utc_pred_resp_mean_df.wsbp.n)
utc.pred.response.mean.grid.wsbp.n <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.wsbp.n$utc.pred.response.mean.wsbp.n)), nrow=utc.x.res*2)

utc_pred_resp_sd.wsbp.n <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.wsbp.n)
utc_pred_resp_sd_df.wsbp.n <- as.data.frame(mapply(c,utc_pred_resp_sd.wsbp.n, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.wsbp.n <- orderBy(~y + x, data=utc_pred_resp_sd_df.wsbp.n)
utc.pred.response.sd.grid.wsbp.n <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.wsbp.n$utc.pred.response.sd.wsbp.n)), nrow=utc.x.res*2)

utc_pred_lat_mean.wsbp.n <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.wsbp.n)
utc_pred_lat_mean_df.wsbp.n <- as.data.frame(mapply(c,utc_pred_lat_mean.wsbp.n, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.wsbp.n <- orderBy(~y + x, data=utc_pred_lat_mean_df.wsbp.n)
utc.pred.latent.mean.grid.wsbp.n <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.wsbp.n$utc.pred.latent.mean.wsbp.n)), nrow=utc.x.res*2)

utc_pred_lat_sd.wsbp.n <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.wsbp.n)
utc_pred_lat_sd_df.wsbp.n <- as.data.frame(mapply(c,utc_pred_lat_sd.wsbp.n, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.wsbp.n <- orderBy(~y + x, data=utc_pred_lat_sd_df.wsbp.n)
utc.pred.latent.sd.grid.wsbp.n <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.wsbp.n$utc.pred.latent.sd.wsbp.n)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wsbp.n-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wsbp.n-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wsbp.n-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wsbp.n-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wsbp.n-utc.pred.latent.mean.grid.wsbp.n)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wsbp.n-utc.pred.latent.mean.grid.wsbp.n)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wsbp.n,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wsbp.n, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wsbp.n,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wsbp.n, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.wsbp.n.mod <- exp(utc.pred.response.mean.wsbp.n)-1
utc.pred.response.mean.wsbp.n.mod[utc.pred.response.mean.wsbp.n.mod<0]<-0
res.15as.wsbp.n<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.wsbp.n.mod/4
RMSE.15as.wsbp.n <- sqrt(mean(res.15as.wsbp.n^2))
RMSE.15as.wsbp.n 
correl.15as.wsbp.n <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.wsbp.n.mod/4) 
correl.15as.wsbp.n 

####################################
# worst case scenario + best prior
# full model w/ log transformation + pc prior est at 3 arcsec full data
# est at 30 arcsec full data + pred at 3 arcsec
####################################
utc.stack.est.30as.wsbp.p <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_30as_fit$pop_utc_15_30as)+1)),
                                        A = list(utc.A.est.30as, 1),
                                        effects = list(c(utc.15as.s.index.pc, list(Intercept=1)), 
                                                       list(ntl=as.numeric(stacked_utc_30as_fit$ntl_utc_15_30as), 
                                                            slope=as.numeric(stacked_utc_30as_fit$slope_utc_30as)
                                                       )),
                                        tag="est")

utc.stack.pred.15as.bp.p <- inla.stack(data = list(logpop=NA),
                                       A = list(utc.A.pred.15as, 1),
                                       effects = list(c(utc.15as.s.index.pc, list(Intercept=1)),
                                                      list(ntl=as.numeric(stacked_utc_pred_15as_fit$ntl_utc_15_15as), 
                                                           slope=as.numeric(stacked_utc_pred_15as_fit$slope_utc_15as)                                                     
                                                      )),
                                       tag="pred.response") 

utc.stack.latent.bp.p <- inla.stack(data = list(xi=NA),
                                    A = list(utc.A.pred.15as),
                                    effects = list(utc.15as.s.index.pc),
                                    tag="pred.latent") 

utc.join.stack.pred.15as.wsbp.p <- inla.stack(utc.stack.est.30as.wsbp.p, utc.stack.latent.bp.p, utc.stack.pred.15as.bp.p)

utc.formula.log.15as.pc <- logpop ~ -1 + Intercept + ntl + slope + f(utc.15as.spatial.field.pc, model=utc.15as.spde.pc)

utc.output.pred.15as.wsbp.p <- inla(utc.formula.log.15as.pc, 
                                    data=inla.stack.data(utc.join.stack.pred.15as.wsbp.p, spde=utc.15as.spde.pc), family="gaussian",
                                    control.fixed=utc.15as.fixed.normal, 
                                    control.family=list(hyper=list(prec=utc.15as.prec.pc)), 
                                    control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.wsbp.p), compute=TRUE),
                                    control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.wsbp.p)

utc.output.pred.15as.wsbp.p$dic$dic
utc.output.pred.15as.wsbp.p$waic$waic
slcpo(utc.output.pred.15as.wsbp.p)

### extract parameter estimates
utc.fixed.out.wsbp.p <- round(utc.output.pred.15as.wsbp.p$summary.fixed[,1:5],3) 
utc.fixed1_marg.wsbp.p  <- utc.output.pred.15as.wsbp.p$marginals.fixed[[1]]
utc.fixed2_marg.wsbp.p  <- utc.output.pred.15as.wsbp.p$marginals.fixed[[2]]
utc.fixed3_marg.wsbp.p  <- utc.output.pred.15as.wsbp.p$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.wsbp.p <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.wsbp.p$marginals.hyperpar[[1]])
utc.sigma2e_m1.wsbp.p <- inla.emarginal(function(x) x, utc.sigma2e_marg.wsbp.p)
utc.sigma2e_m2.wsbp.p <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.wsbp.p)
utc.sigma2e_stdev.wsbp.p <- sqrt(utc.sigma2e_m2.wsbp.p - utc.sigma2e_m1.wsbp.p^2)
utc.sigma2e_quantiles.wsbp.p <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.wsbp.p)

utc.spde.result.wsbp.p <- inla.spde2.result(utc.output.pred.15as.wsbp.p, name="utc.15as.spatial.field.pc", utc.15as.spde.pc)

# variance of spatial effects
utc.var.nom.marg.wsbp.p <- utc.spde.result.wsbp.p$marginals.variance.nominal[[1]]
utc.var.nom.m1.wsbp.p <- inla.emarginal(function(x) x, utc.var.nom.marg.wsbp.p)
utc.var.nom.m2.wsbp.p <- inla.emarginal(function(x) x^2, utc.var.nom.marg.wsbp.p)
utc.var.nom.stdev.wsbp.p <- sqrt(utc.var.nom.m2.wsbp.p - utc.var.nom.m1.wsbp.p^2)
utc.var.nom.quantiles.wsbp.p <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.wsbp.p)

utc.range.nom.marg.wsbp.p <- utc.spde.result.wsbp.p$marginals.range.nominal[[1]]
utc.range.nom.m1.wsbp.p <- inla.emarginal(function(x) x, utc.range.nom.marg.wsbp.p)
utc.range.nom.m2.wsbp.p <- inla.emarginal(function(x) x^2, utc.range.nom.marg.wsbp.p)
utc.range.nom.stdev.wsbp.p <- sqrt(utc.range.nom.m2.wsbp.p - utc.range.nom.m1.wsbp.p^2)
utc.range.nom.quantiles.wsbp.p <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.wsbp.p)

### extract predicted values
utc.index.pred.response.wsbp.p <- inla.stack.index(utc.join.stack.pred.15as.wsbp.p, "pred.response")$data 
utc.index.pred.latent.wsbp.p <- inla.stack.index(utc.join.stack.pred.15as.wsbp.p, "pred.latent")$data 

utc.pred.response.mean.wsbp.p <- utc.output.pred.15as.wsbp.p$summary.linear.predictor[utc.index.pred.response.wsbp.p, "mean"]
utc.pred.response.sd.wsbp.p <- utc.output.pred.15as.wsbp.p$summary.linear.predictor[utc.index.pred.response.wsbp.p, "sd"]

utc.pred.latent.mean.wsbp.p <- utc.output.pred.15as.wsbp.p$summary.linear.predictor[utc.index.pred.latent.wsbp.p, "mean"]
utc.pred.latent.sd.wsbp.p <- utc.output.pred.15as.wsbp.p$summary.linear.predictor[utc.index.pred.latent.wsbp.p, "sd"]

utc_pred_resp_mean.wsbp.p <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.mean.wsbp.p)
utc_pred_resp_mean_df.wsbp.p <- as.data.frame(mapply(c,utc_pred_resp_mean.wsbp.p, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.wsbp.p <- orderBy(~y + x, data=utc_pred_resp_mean_df.wsbp.p)
utc.pred.response.mean.grid.wsbp.p <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.wsbp.p$utc.pred.response.mean.wsbp.p)), nrow=utc.x.res*2)

utc_pred_resp_sd.wsbp.p <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.wsbp.p)
utc_pred_resp_sd_df.wsbp.p <- as.data.frame(mapply(c,utc_pred_resp_sd.wsbp.p, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.wsbp.p <- orderBy(~y + x, data=utc_pred_resp_sd_df.wsbp.p)
utc.pred.response.sd.grid.wsbp.p <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.wsbp.p$utc.pred.response.sd.wsbp.p)), nrow=utc.x.res*2)

utc_pred_lat_mean.wsbp.p <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.wsbp.p)
utc_pred_lat_mean_df.wsbp.p <- as.data.frame(mapply(c,utc_pred_lat_mean.wsbp.p, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.wsbp.p <- orderBy(~y + x, data=utc_pred_lat_mean_df.wsbp.p)
utc.pred.latent.mean.grid.wsbp.p <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.wsbp.p$utc.pred.latent.mean.wsbp.p)), nrow=utc.x.res*2)

utc_pred_lat_sd.wsbp.p <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.wsbp.p)
utc_pred_lat_sd_df.wsbp.p <- as.data.frame(mapply(c,utc_pred_lat_sd.wsbp.p, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.wsbp.p <- orderBy(~y + x, data=utc_pred_lat_sd_df.wsbp.p)
utc.pred.latent.sd.grid.wsbp.p <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.wsbp.p$utc.pred.latent.sd.wsbp.p)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wsbp.p-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wsbp.p-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wsbp.p-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wsbp.p-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wsbp.p-utc.pred.latent.mean.grid.wsbp.p)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wsbp.p-utc.pred.latent.mean.grid.wsbp.p)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wsbp.p,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wsbp.p, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wsbp.p,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wsbp.p, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.wsbp.p.mod <- exp(utc.pred.response.mean.wsbp.p)-1
utc.pred.response.mean.wsbp.p.mod[utc.pred.response.mean.wsbp.p.mod<0]<-0
res.15as.wsbp.p<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.wsbp.p.mod/4
RMSE.15as.wsbp.p <- sqrt(mean(res.15as.wsbp.p^2))
RMSE.15as.wsbp.p 
correl.15as.wsbp.p <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.wsbp.p.mod/4) 
correl.15as.wsbp.p 

####################################
# best case scenario
# full model w/ log transformation + vague prior
# est at 3 arcsec full data + pred at 3 arcsec
####################################
utc.join.stack.pred.15as.bsvp <- inla.stack(utc.stack.est.15as.log, utc.stack.latent, utc.stack.pred.15as.log)

utc.output.pred.15as.bsvp <- inla(utc.formula.log, 
                                  data=inla.stack.data(utc.join.stack.pred.15as.bsvp, spde=utc.spde), family="gaussian",
                                  control.fixed=utc.fixed.vague, 
                                  control.family=list(hyper=list(prec=utc.prec.vague)), 
                                  control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.bsvp), compute=TRUE),
                                  control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.bsvp)

utc.output.pred.15as.bsvp$dic$dic
utc.output.pred.15as.bsvp$waic$waic
slcpo(utc.output.pred.15as.bsvp)

### extract parameter estimates
utc.fixed.out.bsvp <- round(utc.output.pred.15as.bsvp$summary.fixed[,1:5],3) 
utc.fixed1_marg.bsvp  <- utc.output.pred.15as.bsvp$marginals.fixed[[1]]
utc.fixed2_marg.bsvp  <- utc.output.pred.15as.bsvp$marginals.fixed[[2]]
utc.fixed3_marg.bsvp  <- utc.output.pred.15as.bsvp$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.bsvp <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.bsvp$marginals.hyperpar[[1]])
utc.sigma2e_m1.bsvp <- inla.emarginal(function(x) x, utc.sigma2e_marg.bsvp)
utc.sigma2e_m2.bsvp <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.bsvp)
utc.sigma2e_stdev.bsvp <- sqrt(utc.sigma2e_m2.bsvp - utc.sigma2e_m1.bsvp^2)
utc.sigma2e_quantiles.bsvp <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.bsvp)

utc.spde.result.bsvp <- inla.spde2.result(utc.output.pred.15as.bsvp, name="utc.spatial.field", utc.spde)

# variance of spatial effects
utc.var.nom.marg.bsvp <- utc.spde.result.bsvp$marginals.variance.nominal[[1]]
utc.var.nom.m1.bsvp <- inla.emarginal(function(x) x, utc.var.nom.marg.bsvp)
utc.var.nom.m2.bsvp <- inla.emarginal(function(x) x^2, utc.var.nom.marg.bsvp)
utc.var.nom.stdev.bsvp <- sqrt(utc.var.nom.m2.bsvp - utc.var.nom.m1.bsvp^2)
utc.var.nom.quantiles.bsvp <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.bsvp)

utc.range.nom.marg.bsvp <- utc.spde.result.bsvp$marginals.range.nominal[[1]]
utc.range.nom.m1.bsvp <- inla.emarginal(function(x) x, utc.range.nom.marg.bsvp)
utc.range.nom.m2.bsvp <- inla.emarginal(function(x) x^2, utc.range.nom.marg.bsvp)
utc.range.nom.stdev.bsvp <- sqrt(utc.range.nom.m2.bsvp - utc.range.nom.m1.bsvp^2)
utc.range.nom.quantiles.bsvp <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.bsvp)

### extract predicted values
utc.index.pred.response.bsvp <- inla.stack.index(utc.join.stack.pred.15as.bsvp, "pred.response")$data 
utc.index.pred.latent.bsvp <- inla.stack.index(utc.join.stack.pred.15as.bsvp, "pred.latent")$data 

utc.pred.response.mean.bsvp <- utc.output.pred.15as.bsvp$summary.linear.predictor[utc.index.pred.response.bsvp, "mean"]
utc.pred.response.sd.bsvp <- utc.output.pred.15as.bsvp$summary.linear.predictor[utc.index.pred.response.bsvp, "sd"]

utc.pred.latent.mean.bsvp <- utc.output.pred.15as.bsvp$summary.linear.predictor[utc.index.pred.latent.bsvp, "mean"]
utc.pred.latent.sd.bsvp <- utc.output.pred.15as.bsvp$summary.linear.predictor[utc.index.pred.latent.bsvp, "sd"]

utc_pred_resp_mean.bsvp <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.bsvp.mod/4+1))
utc_pred_resp_mean_df.bsvp <- as.data.frame(mapply(c,utc_pred_resp_mean.bsvp, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.bsvp <- orderBy(~y + x, data=utc_pred_resp_mean_df.bsvp)
utc.pred.response.mean.grid.bsvp <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.bsvp$`log(utc.pred.response.mean.bsvp.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.bsvp <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.bsvp)
utc_pred_resp_sd_df.bsvp <- as.data.frame(mapply(c,utc_pred_resp_sd.bsvp, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.bsvp <- orderBy(~y + x, data=utc_pred_resp_sd_df.bsvp)
utc.pred.response.sd.grid.bsvp <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.bsvp$utc.pred.response.sd.bsvp)), nrow=utc.x.res*2)

utc_pred_lat_mean.bsvp <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.bsvp)
utc_pred_lat_mean_df.bsvp <- as.data.frame(mapply(c,utc_pred_lat_mean.bsvp, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.bsvp <- orderBy(~y + x, data=utc_pred_lat_mean_df.bsvp)
utc.pred.latent.mean.grid.bsvp <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.bsvp$utc.pred.latent.mean.bsvp)), nrow=utc.x.res*2)

utc_pred_lat_sd.bsvp <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.bsvp)
utc_pred_lat_sd_df.bsvp <- as.data.frame(mapply(c,utc_pred_lat_sd.bsvp, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.bsvp <- orderBy(~y + x, data=utc_pred_lat_sd_df.bsvp)
utc.pred.latent.sd.grid.bsvp <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.bsvp$utc.pred.latent.sd.bsvp)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.bsvp),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.bsvp), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.bsvp-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-1.5,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.bsvp-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.bsvp-utc.pred.latent.mean.grid.bsvp)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0.2,2.5)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.bsvp-utc.pred.latent.mean.grid.bsvp)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.bsvp,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.bsvp, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.bsvp,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.bsvp, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.bsvp.mod <- exp(utc.pred.response.mean.bsvp)-1
utc.pred.response.mean.bsvp.mod[utc.pred.response.mean.bsvp.mod<0]<-0
res.15as.bsvp<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.bsvp.mod/4
RMSE.15as.bsvp <- sqrt(mean(res.15as.bsvp^2))
RMSE.15as.bsvp 
correl.15as.bsvp <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.bsvp.mod/4) 
correl.15as.bsvp 

####################################
# best case scenario + worst prior
# full model w/ log transformation + normal prior est at 30 arcsec full data
# est at 3 arcsec full data + pred at 3 arcsec
####################################
utc.output.est.30as.log <- inla(utc.formula.log, 
                                data=inla.stack.data(utc.stack.est.30as.log, spde=utc.spde), family="gaussian", 
                                control.fixed=utc.fixed.vague, 
                                control.family=list(hyper=list(prec=utc.prec.vague)), 
                                control.predictor=list(A=inla.stack.A(utc.stack.est.30as.log), compute=TRUE),
                                control.compute=list(cpo=TRUE, dic=TRUE), verbose=TRUE)

summary(utc.output.est.30as.log)

### extract parameter estimates for model w/ ntl and slope and w/ log transformation
utc.30as.fixed.out <- round(utc.output.est.30as.log$summary.fixed[,1:5],3) 

# variance of unstructured residual
utc.30as.sigma2e_marg <- inla.tmarginal(function(x) 1/x,utc.output.est.30as.log$marginals.hyperpar[[1]])
utc.30as.sigma2e_m1 <- inla.emarginal(function(x) x, utc.30as.sigma2e_marg)
utc.30as.sigma2e_m2 <- inla.emarginal(function(x) x^2, utc.30as.sigma2e_marg)
utc.30as.sigma2e_stdev <- sqrt(utc.30as.sigma2e_m2 - utc.30as.sigma2e_m1^2)
utc.30as.sigma2e_quantiles <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.30as.sigma2e_marg)

utc.30as.spde.result <- inla.spde2.result(utc.output.est.30as.log, name="utc.spatial.field", spde=utc.spde)

# variance of spatial effects
utc.30as.var.nom.marg <- utc.30as.spde.result$marginals.variance.nominal[[1]]
utc.30as.var.nom.m1 <- inla.emarginal(function(x) x, utc.30as.var.nom.marg)
utc.30as.var.nom.m2 <- inla.emarginal(function(x) x^2, utc.30as.var.nom.marg)
utc.30as.var.nom.stdev <- sqrt(utc.30as.var.nom.m2 - utc.30as.var.nom.m1^2)
utc.30as.var.nom.quantiles <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.30as.var.nom.marg)

utc.30as.range.nom.marg <- utc.30as.spde.result$marginals.range.nominal[[1]]
utc.30as.range.nom.m1 <- inla.emarginal(function(x) x, utc.30as.range.nom.marg)
utc.30as.range.nom.m2 <- inla.emarginal(function(x) x^2, utc.30as.range.nom.marg)
utc.30as.range.nom.stdev <- sqrt(utc.30as.range.nom.m2 - utc.30as.range.nom.m1^2)
utc.30as.range.nom.quantiles <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.30as.range.nom.marg)

# set up priors
utc.30as.spde.normal <- inla.spde2.matern(mesh=utc.mesh, alpha=2, 
                                          B.tau=matrix(c(0, 1, 0),nrow=1,ncol=3),
                                          B.kappa=matrix(c(0, 0, 1),nrow=1,ncol=3), 
                                          theta.prior.mean=c(utc.30as.spde.result$summary.theta[[2]][1], 
                                                             utc.30as.spde.result$summary.theta[[2]][2]),
                                          theta.prior.prec=c(1/(utc.30as.spde.result$summary.theta[[3]][1]^2), 
                                                             1/(utc.30as.spde.result$summary.theta[[3]][2]^2)))

utc.30as.s.index.normal <- inla.spde.make.index(name="utc.30as.spatial.field.normal", n.spde=utc.30as.spde.normal$n.spde) 

utc.30as.spde.pc <- inla.spde2.pcmatern(utc.mesh, alpha = 2, 
                                        prior.range = c(utc.30as.range.nom.quantiles[1], 0.05),
                                        prior.sigma = c(sqrt(utc.30as.var.nom.quantiles[3]), 0.05))

utc.30as.s.index.pc <- inla.spde.make.index(name = "utc.30as.spatial.field.pc", n.spde = utc.30as.spde.pc$n.spde)

utc.30as.prec.normal <- list(initial=utc.output.est.30as.log$internal.summary.hyperpar[[1]][1], prior="normal", 
                             param=c(utc.output.est.30as.log$internal.summary.hyperpar[[1]][1], 
                                     1/(utc.output.est.30as.log$internal.summary.hyperpar[[2]][1]^2)), fixed=FALSE)

utc.30as.prec.pc <- list(prec=list(prior="pc.prec", param=c(sqrt(utc.30as.sigma2e_quantiles[3]), 0.05)))

utc.30as.fixed.normal <- list(mean.intercept = 0, 
                              prec.intercept = 0, 
                              mean = list(ntl=utc.output.est.30as.log$summary.fixed[[1]][2], 
                                          slope=utc.output.est.30as.log$summary.fixed[[1]][3]),
                              prec = list(ntl=1/(utc.output.est.30as.log$summary.fixed[[2]][2]^2),
                                          slope=1/(utc.output.est.30as.log$summary.fixed[[2]][3]^2)))

utc.stack.est.15as.bswp.n <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as)*4+1)),
                                        A = list(utc.A.pred.15as, 1),
                                        effects = list(c(utc.30as.s.index.normal, list(Intercept=1)), 
                                                       list(ntl=as.numeric(stacked_utc_pred_15as_fit$ntl_utc_15_15as), 
                                                            slope=as.numeric(stacked_utc_pred_15as_fit$slope_utc_15as)
                                                       )),
                                        tag="est")

utc.stack.pred.15as.wp.n <- inla.stack(data = list(logpop=NA),
                                       A = list(utc.A.pred.15as, 1),
                                       effects = list(c(utc.30as.s.index.normal, list(Intercept=1)),
                                                      list(ntl=as.numeric(stacked_utc_pred_15as_fit$ntl_utc_15_15as), 
                                                           slope=as.numeric(stacked_utc_pred_15as_fit$slope_utc_15as)                                                     
                                                      )),
                                       tag="pred.response") 

utc.stack.latent.wp.n <- inla.stack(data = list(xi=NA),
                                    A = list(utc.A.pred.15as),
                                    effects = list(utc.30as.s.index.normal),
                                    tag="pred.latent") 

utc.join.stack.pred.15as.bswp.n <- inla.stack(utc.stack.est.15as.bswp.n, utc.stack.latent.wp.n, utc.stack.pred.15as.wp.n)

utc.formula.log.30as.normal <- logpop ~ -1 + Intercept + ntl + slope + f(utc.30as.spatial.field.normal, model=utc.30as.spde.normal)

utc.output.pred.15as.bswp.n <- inla(utc.formula.log.30as.normal, 
                                    data=inla.stack.data(utc.join.stack.pred.15as.bswp.n, spde=utc.30as.spde.normal), family="gaussian",
                                    control.fixed=utc.30as.fixed.normal, 
                                    control.family=list(hyper=list(prec=utc.30as.prec.normal)), 
                                    control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.bswp.n), compute=TRUE),
                                    control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.bswp.n)

utc.output.pred.15as.bswp.n$dic$dic
utc.output.pred.15as.bswp.n$waic$waic
slcpo(utc.output.pred.15as.bswp.n)

### extract parameter estimates
utc.fixed.out.bswp.n <- round(utc.output.pred.15as.bswp.n$summary.fixed[,1:5],3) 
utc.fixed1_marg.bswp.n  <- utc.output.pred.15as.bswp.n$marginals.fixed[[1]]
utc.fixed2_marg.bswp.n  <- utc.output.pred.15as.bswp.n$marginals.fixed[[2]]
utc.fixed3_marg.bswp.n  <- utc.output.pred.15as.bswp.n$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.bswp.n <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.bswp.n$marginals.hyperpar[[1]])
utc.sigma2e_m1.bswp.n <- inla.emarginal(function(x) x, utc.sigma2e_marg.bswp.n)
utc.sigma2e_m2.bswp.n <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.bswp.n)
utc.sigma2e_stdev.bswp.n <- sqrt(utc.sigma2e_m2.bswp.n - utc.sigma2e_m1.bswp.n^2)
utc.sigma2e_quantiles.bswp.n <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.bswp.n)

utc.spde.result.bswp.n <- inla.spde2.result(utc.output.pred.15as.bswp.n, name="utc.30as.spatial.field.normal", utc.30as.spde.normal)

# variance of spatial effects
utc.var.nom.marg.bswp.n <- utc.spde.result.bswp.n$marginals.variance.nominal[[1]]
utc.var.nom.m1.bswp.n <- inla.emarginal(function(x) x, utc.var.nom.marg.bswp.n)
utc.var.nom.m2.bswp.n <- inla.emarginal(function(x) x^2, utc.var.nom.marg.bswp.n)
utc.var.nom.stdev.bswp.n <- sqrt(utc.var.nom.m2.bswp.n - utc.var.nom.m1.bswp.n^2)
utc.var.nom.quantiles.bswp.n <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.bswp.n)

utc.range.nom.marg.bswp.n <- utc.spde.result.bswp.n$marginals.range.nominal[[1]]
utc.range.nom.m1.bswp.n <- inla.emarginal(function(x) x, utc.range.nom.marg.bswp.n)
utc.range.nom.m2.bswp.n <- inla.emarginal(function(x) x^2, utc.range.nom.marg.bswp.n)
utc.range.nom.stdev.bswp.n <- sqrt(utc.range.nom.m2.bswp.n - utc.range.nom.m1.bswp.n^2)
utc.range.nom.quantiles.bswp.n <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.bswp.n)

### extract predicted values
utc.index.pred.response.bswp.n <- inla.stack.index(utc.join.stack.pred.15as.bswp.n, "pred.response")$data 
utc.index.pred.latent.bswp.n <- inla.stack.index(utc.join.stack.pred.15as.bswp.n, "pred.latent")$data 

utc.pred.response.mean.bswp.n <- utc.output.pred.15as.bswp.n$summary.linear.predictor[utc.index.pred.response.bswp.n, "mean"]
utc.pred.response.sd.bswp.n <- utc.output.pred.15as.bswp.n$summary.linear.predictor[utc.index.pred.response.bswp.n, "sd"]

utc.pred.latent.mean.bswp.n <- utc.output.pred.15as.bswp.n$summary.linear.predictor[utc.index.pred.latent.bswp.n, "mean"]
utc.pred.latent.sd.bswp.n <- utc.output.pred.15as.bswp.n$summary.linear.predictor[utc.index.pred.latent.bswp.n, "sd"]

utc_pred_resp_mean.bswp.n <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.bswp.n.mod/4+1))
utc_pred_resp_mean_df.bswp.n <- as.data.frame(mapply(c,utc_pred_resp_mean.bswp.n, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.bswp.n <- orderBy(~y + x, data=utc_pred_resp_mean_df.bswp.n)
utc.pred.response.mean.grid.bswp.n <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.bswp.n$`log(utc.pred.response.mean.bswp.n.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.bswp.n <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.bswp.n)
utc_pred_resp_sd_df.bswp.n <- as.data.frame(mapply(c,utc_pred_resp_sd.bswp.n, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.bswp.n <- orderBy(~y + x, data=utc_pred_resp_sd_df.bswp.n)
utc.pred.response.sd.grid.bswp.n <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.bswp.n$utc.pred.response.sd.bswp.n)), nrow=utc.x.res*2)

utc_pred_lat_mean.bswp.n <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.bswp.n)
utc_pred_lat_mean_df.bswp.n <- as.data.frame(mapply(c,utc_pred_lat_mean.bswp.n, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.bswp.n <- orderBy(~y + x, data=utc_pred_lat_mean_df.bswp.n)
utc.pred.latent.mean.grid.bswp.n <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.bswp.n$utc.pred.latent.mean.bswp.n)), nrow=utc.x.res*2)

utc_pred_lat_sd.bswp.n <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.bswp.n)
utc_pred_lat_sd_df.bswp.n <- as.data.frame(mapply(c,utc_pred_lat_sd.bswp.n, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.bswp.n <- orderBy(~y + x, data=utc_pred_lat_sd_df.bswp.n)
utc.pred.latent.sd.grid.bswp.n <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.bswp.n$utc.pred.latent.sd.bswp.n)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.bswp.n),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.bswp.n), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.bswp.n-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.bswp.n-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.bswp.n-utc.pred.latent.mean.grid.bswp.n)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.bswp.n-utc.pred.latent.mean.grid.bswp.n)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.bswp.n,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.bswp.n, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.bswp.n,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.bswp.n, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.bswp.n.mod <- exp(utc.pred.response.mean.bswp.n)-1
utc.pred.response.mean.bswp.n.mod[utc.pred.response.mean.bswp.n.mod<0]<-0
res.15as.bswp.n<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.bswp.n.mod/4
RMSE.15as.bswp.n <- sqrt(mean(res.15as.bswp.n^2))
RMSE.15as.bswp.n 
correl.15as.bswp.n <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.bswp.n.mod/4) 
correl.15as.bswp.n 

####################################
# best case scenario + worst prior
# full model w/ log transformation + pc prior est at 30 arcsec full data
# est at 3 arcsec full data + pred at 3 arcsec
####################################
utc.stack.est.15as.bswp.p <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as)*4+1)),
                                        A = list(utc.A.pred.15as, 1),
                                        effects = list(c(utc.30as.s.index.pc, list(Intercept=1)), 
                                                       list(ntl=as.numeric(stacked_utc_pred_15as_fit$ntl_utc_15_15as), 
                                                            slope=as.numeric(stacked_utc_pred_15as_fit$slope_utc_15as)
                                                       )),
                                        tag="est")

utc.stack.pred.15as.wp.p <- inla.stack(data = list(logpop=NA),
                                       A = list(utc.A.pred.15as, 1),
                                       effects = list(c(utc.30as.s.index.pc, list(Intercept=1)),
                                                      list(ntl=as.numeric(stacked_utc_pred_15as_fit$ntl_utc_15_15as), 
                                                           slope=as.numeric(stacked_utc_pred_15as_fit$slope_utc_15as)                                                     
                                                      )),
                                       tag="pred.response") 

utc.stack.latent.wp.p <- inla.stack(data = list(xi=NA),
                                    A = list(utc.A.pred.15as),
                                    effects = list(utc.30as.s.index.pc),
                                    tag="pred.latent") 

utc.join.stack.pred.15as.bswp.p <- inla.stack(utc.stack.est.15as.bswp.p, utc.stack.latent.wp.p, utc.stack.pred.15as.wp.p)

utc.formula.log.30as.pc <- logpop ~ -1 + Intercept + ntl + slope + f(utc.30as.spatial.field.pc, model=utc.30as.spde.pc)

utc.output.pred.15as.bswp.p <- inla(utc.formula.log.30as.pc, 
                                    data=inla.stack.data(utc.join.stack.pred.15as.bswp.p, spde=utc.30as.spde.pc), family="gaussian",
                                    control.fixed=utc.30as.fixed.normal, 
                                    control.family=list(hyper=list(prec=utc.30as.prec.pc)), 
                                    control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.bswp.p), compute=TRUE),
                                    control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.bswp.p)

utc.output.pred.15as.bswp.p$dic$dic
utc.output.pred.15as.bswp.p$waic$waic
slcpo(utc.output.pred.15as.bswp.p)

### extract parameter estimates
utc.fixed.out.bswp.p <- round(utc.output.pred.15as.bswp.p$summary.fixed[,1:5],3) 
utc.fixed1_marg.bswp.p  <- utc.output.pred.15as.bswp.p$marginals.fixed[[1]]
utc.fixed2_marg.bswp.p  <- utc.output.pred.15as.bswp.p$marginals.fixed[[2]]
utc.fixed3_marg.bswp.p  <- utc.output.pred.15as.bswp.p$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.bswp.p <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.bswp.p$marginals.hyperpar[[1]])
utc.sigma2e_m1.bswp.p <- inla.emarginal(function(x) x, utc.sigma2e_marg.bswp.p)
utc.sigma2e_m2.bswp.p <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.bswp.p)
utc.sigma2e_stdev.bswp.p <- sqrt(utc.sigma2e_m2.bswp.p - utc.sigma2e_m1.bswp.p^2)
utc.sigma2e_quantiles.bswp.p <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.bswp.p)

utc.spde.result.bswp.p <- inla.spde2.result(utc.output.pred.15as.bswp.p, name="utc.30as.spatial.field.pc", utc.30as.spde.pc)

# variance of spatial effects
utc.var.nom.marg.bswp.p <- utc.spde.result.bswp.p$marginals.variance.nominal[[1]]
utc.var.nom.m1.bswp.p <- inla.emarginal(function(x) x, utc.var.nom.marg.bswp.p)
utc.var.nom.m2.bswp.p <- inla.emarginal(function(x) x^2, utc.var.nom.marg.bswp.p)
utc.var.nom.stdev.bswp.p <- sqrt(utc.var.nom.m2.bswp.p - utc.var.nom.m1.bswp.p^2)
utc.var.nom.quantiles.bswp.p <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.bswp.p)

utc.range.nom.marg.bswp.p <- utc.spde.result.bswp.p$marginals.range.nominal[[1]]
utc.range.nom.m1.bswp.p <- inla.emarginal(function(x) x, utc.range.nom.marg.bswp.p)
utc.range.nom.m2.bswp.p <- inla.emarginal(function(x) x^2, utc.range.nom.marg.bswp.p)
utc.range.nom.stdev.bswp.p <- sqrt(utc.range.nom.m2.bswp.p - utc.range.nom.m1.bswp.p^2)
utc.range.nom.quantiles.bswp.p <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.bswp.p)

### extract predicted values
utc.index.pred.response.bswp.p <- inla.stack.index(utc.join.stack.pred.15as.bswp.p, "pred.response")$data 
utc.index.pred.latent.bswp.p <- inla.stack.index(utc.join.stack.pred.15as.bswp.p, "pred.latent")$data 

utc.pred.response.mean.bswp.p <- utc.output.pred.15as.bswp.p$summary.linear.predictor[utc.index.pred.response.bswp.p, "mean"]
utc.pred.response.sd.bswp.p <- utc.output.pred.15as.bswp.p$summary.linear.predictor[utc.index.pred.response.bswp.p, "sd"]

utc.pred.latent.mean.bswp.p <- utc.output.pred.15as.bswp.p$summary.linear.predictor[utc.index.pred.latent.bswp.p, "mean"]
utc.pred.latent.sd.bswp.p <- utc.output.pred.15as.bswp.p$summary.linear.predictor[utc.index.pred.latent.bswp.p, "sd"]

utc_pred_resp_mean.bswp.p <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.bswp.p.mod/4+1))
utc_pred_resp_mean_df.bswp.p <- as.data.frame(mapply(c,utc_pred_resp_mean.bswp.p, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.bswp.p <- orderBy(~y + x, data=utc_pred_resp_mean_df.bswp.p)
utc.pred.response.mean.grid.bswp.p <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.bswp.p$`log(utc.pred.response.mean.bswp.p.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.bswp.p <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.bswp.p)
utc_pred_resp_sd_df.bswp.p <- as.data.frame(mapply(c,utc_pred_resp_sd.bswp.p, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.bswp.p <- orderBy(~y + x, data=utc_pred_resp_sd_df.bswp.p)
utc.pred.response.sd.grid.bswp.p <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.bswp.p$utc.pred.response.sd.bswp.p)), nrow=utc.x.res*2)

utc_pred_lat_mean.bswp.p <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.bswp.p)
utc_pred_lat_mean_df.bswp.p <- as.data.frame(mapply(c,utc_pred_lat_mean.bswp.p, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.bswp.p <- orderBy(~y + x, data=utc_pred_lat_mean_df.bswp.p)
utc.pred.latent.mean.grid.bswp.p <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.bswp.p$utc.pred.latent.mean.bswp.p)), nrow=utc.x.res*2)

utc_pred_lat_sd.bswp.p <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.bswp.p)
utc_pred_lat_sd_df.bswp.p <- as.data.frame(mapply(c,utc_pred_lat_sd.bswp.p, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.bswp.p <- orderBy(~y + x, data=utc_pred_lat_sd_df.bswp.p)
utc.pred.latent.sd.grid.bswp.p <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.bswp.p$utc.pred.latent.sd.bswp.p)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.bswp.p),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.bswp.p), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.bswp.p-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.bswp.p-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.bswp.p-utc.pred.latent.mean.grid.bswp.p)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.bswp.p-utc.pred.latent.mean.grid.bswp.p)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.bswp.p,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.bswp.p, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.bswp.p,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.bswp.p, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.bswp.p.mod <- exp(utc.pred.response.mean.bswp.p)-1
utc.pred.response.mean.bswp.p.mod[utc.pred.response.mean.bswp.p.mod<0]<-0
res.15as.bswp.p<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.bswp.p.mod/4
RMSE.15as.bswp.p <- sqrt(mean(res.15as.bswp.p^2))
RMSE.15as.bswp.p 
correl.15as.bswp.p <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.bswp.p.mod/4) 
correl.15as.bswp.p 

####################################
# real case sampled scenario + worst prior
# full model w/ log transformation + normal prior est at 30 arcsec full data
# est at 3 arcsec sampled 1% data (not weighted) + pred at 3 arcsec
####################################
utc.loc_sampled_15as_u1 <- cbind(stacked_utc_sampled_u1$x, stacked_utc_sampled_u1$y)
utc.A.sampled.15as_u1 <- inla.spde.make.A(mesh=utc.mesh, loc=utc.loc_sampled_15as_u1)

utc.stack.est.15as.un1 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_u1$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_u1, 1),
                                     effects = list(c(utc.30as.s.index.normal, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_u1$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_u1$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.un1 <- inla.stack(utc.stack.est.15as.un1, utc.stack.latent.wp.n, utc.stack.pred.15as.wp.n)

utc.output.pred.15as.un1 <- inla(utc.formula.log.30as.normal, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.un1, spde=utc.30as.spde.normal), family="gaussian",
                                 control.fixed=utc.30as.fixed.normal, 
                                 control.family=list(hyper=list(prec=utc.30as.prec.normal)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.un1), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.un1)

utc.output.pred.15as.un1$dic$dic
utc.output.pred.15as.un1$waic$waic
slcpo(utc.output.pred.15as.un1)

### extract parameter estimates
utc.fixed.out.un1 <- round(utc.output.pred.15as.un1$summary.fixed[,1:5],3)
utc.fixed1_marg.un1  <- utc.output.pred.15as.un1$marginals.fixed[[1]]
utc.fixed2_marg.un1  <- utc.output.pred.15as.un1$marginals.fixed[[2]]
utc.fixed3_marg.un1  <- utc.output.pred.15as.un1$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.un1 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.un1$marginals.hyperpar[[1]])
utc.sigma2e_m1.un1 <- inla.emarginal(function(x) x, utc.sigma2e_marg.un1)
utc.sigma2e_m2.un1 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.un1)
utc.sigma2e_stdev.un1 <- sqrt(utc.sigma2e_m2.un1 - utc.sigma2e_m1.un1^2)
utc.sigma2e_quantiles.un1 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.un1)

utc.spde.result.un1 <- inla.spde2.result(utc.output.pred.15as.un1, name="utc.30as.spatial.field.normal", utc.30as.spde.normal)

# variance of spatial effects
utc.var.nom.marg.un1 <- utc.spde.result.un1$marginals.variance.nominal[[1]]
utc.var.nom.m1.un1 <- inla.emarginal(function(x) x, utc.var.nom.marg.un1)
utc.var.nom.m2.un1 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.un1)
utc.var.nom.stdev.un1 <- sqrt(utc.var.nom.m2.un1 - utc.var.nom.m1.un1^2)
utc.var.nom.quantiles.un1 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.un1)

utc.range.nom.marg.un1 <- utc.spde.result.un1$marginals.range.nominal[[1]]
utc.range.nom.m1.un1 <- inla.emarginal(function(x) x, utc.range.nom.marg.un1)
utc.range.nom.m2.un1 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.un1)
utc.range.nom.stdev.un1 <- sqrt(utc.range.nom.m2.un1 - utc.range.nom.m1.un1^2)
utc.range.nom.quantiles.un1 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.un1)

### extract predicted values
utc.index.pred.response.un1 <- inla.stack.index(utc.join.stack.pred.15as.un1, "pred.response")$data 
utc.index.pred.latent.un1 <- inla.stack.index(utc.join.stack.pred.15as.un1, "pred.latent")$data 

utc.pred.response.mean.un1 <- utc.output.pred.15as.un1$summary.linear.predictor[utc.index.pred.response.un1, "mean"]
utc.pred.response.sd.un1 <- utc.output.pred.15as.un1$summary.linear.predictor[utc.index.pred.response.un1, "sd"]

utc.pred.latent.mean.un1 <- utc.output.pred.15as.un1$summary.linear.predictor[utc.index.pred.latent.un1, "mean"]
utc.pred.latent.sd.un1 <- utc.output.pred.15as.un1$summary.linear.predictor[utc.index.pred.latent.un1, "sd"]

utc_pred_resp_mean.un1 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.un1.mod/4+1))
utc_pred_resp_mean_df.un1 <- as.data.frame(mapply(c,utc_pred_resp_mean.un1, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.un1 <- orderBy(~y + x, data=utc_pred_resp_mean_df.un1)
utc.pred.response.mean.grid.un1 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.un1$`log(utc.pred.response.mean.un1.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.un1 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.un1)
utc_pred_resp_sd_df.un1 <- as.data.frame(mapply(c,utc_pred_resp_sd.un1, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.un1 <- orderBy(~y + x, data=utc_pred_resp_sd_df.un1)
utc.pred.response.sd.grid.un1 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.un1$utc.pred.response.sd.un1)), nrow=utc.x.res*2)

utc_pred_lat_mean.un1 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.un1)
utc_pred_lat_mean_df.un1 <- as.data.frame(mapply(c,utc_pred_lat_mean.un1, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.un1 <- orderBy(~y + x, data=utc_pred_lat_mean_df.un1)
utc.pred.latent.mean.grid.un1 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.un1$utc.pred.latent.mean.un1)), nrow=utc.x.res*2)

utc_pred_lat_sd.un1 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.un1)
utc_pred_lat_sd_df.un1 <- as.data.frame(mapply(c,utc_pred_lat_sd.un1, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.un1 <- orderBy(~y + x, data=utc_pred_lat_sd_df.un1)
utc.pred.latent.sd.grid.un1 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.un1$utc.pred.latent.sd.un1)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.un1),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.un1), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.un1-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.un1-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.un1-utc.pred.latent.mean.grid.un1)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.un1-utc.pred.latent.mean.grid.un1)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.un1,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.un1, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.un1,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.un1, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.un1.mod <- exp(utc.pred.response.mean.un1)-1
utc.pred.response.mean.un1.mod[utc.pred.response.mean.un1.mod<0]<-0
res.15as.un1<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.un1.mod/4
RMSE.15as.un1 <- sqrt(mean(res.15as.un1^2))
RMSE.15as.un1 
correl.15as.un1 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.un1.mod/4) 
correl.15as.un1 

####################################
# real case sampled scenario + worst prior
# full model w/ log transformation + normal prior est at 30 arcsec full data
# est at 3 arcsec sampled 2% data (not weighted) + pred at 3 arcsec
####################################
utc.loc_sampled_15as_u2 <- cbind(stacked_utc_sampled_u2$x, stacked_utc_sampled_u2$y)
utc.A.sampled.15as_u2 <- inla.spde.make.A(mesh=utc.mesh, loc=utc.loc_sampled_15as_u2)

utc.stack.est.15as.un2 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_u2$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_u2, 1),
                                     effects = list(c(utc.30as.s.index.normal, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_u2$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_u2$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.un2 <- inla.stack(utc.stack.est.15as.un2, utc.stack.latent.wp.n, utc.stack.pred.15as.wp.n)

utc.output.pred.15as.un2 <- inla(utc.formula.log.30as.normal, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.un2, spde=utc.30as.spde.normal), family="gaussian",
                                 control.fixed=utc.30as.fixed.normal, 
                                 control.family=list(hyper=list(prec=utc.30as.prec.normal)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.un2), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.un2)

utc.output.pred.15as.un2$dic$dic
utc.output.pred.15as.un2$waic$waic
slcpo(utc.output.pred.15as.un2)

### extract parameter estimates
utc.fixed.out.un2 <- round(utc.output.pred.15as.un2$summary.fixed[,1:5],3) 
utc.fixed1_marg.un2  <- utc.output.pred.15as.un2$marginals.fixed[[1]]
utc.fixed2_marg.un2  <- utc.output.pred.15as.un2$marginals.fixed[[2]]
utc.fixed3_marg.un2  <- utc.output.pred.15as.un2$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.un2 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.un2$marginals.hyperpar[[1]])
utc.sigma2e_m1.un2 <- inla.emarginal(function(x) x, utc.sigma2e_marg.un2)
utc.sigma2e_m2.un2 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.un2)
utc.sigma2e_stdev.un2 <- sqrt(utc.sigma2e_m2.un2 - utc.sigma2e_m1.un2^2)
utc.sigma2e_quantiles.un2 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.un2)

utc.spde.result.un2 <- inla.spde2.result(utc.output.pred.15as.un2, name="utc.30as.spatial.field.normal", utc.30as.spde.normal)

# variance of spatial effects
utc.var.nom.marg.un2 <- utc.spde.result.un2$marginals.variance.nominal[[1]]
utc.var.nom.m1.un2 <- inla.emarginal(function(x) x, utc.var.nom.marg.un2)
utc.var.nom.m2.un2 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.un2)
utc.var.nom.stdev.un2 <- sqrt(utc.var.nom.m2.un2 - utc.var.nom.m1.un2^2)
utc.var.nom.quantiles.un2 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.un2)

utc.range.nom.marg.un2 <- utc.spde.result.un2$marginals.range.nominal[[1]]
utc.range.nom.m1.un2 <- inla.emarginal(function(x) x, utc.range.nom.marg.un2)
utc.range.nom.m2.un2 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.un2)
utc.range.nom.stdev.un2 <- sqrt(utc.range.nom.m2.un2 - utc.range.nom.m1.un2^2)
utc.range.nom.quantiles.un2 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.un2)

### extract predicted values
utc.index.pred.response.un2 <- inla.stack.index(utc.join.stack.pred.15as.un2, "pred.response")$data 
utc.index.pred.latent.un2 <- inla.stack.index(utc.join.stack.pred.15as.un2, "pred.latent")$data 

utc.pred.response.mean.un2 <- utc.output.pred.15as.un2$summary.linear.predictor[utc.index.pred.response.un2, "mean"]
utc.pred.response.sd.un2 <- utc.output.pred.15as.un2$summary.linear.predictor[utc.index.pred.response.un2, "sd"]

utc.pred.latent.mean.un2 <- utc.output.pred.15as.un2$summary.linear.predictor[utc.index.pred.latent.un2, "mean"]
utc.pred.latent.sd.un2 <- utc.output.pred.15as.un2$summary.linear.predictor[utc.index.pred.latent.un2, "sd"]

utc_pred_resp_mean.un2 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.un2.mod/4+1))
utc_pred_resp_mean_df.un2 <- as.data.frame(mapply(c,utc_pred_resp_mean.un2, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.un2 <- orderBy(~y + x, data=utc_pred_resp_mean_df.un2)
utc.pred.response.mean.grid.un2 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.un2$`log(utc.pred.response.mean.un2.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.un2 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.un2)
utc_pred_resp_sd_df.un2 <- as.data.frame(mapply(c,utc_pred_resp_sd.un2, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.un2 <- orderBy(~y + x, data=utc_pred_resp_sd_df.un2)
utc.pred.response.sd.grid.un2 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.un2$utc.pred.response.sd.un2)), nrow=utc.x.res*2)

utc_pred_lat_mean.un2 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.un2)
utc_pred_lat_mean_df.un2 <- as.data.frame(mapply(c,utc_pred_lat_mean.un2, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.un2 <- orderBy(~y + x, data=utc_pred_lat_mean_df.un2)
utc.pred.latent.mean.grid.un2 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.un2$utc.pred.latent.mean.un2)), nrow=utc.x.res*2)

utc_pred_lat_sd.un2 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.un2)
utc_pred_lat_sd_df.un2 <- as.data.frame(mapply(c,utc_pred_lat_sd.un2, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.un2 <- orderBy(~y + x, data=utc_pred_lat_sd_df.un2)
utc.pred.latent.sd.grid.un2 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.un2$utc.pred.latent.sd.un2)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.un2),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.un2), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.un2-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.un2-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.un2-utc.pred.latent.mean.grid.un2)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.un2-utc.pred.latent.mean.grid.un2)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.un2,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.un2, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.un2,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.un2, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.un2.mod <- exp(utc.pred.response.mean.un2)-1
utc.pred.response.mean.un2.mod[utc.pred.response.mean.un2.mod<0]<-0
res.15as.un2<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.un2.mod/4
RMSE.15as.un2 <- sqrt(mean(res.15as.un2^2))
RMSE.15as.un2 
correl.15as.un2 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.un2.mod/4) 
correl.15as.un2 

####################################
# real case sampled scenario + worst prior
# full model w/ log transformation + normal prior est at 30 arcsec full data
# est at 3 arcsec sampled 5% data (not weighted) + pred at 3 arcsec
####################################
utc.loc_sampled_15as_u3 <- cbind(stacked_utc_sampled_u3$x, stacked_utc_sampled_u3$y)
utc.A.sampled.15as_u3 <- inla.spde.make.A(mesh=utc.mesh, loc=utc.loc_sampled_15as_u3)

utc.stack.est.15as.un3 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_u3$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_u3, 1),
                                     effects = list(c(utc.30as.s.index.normal, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_u3$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_u3$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.un3 <- inla.stack(utc.stack.est.15as.un3, utc.stack.latent.wp.n, utc.stack.pred.15as.wp.n)

utc.output.pred.15as.un3 <- inla(utc.formula.log.30as.normal, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.un3, spde=utc.30as.spde.normal), family="gaussian",
                                 control.fixed=utc.30as.fixed.normal, 
                                 control.family=list(hyper=list(prec=utc.30as.prec.normal)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.un3), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.un3)

utc.output.pred.15as.un3$dic$dic
utc.output.pred.15as.un3$waic$waic
slcpo(utc.output.pred.15as.un3)

### extract parameter estimates
utc.fixed.out.un3 <- round(utc.output.pred.15as.un3$summary.fixed[,1:5],3) 
utc.fixed1_marg.un3  <- utc.output.pred.15as.un3$marginals.fixed[[1]]
utc.fixed2_marg.un3  <- utc.output.pred.15as.un3$marginals.fixed[[2]]
utc.fixed3_marg.un3  <- utc.output.pred.15as.un3$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.un3 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.un3$marginals.hyperpar[[1]])
utc.sigma2e_m1.un3 <- inla.emarginal(function(x) x, utc.sigma2e_marg.un3)
utc.sigma2e_m2.un3 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.un3)
utc.sigma2e_stdev.un3 <- sqrt(utc.sigma2e_m2.un3 - utc.sigma2e_m1.un3^2)
utc.sigma2e_quantiles.un3 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.un3)

utc.spde.result.un3 <- inla.spde2.result(utc.output.pred.15as.un3, name="utc.30as.spatial.field.normal", utc.30as.spde.normal)

# variance of spatial effects
utc.var.nom.marg.un3 <- utc.spde.result.un3$marginals.variance.nominal[[1]]
utc.var.nom.m1.un3 <- inla.emarginal(function(x) x, utc.var.nom.marg.un3)
utc.var.nom.m2.un3 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.un3)
utc.var.nom.stdev.un3 <- sqrt(utc.var.nom.m2.un3 - utc.var.nom.m1.un3^2)
utc.var.nom.quantiles.un3 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.un3)

utc.range.nom.marg.un3 <- utc.spde.result.un3$marginals.range.nominal[[1]]
utc.range.nom.m1.un3 <- inla.emarginal(function(x) x, utc.range.nom.marg.un3)
utc.range.nom.m2.un3 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.un3)
utc.range.nom.stdev.un3 <- sqrt(utc.range.nom.m2.un3 - utc.range.nom.m1.un3^2)
utc.range.nom.quantiles.un3 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.un3)

### extract predicted values
utc.index.pred.response.un3 <- inla.stack.index(utc.join.stack.pred.15as.un3, "pred.response")$data 
utc.index.pred.latent.un3 <- inla.stack.index(utc.join.stack.pred.15as.un3, "pred.latent")$data 

utc.pred.response.mean.un3 <- utc.output.pred.15as.un3$summary.linear.predictor[utc.index.pred.response.un3, "mean"]
utc.pred.response.sd.un3 <- utc.output.pred.15as.un3$summary.linear.predictor[utc.index.pred.response.un3, "sd"]

utc.pred.latent.mean.un3 <- utc.output.pred.15as.un3$summary.linear.predictor[utc.index.pred.latent.un3, "mean"]
utc.pred.latent.sd.un3 <- utc.output.pred.15as.un3$summary.linear.predictor[utc.index.pred.latent.un3, "sd"]

utc_pred_resp_mean.un3 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.un3.mod/4+1))
utc_pred_resp_mean_df.un3 <- as.data.frame(mapply(c,utc_pred_resp_mean.un3, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.un3 <- orderBy(~y + x, data=utc_pred_resp_mean_df.un3)
utc.pred.response.mean.grid.un3 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.un3$`log(utc.pred.response.mean.un3.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.un3 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.un3)
utc_pred_resp_sd_df.un3 <- as.data.frame(mapply(c,utc_pred_resp_sd.un3, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.un3 <- orderBy(~y + x, data=utc_pred_resp_sd_df.un3)
utc.pred.response.sd.grid.un3 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.un3$utc.pred.response.sd.un3)), nrow=utc.x.res*2)

utc_pred_lat_mean.un3 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.un3)
utc_pred_lat_mean_df.un3 <- as.data.frame(mapply(c,utc_pred_lat_mean.un3, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.un3 <- orderBy(~y + x, data=utc_pred_lat_mean_df.un3)
utc.pred.latent.mean.grid.un3 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.un3$utc.pred.latent.mean.un3)), nrow=utc.x.res*2)

utc_pred_lat_sd.un3 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.un3)
utc_pred_lat_sd_df.un3 <- as.data.frame(mapply(c,utc_pred_lat_sd.un3, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.un3 <- orderBy(~y + x, data=utc_pred_lat_sd_df.un3)
utc.pred.latent.sd.grid.un3 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.un3$utc.pred.latent.sd.un3)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.un3),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.un3), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.un3-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.un3-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.un3-utc.pred.latent.mean.grid.un3)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.un3-utc.pred.latent.mean.grid.un3)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.un3,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.un3, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.un3,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.un3, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.un3.mod <- exp(utc.pred.response.mean.un3)-1
utc.pred.response.mean.un3.mod[utc.pred.response.mean.un3.mod<0]<-0
res.15as.un3<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.un3.mod/4
RMSE.15as.un3 <- sqrt(mean(res.15as.un3^2))
RMSE.15as.un3 
correl.15as.un3 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.un3.mod/4) 
correl.15as.un3 

####################################
# real case sampled scenario + worst prior
# full model w/ log transformation + normal prior est at 30 arcsec full data
# est at 3 arcsec sampled 10% data (not weighted) + pred at 3 arcsec
####################################
utc.loc_sampled_15as_u4 <- cbind(stacked_utc_sampled_u4$x, stacked_utc_sampled_u4$y)
utc.A.sampled.15as_u4 <- inla.spde.make.A(mesh=utc.mesh, loc=utc.loc_sampled_15as_u4)

utc.stack.est.15as.un4 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_u4$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_u4, 1),
                                     effects = list(c(utc.30as.s.index.normal, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_u4$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_u4$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.un4 <- inla.stack(utc.stack.est.15as.un4, utc.stack.latent.wp.n, utc.stack.pred.15as.wp.n)

utc.output.pred.15as.un4 <- inla(utc.formula.log.30as.normal, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.un4, spde=utc.30as.spde.normal), family="gaussian",
                                 control.fixed=utc.30as.fixed.normal, 
                                 control.family=list(hyper=list(prec=utc.30as.prec.normal)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.un4), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.un4)

utc.output.pred.15as.un4$dic$dic
utc.output.pred.15as.un4$waic$waic
slcpo(utc.output.pred.15as.un4)

### extract parameter estimates
utc.fixed.out.un4 <- round(utc.output.pred.15as.un4$summary.fixed[,1:5],3) 
utc.fixed1_marg.un4  <- utc.output.pred.15as.un4$marginals.fixed[[1]]
utc.fixed2_marg.un4  <- utc.output.pred.15as.un4$marginals.fixed[[2]]
utc.fixed3_marg.un4  <- utc.output.pred.15as.un4$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.un4 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.un4$marginals.hyperpar[[1]])
utc.sigma2e_m1.un4 <- inla.emarginal(function(x) x, utc.sigma2e_marg.un4)
utc.sigma2e_m2.un4 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.un4)
utc.sigma2e_stdev.un4 <- sqrt(utc.sigma2e_m2.un4 - utc.sigma2e_m1.un4^2)
utc.sigma2e_quantiles.un4 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.un4)

utc.spde.result.un4 <- inla.spde2.result(utc.output.pred.15as.un4, name="utc.30as.spatial.field.normal", utc.30as.spde.normal)

# variance of spatial effects
utc.var.nom.marg.un4 <- utc.spde.result.un4$marginals.variance.nominal[[1]]
utc.var.nom.m1.un4 <- inla.emarginal(function(x) x, utc.var.nom.marg.un4)
utc.var.nom.m2.un4 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.un4)
utc.var.nom.stdev.un4 <- sqrt(utc.var.nom.m2.un4 - utc.var.nom.m1.un4^2)
utc.var.nom.quantiles.un4 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.un4)

utc.range.nom.marg.un4 <- utc.spde.result.un4$marginals.range.nominal[[1]]
utc.range.nom.m1.un4 <- inla.emarginal(function(x) x, utc.range.nom.marg.un4)
utc.range.nom.m2.un4 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.un4)
utc.range.nom.stdev.un4 <- sqrt(utc.range.nom.m2.un4 - utc.range.nom.m1.un4^2)
utc.range.nom.quantiles.un4 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.un4)

### extract predicted values
utc.index.pred.response.un4 <- inla.stack.index(utc.join.stack.pred.15as.un4, "pred.response")$data 
utc.index.pred.latent.un4 <- inla.stack.index(utc.join.stack.pred.15as.un4, "pred.latent")$data 

utc.pred.response.mean.un4 <- utc.output.pred.15as.un4$summary.linear.predictor[utc.index.pred.response.un4, "mean"]
utc.pred.response.sd.un4 <- utc.output.pred.15as.un4$summary.linear.predictor[utc.index.pred.response.un4, "sd"]

utc.pred.latent.mean.un4 <- utc.output.pred.15as.un4$summary.linear.predictor[utc.index.pred.latent.un4, "mean"]
utc.pred.latent.sd.un4 <- utc.output.pred.15as.un4$summary.linear.predictor[utc.index.pred.latent.un4, "sd"]

utc_pred_resp_mean.un4 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.un4.mod/4+1))
utc_pred_resp_mean_df.un4 <- as.data.frame(mapply(c,utc_pred_resp_mean.un4, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.un4 <- orderBy(~y + x, data=utc_pred_resp_mean_df.un4)
utc.pred.response.mean.grid.un4 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.un4$`log(utc.pred.response.mean.un4.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.un4 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.un4)
utc_pred_resp_sd_df.un4 <- as.data.frame(mapply(c,utc_pred_resp_sd.un4, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.un4 <- orderBy(~y + x, data=utc_pred_resp_sd_df.un4)
utc.pred.response.sd.grid.un4 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.un4$utc.pred.response.sd.un4)), nrow=utc.x.res*2)

utc_pred_lat_mean.un4 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.un4)
utc_pred_lat_mean_df.un4 <- as.data.frame(mapply(c,utc_pred_lat_mean.un4, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.un4 <- orderBy(~y + x, data=utc_pred_lat_mean_df.un4)
utc.pred.latent.mean.grid.un4 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.un4$utc.pred.latent.mean.un4)), nrow=utc.x.res*2)

utc_pred_lat_sd.un4 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.un4)
utc_pred_lat_sd_df.un4 <- as.data.frame(mapply(c,utc_pred_lat_sd.un4, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.un4 <- orderBy(~y + x, data=utc_pred_lat_sd_df.un4)
utc.pred.latent.sd.grid.un4 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.un4$utc.pred.latent.sd.un4)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.un4),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.un4), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.un4-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.un4-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.un4-utc.pred.latent.mean.grid.un4)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.un4-utc.pred.latent.mean.grid.un4)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.un4,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.un4, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.un4,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.un4, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.un4.mod <- exp(utc.pred.response.mean.un4)-1
utc.pred.response.mean.un4.mod[utc.pred.response.mean.un4.mod<0]<-0
res.15as.un4<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.un4.mod/4
RMSE.15as.un4 <- sqrt(mean(res.15as.un4^2))
RMSE.15as.un4 
correl.15as.un4 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.un4.mod/4) 
correl.15as.un4 

####################################
# real case sampled scenario + worst prior
# full model w/ log transformation + normal prior est at 30 arcsec full data
# est at 3 arcsec sampled 20% data (not weighted) + pred at 3 arcsec
####################################
utc.loc_sampled_15as_u5 <- cbind(stacked_utc_sampled_u5$x, stacked_utc_sampled_u5$y)
utc.A.sampled.15as_u5 <- inla.spde.make.A(mesh=utc.mesh, loc=utc.loc_sampled_15as_u5)

utc.stack.est.15as.un5 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_u5$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_u5, 1),
                                     effects = list(c(utc.30as.s.index.normal, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_u5$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_u5$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.un5 <- inla.stack(utc.stack.est.15as.un5, utc.stack.latent.wp.n, utc.stack.pred.15as.wp.n)

utc.output.pred.15as.un5 <- inla(utc.formula.log.30as.normal, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.un5, spde=utc.30as.spde.normal), family="gaussian",
                                 control.fixed=utc.30as.fixed.normal, 
                                 control.family=list(hyper=list(prec=utc.30as.prec.normal)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.un5), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.un5)

utc.output.pred.15as.un5$dic$dic
utc.output.pred.15as.un5$waic$waic
slcpo(utc.output.pred.15as.un5)

### extract parameter estimates
utc.fixed.out.un5 <- round(utc.output.pred.15as.un5$summary.fixed[,1:5],3) 
utc.fixed1_marg.un5  <- utc.output.pred.15as.un5$marginals.fixed[[1]]
utc.fixed2_marg.un5  <- utc.output.pred.15as.un5$marginals.fixed[[2]]
utc.fixed3_marg.un5  <- utc.output.pred.15as.un5$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.un5 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.un5$marginals.hyperpar[[1]])
utc.sigma2e_m1.un5 <- inla.emarginal(function(x) x, utc.sigma2e_marg.un5)
utc.sigma2e_m2.un5 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.un5)
utc.sigma2e_stdev.un5 <- sqrt(utc.sigma2e_m2.un5 - utc.sigma2e_m1.un5^2)
utc.sigma2e_quantiles.un5 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.un5)

utc.spde.result.un5 <- inla.spde2.result(utc.output.pred.15as.un5, name="utc.30as.spatial.field.normal", utc.30as.spde.normal)

# variance of spatial effects
utc.var.nom.marg.un5 <- utc.spde.result.un5$marginals.variance.nominal[[1]]
utc.var.nom.m1.un5 <- inla.emarginal(function(x) x, utc.var.nom.marg.un5)
utc.var.nom.m2.un5 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.un5)
utc.var.nom.stdev.un5 <- sqrt(utc.var.nom.m2.un5 - utc.var.nom.m1.un5^2)
utc.var.nom.quantiles.un5 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.un5)

utc.range.nom.marg.un5 <- utc.spde.result.un5$marginals.range.nominal[[1]]
utc.range.nom.m1.un5 <- inla.emarginal(function(x) x, utc.range.nom.marg.un5)
utc.range.nom.m2.un5 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.un5)
utc.range.nom.stdev.un5 <- sqrt(utc.range.nom.m2.un5 - utc.range.nom.m1.un5^2)
utc.range.nom.quantiles.un5 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.un5)

### extract predicted values
utc.index.pred.response.un5 <- inla.stack.index(utc.join.stack.pred.15as.un5, "pred.response")$data 
utc.index.pred.latent.un5 <- inla.stack.index(utc.join.stack.pred.15as.un5, "pred.latent")$data 

utc.pred.response.mean.un5 <- utc.output.pred.15as.un5$summary.linear.predictor[utc.index.pred.response.un5, "mean"]
utc.pred.response.sd.un5 <- utc.output.pred.15as.un5$summary.linear.predictor[utc.index.pred.response.un5, "sd"]

utc.pred.latent.mean.un5 <- utc.output.pred.15as.un5$summary.linear.predictor[utc.index.pred.latent.un5, "mean"]
utc.pred.latent.sd.un5 <- utc.output.pred.15as.un5$summary.linear.predictor[utc.index.pred.latent.un5, "sd"]

utc_pred_resp_mean.un5 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.un5.mod/4+1))
utc_pred_resp_mean_df.un5 <- as.data.frame(mapply(c,utc_pred_resp_mean.un5, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.un5 <- orderBy(~y + x, data=utc_pred_resp_mean_df.un5)
utc.pred.response.mean.grid.un5 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.un5$`log(utc.pred.response.mean.un5.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.un5 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.un5)
utc_pred_resp_sd_df.un5 <- as.data.frame(mapply(c,utc_pred_resp_sd.un5, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.un5 <- orderBy(~y + x, data=utc_pred_resp_sd_df.un5)
utc.pred.response.sd.grid.un5 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.un5$utc.pred.response.sd.un5)), nrow=utc.x.res*2)

utc_pred_lat_mean.un5 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.un5)
utc_pred_lat_mean_df.un5 <- as.data.frame(mapply(c,utc_pred_lat_mean.un5, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.un5 <- orderBy(~y + x, data=utc_pred_lat_mean_df.un5)
utc.pred.latent.mean.grid.un5 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.un5$utc.pred.latent.mean.un5)), nrow=utc.x.res*2)

utc_pred_lat_sd.un5 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.un5)
utc_pred_lat_sd_df.un5 <- as.data.frame(mapply(c,utc_pred_lat_sd.un5, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.un5 <- orderBy(~y + x, data=utc_pred_lat_sd_df.un5)
utc.pred.latent.sd.grid.un5 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.un5$utc.pred.latent.sd.un5)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.un5),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.un5), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.un5-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.un5-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.un5-utc.pred.latent.mean.grid.un5)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.un5-utc.pred.latent.mean.grid.un5)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.un5,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.un5, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.un5,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.un5, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.un5.mod <- exp(utc.pred.response.mean.un5)-1
utc.pred.response.mean.un5.mod[utc.pred.response.mean.un5.mod<0]<-0
res.15as.un5<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.un5.mod/4
RMSE.15as.un5 <- sqrt(mean(res.15as.un5^2))
RMSE.15as.un5 
correl.15as.un5 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.un5.mod/4) 
correl.15as.un5 

####################################
# real case sampled scenario + worst prior
# full model w/ log transformation + normal prior est at 30 arcsec full data
# est at 3 arcsec sampled 50% data (not weighted) + pred at 3 arcsec
####################################
utc.loc_sampled_15as_u6 <- cbind(stacked_utc_sampled_u6$x, stacked_utc_sampled_u6$y)
utc.A.sampled.15as_u6 <- inla.spde.make.A(mesh=utc.mesh, loc=utc.loc_sampled_15as_u6)

utc.stack.est.15as.un6 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_u6$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_u6, 1),
                                     effects = list(c(utc.30as.s.index.normal, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_u6$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_u6$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.un6 <- inla.stack(utc.stack.est.15as.un6, utc.stack.latent.wp.n, utc.stack.pred.15as.wp.n)

utc.output.pred.15as.un6 <- inla(utc.formula.log.30as.normal, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.un6, spde=utc.30as.spde.normal), family="gaussian",
                                 control.fixed=utc.30as.fixed.normal, 
                                 control.family=list(hyper=list(prec=utc.30as.prec.normal)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.un6), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.un6)

utc.output.pred.15as.un6$dic$dic
utc.output.pred.15as.un6$waic$waic
slcpo(utc.output.pred.15as.un6)

### extract parameter estimates
utc.fixed.out.un6 <- round(utc.output.pred.15as.un6$summary.fixed[,1:5],3) 
utc.fixed1_marg.un6  <- utc.output.pred.15as.un6$marginals.fixed[[1]]
utc.fixed2_marg.un6  <- utc.output.pred.15as.un6$marginals.fixed[[2]]
utc.fixed3_marg.un6  <- utc.output.pred.15as.un6$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.un6 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.un6$marginals.hyperpar[[1]])
utc.sigma2e_m1.un6 <- inla.emarginal(function(x) x, utc.sigma2e_marg.un6)
utc.sigma2e_m2.un6 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.un6)
utc.sigma2e_stdev.un6 <- sqrt(utc.sigma2e_m2.un6 - utc.sigma2e_m1.un6^2)
utc.sigma2e_quantiles.un6 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.un6)

utc.spde.result.un6 <- inla.spde2.result(utc.output.pred.15as.un6, name="utc.30as.spatial.field.normal", utc.30as.spde.normal)

# variance of spatial effects
utc.var.nom.marg.un6 <- utc.spde.result.un6$marginals.variance.nominal[[1]]
utc.var.nom.m1.un6 <- inla.emarginal(function(x) x, utc.var.nom.marg.un6)
utc.var.nom.m2.un6 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.un6)
utc.var.nom.stdev.un6 <- sqrt(utc.var.nom.m2.un6 - utc.var.nom.m1.un6^2)
utc.var.nom.quantiles.un6 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.un6)

utc.range.nom.marg.un6 <- utc.spde.result.un6$marginals.range.nominal[[1]]
utc.range.nom.m1.un6 <- inla.emarginal(function(x) x, utc.range.nom.marg.un6)
utc.range.nom.m2.un6 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.un6)
utc.range.nom.stdev.un6 <- sqrt(utc.range.nom.m2.un6 - utc.range.nom.m1.un6^2)
utc.range.nom.quantiles.un6 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.un6)

### extract predicted values
utc.index.pred.response.un6 <- inla.stack.index(utc.join.stack.pred.15as.un6, "pred.response")$data 
utc.index.pred.latent.un6 <- inla.stack.index(utc.join.stack.pred.15as.un6, "pred.latent")$data 

utc.pred.response.mean.un6 <- utc.output.pred.15as.un6$summary.linear.predictor[utc.index.pred.response.un6, "mean"]
utc.pred.response.sd.un6 <- utc.output.pred.15as.un6$summary.linear.predictor[utc.index.pred.response.un6, "sd"]

utc.pred.latent.mean.un6 <- utc.output.pred.15as.un6$summary.linear.predictor[utc.index.pred.latent.un6, "mean"]
utc.pred.latent.sd.un6 <- utc.output.pred.15as.un6$summary.linear.predictor[utc.index.pred.latent.un6, "sd"]

utc_pred_resp_mean.un6 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.un6.mod/4+1))
utc_pred_resp_mean_df.un6 <- as.data.frame(mapply(c,utc_pred_resp_mean.un6, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.un6 <- orderBy(~y + x, data=utc_pred_resp_mean_df.un6)
utc.pred.response.mean.grid.un6 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.un6$`log(utc.pred.response.mean.un6.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.un6 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.un6)
utc_pred_resp_sd_df.un6 <- as.data.frame(mapply(c,utc_pred_resp_sd.un6, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.un6 <- orderBy(~y + x, data=utc_pred_resp_sd_df.un6)
utc.pred.response.sd.grid.un6 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.un6$utc.pred.response.sd.un6)), nrow=utc.x.res*2)

utc_pred_lat_mean.un6 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.un6)
utc_pred_lat_mean_df.un6 <- as.data.frame(mapply(c,utc_pred_lat_mean.un6, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.un6 <- orderBy(~y + x, data=utc_pred_lat_mean_df.un6)
utc.pred.latent.mean.grid.un6 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.un6$utc.pred.latent.mean.un6)), nrow=utc.x.res*2)

utc_pred_lat_sd.un6 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.un6)
utc_pred_lat_sd_df.un6 <- as.data.frame(mapply(c,utc_pred_lat_sd.un6, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.un6 <- orderBy(~y + x, data=utc_pred_lat_sd_df.un6)
utc.pred.latent.sd.grid.un6 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.un6$utc.pred.latent.sd.un6)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.un6),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.un6), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.un6-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.un6-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.un6-utc.pred.latent.mean.grid.un6)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.un6-utc.pred.latent.mean.grid.un6)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.un6,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.un6, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.un6,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.un6, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.un6.mod <- exp(utc.pred.response.mean.un6)-1
utc.pred.response.mean.un6.mod[utc.pred.response.mean.un6.mod<0]<-0
res.15as.un6<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.un6.mod/4
RMSE.15as.un6 <- sqrt(mean(res.15as.un6^2))
RMSE.15as.un6 
correl.15as.un6 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.un6.mod/4) 
correl.15as.un6 

####################################
# real case sampled scenario + worst prior
# full model w/ log transformation + pc prior est at 30 arcsec full data
# est at 3 arcsec sampled 1% data (not weighted) + pred at 3 arcsec
####################################
utc.stack.est.15as.up1 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_u1$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_u1, 1),
                                     effects = list(c(utc.30as.s.index.pc, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_u1$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_u1$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.up1 <- inla.stack(utc.stack.est.15as.up1, utc.stack.latent.wp.p, utc.stack.pred.15as.wp.p)

utc.output.pred.15as.up1 <- inla(utc.formula.log.30as.pc, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.up1, spde=utc.30as.spde.pc), family="gaussian",
                                 control.fixed=utc.30as.fixed.normal, 
                                 control.family=list(hyper=list(prec=utc.30as.prec.pc)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.up1), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.up1)

utc.output.pred.15as.up1$dic$dic
utc.output.pred.15as.up1$waic$waic
slcpo(utc.output.pred.15as.up1)

### extract parameter estimates
utc.fixed.out.up1 <- round(utc.output.pred.15as.up1$summary.fixed[,1:5],3) 
utc.fixed1_marg.up1  <- utc.output.pred.15as.up1$marginals.fixed[[1]]
utc.fixed2_marg.up1  <- utc.output.pred.15as.up1$marginals.fixed[[2]]
utc.fixed3_marg.up1  <- utc.output.pred.15as.up1$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.up1 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.up1$marginals.hyperpar[[1]])
utc.sigma2e_m1.up1 <- inla.emarginal(function(x) x, utc.sigma2e_marg.up1)
utc.sigma2e_m2.up1 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.up1)
utc.sigma2e_stdev.up1 <- sqrt(utc.sigma2e_m2.up1 - utc.sigma2e_m1.up1^2)
utc.sigma2e_quantiles.up1 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.up1)

utc.spde.result.up1 <- inla.spde2.result(utc.output.pred.15as.up1, name="utc.30as.spatial.field.pc", utc.30as.spde.pc)

# variance of spatial effects
utc.var.nom.marg.up1 <- utc.spde.result.up1$marginals.variance.nominal[[1]]
utc.var.nom.m1.up1 <- inla.emarginal(function(x) x, utc.var.nom.marg.up1)
utc.var.nom.m2.up1 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.up1)
utc.var.nom.stdev.up1 <- sqrt(utc.var.nom.m2.up1 - utc.var.nom.m1.up1^2)
utc.var.nom.quantiles.up1 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.up1)

utc.range.nom.marg.up1 <- utc.spde.result.up1$marginals.range.nominal[[1]]
utc.range.nom.m1.up1 <- inla.emarginal(function(x) x, utc.range.nom.marg.up1)
utc.range.nom.m2.up1 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.up1)
utc.range.nom.stdev.up1 <- sqrt(utc.range.nom.m2.up1 - utc.range.nom.m1.up1^2)
utc.range.nom.quantiles.up1 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.up1)

### extract predicted values
utc.index.pred.response.up1 <- inla.stack.index(utc.join.stack.pred.15as.up1, "pred.response")$data 
utc.index.pred.latent.up1 <- inla.stack.index(utc.join.stack.pred.15as.up1, "pred.latent")$data 

utc.pred.response.mean.up1 <- utc.output.pred.15as.up1$summary.linear.predictor[utc.index.pred.response.up1, "mean"]
utc.pred.response.sd.up1 <- utc.output.pred.15as.up1$summary.linear.predictor[utc.index.pred.response.up1, "sd"]

utc.pred.latent.mean.up1 <- utc.output.pred.15as.up1$summary.linear.predictor[utc.index.pred.latent.up1, "mean"]
utc.pred.latent.sd.up1 <- utc.output.pred.15as.up1$summary.linear.predictor[utc.index.pred.latent.up1, "sd"]

utc_pred_resp_mean.up1 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.up1.mod/4+1))
utc_pred_resp_mean_df.up1 <- as.data.frame(mapply(c,utc_pred_resp_mean.up1, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.up1 <- orderBy(~y + x, data=utc_pred_resp_mean_df.up1)
utc.pred.response.mean.grid.up1 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.up1$`log(utc.pred.response.mean.up1.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.up1 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.up1)
utc_pred_resp_sd_df.up1 <- as.data.frame(mapply(c,utc_pred_resp_sd.up1, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.up1 <- orderBy(~y + x, data=utc_pred_resp_sd_df.up1)
utc.pred.response.sd.grid.up1 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.up1$utc.pred.response.sd.up1)), nrow=utc.x.res*2)

utc_pred_lat_mean.up1 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.up1)
utc_pred_lat_mean_df.up1 <- as.data.frame(mapply(c,utc_pred_lat_mean.up1, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.up1 <- orderBy(~y + x, data=utc_pred_lat_mean_df.up1)
utc.pred.latent.mean.grid.up1 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.up1$utc.pred.latent.mean.up1)), nrow=utc.x.res*2)

utc_pred_lat_sd.up1 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.up1)
utc_pred_lat_sd_df.up1 <- as.data.frame(mapply(c,utc_pred_lat_sd.up1, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.up1 <- orderBy(~y + x, data=utc_pred_lat_sd_df.up1)
utc.pred.latent.sd.grid.up1 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.up1$utc.pred.latent.sd.up1)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.up1),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.up1), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.up1-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.up1-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.up1-utc.pred.latent.mean.grid.up1)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.up1-utc.pred.latent.mean.grid.up1)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.up1,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.up1, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.up1,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.up1, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.up1.mod <- exp(utc.pred.response.mean.up1)-1
utc.pred.response.mean.up1.mod[utc.pred.response.mean.up1.mod<0]<-0
res.15as.up1<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.up1.mod/4
RMSE.15as.up1 <- sqrt(mean(res.15as.up1^2))
RMSE.15as.up1 
correl.15as.up1 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.up1.mod/4) 
correl.15as.up1 

####################################
# real case sampled scenario + worst prior
# full model w/ log transformation + pc prior est at 30 arcsec full data
# est at 3 arcsec sampled 2% data (not weighted) + pred at 3 arcsec
####################################
utc.stack.est.15as.up2 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_u2$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_u2, 1),
                                     effects = list(c(utc.30as.s.index.pc, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_u2$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_u2$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.up2 <- inla.stack(utc.stack.est.15as.up2, utc.stack.latent.wp.p, utc.stack.pred.15as.wp.p)

utc.output.pred.15as.up2 <- inla(utc.formula.log.30as.pc, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.up2, spde=utc.30as.spde.pc), family="gaussian",
                                 control.fixed=utc.30as.fixed.normal, 
                                 control.family=list(hyper=list(prec=utc.30as.prec.pc)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.up2), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.up2)

utc.output.pred.15as.up2$dic$dic
utc.output.pred.15as.up2$waic$waic
slcpo(utc.output.pred.15as.up2)

### extract parameter estimates
utc.fixed.out.up2 <- round(utc.output.pred.15as.up2$summary.fixed[,1:5],3)
utc.fixed1_marg.up2  <- utc.output.pred.15as.up2$marginals.fixed[[1]]
utc.fixed2_marg.up2  <- utc.output.pred.15as.up2$marginals.fixed[[2]]
utc.fixed3_marg.up2  <- utc.output.pred.15as.up2$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.up2 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.up2$marginals.hyperpar[[1]])
utc.sigma2e_m1.up2 <- inla.emarginal(function(x) x, utc.sigma2e_marg.up2)
utc.sigma2e_m2.up2 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.up2)
utc.sigma2e_stdev.up2 <- sqrt(utc.sigma2e_m2.up2 - utc.sigma2e_m1.up2^2)
utc.sigma2e_quantiles.up2 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.up2)

utc.spde.result.up2 <- inla.spde2.result(utc.output.pred.15as.up2, name="utc.30as.spatial.field.pc", utc.30as.spde.pc)

# variance of spatial effects
utc.var.nom.marg.up2 <- utc.spde.result.up2$marginals.variance.nominal[[1]]
utc.var.nom.m1.up2 <- inla.emarginal(function(x) x, utc.var.nom.marg.up2)
utc.var.nom.m2.up2 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.up2)
utc.var.nom.stdev.up2 <- sqrt(utc.var.nom.m2.up2 - utc.var.nom.m1.up2^2)
utc.var.nom.quantiles.up2 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.up2)

utc.range.nom.marg.up2 <- utc.spde.result.up2$marginals.range.nominal[[1]]
utc.range.nom.m1.up2 <- inla.emarginal(function(x) x, utc.range.nom.marg.up2)
utc.range.nom.m2.up2 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.up2)
utc.range.nom.stdev.up2 <- sqrt(utc.range.nom.m2.up2 - utc.range.nom.m1.up2^2)
utc.range.nom.quantiles.up2 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.up2)

### extract predicted values
utc.index.pred.response.up2 <- inla.stack.index(utc.join.stack.pred.15as.up2, "pred.response")$data 
utc.index.pred.latent.up2 <- inla.stack.index(utc.join.stack.pred.15as.up2, "pred.latent")$data 

utc.pred.response.mean.up2 <- utc.output.pred.15as.up2$summary.linear.predictor[utc.index.pred.response.up2, "mean"]
utc.pred.response.sd.up2 <- utc.output.pred.15as.up2$summary.linear.predictor[utc.index.pred.response.up2, "sd"]

utc.pred.latent.mean.up2 <- utc.output.pred.15as.up2$summary.linear.predictor[utc.index.pred.latent.up2, "mean"]
utc.pred.latent.sd.up2 <- utc.output.pred.15as.up2$summary.linear.predictor[utc.index.pred.latent.up2, "sd"]

utc_pred_resp_mean.up2 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.up2.mod/4+1))
utc_pred_resp_mean_df.up2 <- as.data.frame(mapply(c,utc_pred_resp_mean.up2, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.up2 <- orderBy(~y + x, data=utc_pred_resp_mean_df.up2)
utc.pred.response.mean.grid.up2 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.up2$`log(utc.pred.response.mean.up2.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.up2 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.up2)
utc_pred_resp_sd_df.up2 <- as.data.frame(mapply(c,utc_pred_resp_sd.up2, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.up2 <- orderBy(~y + x, data=utc_pred_resp_sd_df.up2)
utc.pred.response.sd.grid.up2 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.up2$utc.pred.response.sd.up2)), nrow=utc.x.res*2)

utc_pred_lat_mean.up2 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.up2)
utc_pred_lat_mean_df.up2 <- as.data.frame(mapply(c,utc_pred_lat_mean.up2, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.up2 <- orderBy(~y + x, data=utc_pred_lat_mean_df.up2)
utc.pred.latent.mean.grid.up2 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.up2$utc.pred.latent.mean.up2)), nrow=utc.x.res*2)

utc_pred_lat_sd.up2 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.up2)
utc_pred_lat_sd_df.up2 <- as.data.frame(mapply(c,utc_pred_lat_sd.up2, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.up2 <- orderBy(~y + x, data=utc_pred_lat_sd_df.up2)
utc.pred.latent.sd.grid.up2 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.up2$utc.pred.latent.sd.up2)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.up2),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.up2), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.up2-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.up2-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.up2-utc.pred.latent.mean.grid.up2)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.up2-utc.pred.latent.mean.grid.up2)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.up2,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.up2, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.up2,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.up2, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.up2.mod <- exp(utc.pred.response.mean.up2)-1
utc.pred.response.mean.up2.mod[utc.pred.response.mean.up2.mod<0]<-0
res.15as.up2<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.up2.mod/4
RMSE.15as.up2 <- sqrt(mean(res.15as.up2^2))
RMSE.15as.up2 
correl.15as.up2 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.up2.mod/4) 
correl.15as.up2 

####################################
# real case sampled scenario + worst prior
# full model w/ log transformation + pc prior est at 30 arcsec full data
# est at 3 arcsec sampled 5% data (not weighted) + pred at 3 arcsec
####################################
utc.stack.est.15as.up3 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_u3$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_u3, 1),
                                     effects = list(c(utc.30as.s.index.pc, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_u3$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_u3$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.up3 <- inla.stack(utc.stack.est.15as.up3, utc.stack.latent.wp.p, utc.stack.pred.15as.wp.p)

utc.output.pred.15as.up3 <- inla(utc.formula.log.30as.pc, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.up3, spde=utc.30as.spde.pc), family="gaussian",
                                 control.fixed=utc.30as.fixed.normal, 
                                 control.family=list(hyper=list(prec=utc.30as.prec.pc)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.up3), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.up3)

utc.output.pred.15as.up3$dic$dic
utc.output.pred.15as.up3$waic$waic
slcpo(utc.output.pred.15as.up3)

### extract parameter estimates
utc.fixed.out.up3 <- round(utc.output.pred.15as.up3$summary.fixed[,1:5],3) 
utc.fixed1_marg.up3  <- utc.output.pred.15as.up3$marginals.fixed[[1]]
utc.fixed2_marg.up3  <- utc.output.pred.15as.up3$marginals.fixed[[2]]
utc.fixed3_marg.up3  <- utc.output.pred.15as.up3$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.up3 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.up3$marginals.hyperpar[[1]])
utc.sigma2e_m1.up3 <- inla.emarginal(function(x) x, utc.sigma2e_marg.up3)
utc.sigma2e_m2.up3 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.up3)
utc.sigma2e_stdev.up3 <- sqrt(utc.sigma2e_m2.up3 - utc.sigma2e_m1.up3^2)
utc.sigma2e_quantiles.up3 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.up3)

utc.spde.result.up3 <- inla.spde2.result(utc.output.pred.15as.up3, name="utc.30as.spatial.field.pc", utc.30as.spde.pc)

# variance of spatial effects
utc.var.nom.marg.up3 <- utc.spde.result.up3$marginals.variance.nominal[[1]]
utc.var.nom.m1.up3 <- inla.emarginal(function(x) x, utc.var.nom.marg.up3)
utc.var.nom.m2.up3 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.up3)
utc.var.nom.stdev.up3 <- sqrt(utc.var.nom.m2.up3 - utc.var.nom.m1.up3^2)
utc.var.nom.quantiles.up3 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.up3)

utc.range.nom.marg.up3 <- utc.spde.result.up3$marginals.range.nominal[[1]]
utc.range.nom.m1.up3 <- inla.emarginal(function(x) x, utc.range.nom.marg.up3)
utc.range.nom.m2.up3 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.up3)
utc.range.nom.stdev.up3 <- sqrt(utc.range.nom.m2.up3 - utc.range.nom.m1.up3^2)
utc.range.nom.quantiles.up3 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.up3)

### extract predicted values
utc.index.pred.response.up3 <- inla.stack.index(utc.join.stack.pred.15as.up3, "pred.response")$data 
utc.index.pred.latent.up3 <- inla.stack.index(utc.join.stack.pred.15as.up3, "pred.latent")$data 

utc.pred.response.mean.up3 <- utc.output.pred.15as.up3$summary.linear.predictor[utc.index.pred.response.up3, "mean"]
utc.pred.response.sd.up3 <- utc.output.pred.15as.up3$summary.linear.predictor[utc.index.pred.response.up3, "sd"]

utc.pred.latent.mean.up3 <- utc.output.pred.15as.up3$summary.linear.predictor[utc.index.pred.latent.up3, "mean"]
utc.pred.latent.sd.up3 <- utc.output.pred.15as.up3$summary.linear.predictor[utc.index.pred.latent.up3, "sd"]

utc_pred_resp_mean.up3 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.up3.mod/4+1))
utc_pred_resp_mean_df.up3 <- as.data.frame(mapply(c,utc_pred_resp_mean.up3, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.up3 <- orderBy(~y + x, data=utc_pred_resp_mean_df.up3)
utc.pred.response.mean.grid.up3 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.up3$`log(utc.pred.response.mean.up3.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.up3 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.up3)
utc_pred_resp_sd_df.up3 <- as.data.frame(mapply(c,utc_pred_resp_sd.up3, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.up3 <- orderBy(~y + x, data=utc_pred_resp_sd_df.up3)
utc.pred.response.sd.grid.up3 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.up3$utc.pred.response.sd.up3)), nrow=utc.x.res*2)

utc_pred_lat_mean.up3 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.up3)
utc_pred_lat_mean_df.up3 <- as.data.frame(mapply(c,utc_pred_lat_mean.up3, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.up3 <- orderBy(~y + x, data=utc_pred_lat_mean_df.up3)
utc.pred.latent.mean.grid.up3 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.up3$utc.pred.latent.mean.up3)), nrow=utc.x.res*2)

utc_pred_lat_sd.up3 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.up3)
utc_pred_lat_sd_df.up3 <- as.data.frame(mapply(c,utc_pred_lat_sd.up3, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.up3 <- orderBy(~y + x, data=utc_pred_lat_sd_df.up3)
utc.pred.latent.sd.grid.up3 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.up3$utc.pred.latent.sd.up3)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.up3),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.up3), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.up3-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.up3-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.up3-utc.pred.latent.mean.grid.up3)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.up3-utc.pred.latent.mean.grid.up3)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.up3,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.up3, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.up3,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.up3, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.up3.mod <- exp(utc.pred.response.mean.up3)-1
utc.pred.response.mean.up3.mod[utc.pred.response.mean.up3.mod<0]<-0
res.15as.up3<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.up3.mod/4
RMSE.15as.up3 <- sqrt(mean(res.15as.up3^2))
RMSE.15as.up3 
correl.15as.up3 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.up3.mod/4) 
correl.15as.up3 

####################################
# real case sampled scenario + worst prior
# full model w/ log transformation + pc prior est at 30 arcsec full data
# est at 3 arcsec sampled 10% data (not weighted) + pred at 3 arcsec
####################################
utc.stack.est.15as.up4 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_u4$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_u4, 1),
                                     effects = list(c(utc.30as.s.index.pc, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_u4$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_u4$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.up4 <- inla.stack(utc.stack.est.15as.up4, utc.stack.latent.wp.p, utc.stack.pred.15as.wp.p)

utc.output.pred.15as.up4 <- inla(utc.formula.log.30as.pc, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.up4, spde=utc.30as.spde.pc), family="gaussian",
                                 control.fixed=utc.30as.fixed.normal, 
                                 control.family=list(hyper=list(prec=utc.30as.prec.pc)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.up4), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.up4)

utc.output.pred.15as.up4$dic$dic
utc.output.pred.15as.up4$waic$waic
slcpo(utc.output.pred.15as.up4)

### extract parameter estimates
utc.fixed.out.up4 <- round(utc.output.pred.15as.up4$summary.fixed[,1:5],3) 
utc.fixed1_marg.up4  <- utc.output.pred.15as.up4$marginals.fixed[[1]]
utc.fixed2_marg.up4  <- utc.output.pred.15as.up4$marginals.fixed[[2]]
utc.fixed3_marg.up4  <- utc.output.pred.15as.up4$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.up4 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.up4$marginals.hyperpar[[1]])
utc.sigma2e_m1.up4 <- inla.emarginal(function(x) x, utc.sigma2e_marg.up4)
utc.sigma2e_m2.up4 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.up4)
utc.sigma2e_stdev.up4 <- sqrt(utc.sigma2e_m2.up4 - utc.sigma2e_m1.up4^2)
utc.sigma2e_quantiles.up4 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.up4)

utc.spde.result.up4 <- inla.spde2.result(utc.output.pred.15as.up4, name="utc.30as.spatial.field.pc", utc.30as.spde.pc)

# variance of spatial effects
utc.var.nom.marg.up4 <- utc.spde.result.up4$marginals.variance.nominal[[1]]
utc.var.nom.m1.up4 <- inla.emarginal(function(x) x, utc.var.nom.marg.up4)
utc.var.nom.m2.up4 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.up4)
utc.var.nom.stdev.up4 <- sqrt(utc.var.nom.m2.up4 - utc.var.nom.m1.up4^2)
utc.var.nom.quantiles.up4 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.up4)

utc.range.nom.marg.up4 <- utc.spde.result.up4$marginals.range.nominal[[1]]
utc.range.nom.m1.up4 <- inla.emarginal(function(x) x, utc.range.nom.marg.up4)
utc.range.nom.m2.up4 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.up4)
utc.range.nom.stdev.up4 <- sqrt(utc.range.nom.m2.up4 - utc.range.nom.m1.up4^2)
utc.range.nom.quantiles.up4 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.up4)

### extract predicted values
utc.index.pred.response.up4 <- inla.stack.index(utc.join.stack.pred.15as.up4, "pred.response")$data 
utc.index.pred.latent.up4 <- inla.stack.index(utc.join.stack.pred.15as.up4, "pred.latent")$data 

utc.pred.response.mean.up4 <- utc.output.pred.15as.up4$summary.linear.predictor[utc.index.pred.response.up4, "mean"]
utc.pred.response.sd.up4 <- utc.output.pred.15as.up4$summary.linear.predictor[utc.index.pred.response.up4, "sd"]

utc.pred.latent.mean.up4 <- utc.output.pred.15as.up4$summary.linear.predictor[utc.index.pred.latent.up4, "mean"]
utc.pred.latent.sd.up4 <- utc.output.pred.15as.up4$summary.linear.predictor[utc.index.pred.latent.up4, "sd"]

utc_pred_resp_mean.up4 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.up4.mod/4+1))
utc_pred_resp_mean_df.up4 <- as.data.frame(mapply(c,utc_pred_resp_mean.up4, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.up4 <- orderBy(~y + x, data=utc_pred_resp_mean_df.up4)
utc.pred.response.mean.grid.up4 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.up4$`log(utc.pred.response.mean.up4.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.up4 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.up4)
utc_pred_resp_sd_df.up4 <- as.data.frame(mapply(c,utc_pred_resp_sd.up4, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.up4 <- orderBy(~y + x, data=utc_pred_resp_sd_df.up4)
utc.pred.response.sd.grid.up4 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.up4$utc.pred.response.sd.up4)), nrow=utc.x.res*2)

utc_pred_lat_mean.up4 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.up4)
utc_pred_lat_mean_df.up4 <- as.data.frame(mapply(c,utc_pred_lat_mean.up4, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.up4 <- orderBy(~y + x, data=utc_pred_lat_mean_df.up4)
utc.pred.latent.mean.grid.up4 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.up4$utc.pred.latent.mean.up4)), nrow=utc.x.res*2)

utc_pred_lat_sd.up4 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.up4)
utc_pred_lat_sd_df.up4 <- as.data.frame(mapply(c,utc_pred_lat_sd.up4, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.up4 <- orderBy(~y + x, data=utc_pred_lat_sd_df.up4)
utc.pred.latent.sd.grid.up4 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.up4$utc.pred.latent.sd.up4)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.up4),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.up4), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.up4-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.up4-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.up4-utc.pred.latent.mean.grid.up4)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.up4-utc.pred.latent.mean.grid.up4)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.up4,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.up4, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.up4,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.up4, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.up4.mod <- exp(utc.pred.response.mean.up4)-1
utc.pred.response.mean.up4.mod[utc.pred.response.mean.up4.mod<0]<-0
res.15as.up4<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.up4.mod/4
RMSE.15as.up4 <- sqrt(mean(res.15as.up4^2))
RMSE.15as.up4 
correl.15as.up4 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.up4.mod/4) 
correl.15as.up4 

####################################
# real case sampled scenario + worst prior
# full model w/ log transformation + pc prior est at 30 arcsec full data
# est at 3 arcsec sampled 20% data (not weighted) + pred at 3 arcsec
####################################
utc.stack.est.15as.up5 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_u5$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_u5, 1),
                                     effects = list(c(utc.30as.s.index.pc, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_u5$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_u5$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.up5 <- inla.stack(utc.stack.est.15as.up5, utc.stack.latent.wp.p, utc.stack.pred.15as.wp.p)

utc.output.pred.15as.up5 <- inla(utc.formula.log.30as.pc, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.up5, spde=utc.30as.spde.pc), family="gaussian",
                                 control.fixed=utc.30as.fixed.normal, 
                                 control.family=list(hyper=list(prec=utc.30as.prec.pc)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.up5), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.up5)

utc.output.pred.15as.up5$dic$dic
utc.output.pred.15as.up5$waic$waic
slcpo(utc.output.pred.15as.up5)

### extract parameter estimates
utc.fixed.out.up5 <- round(utc.output.pred.15as.up5$summary.fixed[,1:5],3) 
utc.fixed1_marg.up5  <- utc.output.pred.15as.up5$marginals.fixed[[1]]
utc.fixed2_marg.up5  <- utc.output.pred.15as.up5$marginals.fixed[[2]]
utc.fixed3_marg.up5  <- utc.output.pred.15as.up5$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.up5 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.up5$marginals.hyperpar[[1]])
utc.sigma2e_m1.up5 <- inla.emarginal(function(x) x, utc.sigma2e_marg.up5)
utc.sigma2e_m2.up5 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.up5)
utc.sigma2e_stdev.up5 <- sqrt(utc.sigma2e_m2.up5 - utc.sigma2e_m1.up5^2)
utc.sigma2e_quantiles.up5 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.up5)

utc.spde.result.up5 <- inla.spde2.result(utc.output.pred.15as.up5, name="utc.30as.spatial.field.pc", utc.30as.spde.pc)

# variance of spatial effects
utc.var.nom.marg.up5 <- utc.spde.result.up5$marginals.variance.nominal[[1]]
utc.var.nom.m1.up5 <- inla.emarginal(function(x) x, utc.var.nom.marg.up5)
utc.var.nom.m2.up5 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.up5)
utc.var.nom.stdev.up5 <- sqrt(utc.var.nom.m2.up5 - utc.var.nom.m1.up5^2)
utc.var.nom.quantiles.up5 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.up5)

utc.range.nom.marg.up5 <- utc.spde.result.up5$marginals.range.nominal[[1]]
utc.range.nom.m1.up5 <- inla.emarginal(function(x) x, utc.range.nom.marg.up5)
utc.range.nom.m2.up5 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.up5)
utc.range.nom.stdev.up5 <- sqrt(utc.range.nom.m2.up5 - utc.range.nom.m1.up5^2)
utc.range.nom.quantiles.up5 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.up5)

### extract predicted values
utc.index.pred.response.up5 <- inla.stack.index(utc.join.stack.pred.15as.up5, "pred.response")$data 
utc.index.pred.latent.up5 <- inla.stack.index(utc.join.stack.pred.15as.up5, "pred.latent")$data 

utc.pred.response.mean.up5 <- utc.output.pred.15as.up5$summary.linear.predictor[utc.index.pred.response.up5, "mean"]
utc.pred.response.sd.up5 <- utc.output.pred.15as.up5$summary.linear.predictor[utc.index.pred.response.up5, "sd"]

utc.pred.latent.mean.up5 <- utc.output.pred.15as.up5$summary.linear.predictor[utc.index.pred.latent.up5, "mean"]
utc.pred.latent.sd.up5 <- utc.output.pred.15as.up5$summary.linear.predictor[utc.index.pred.latent.up5, "sd"]

utc_pred_resp_mean.up5 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.up5.mod/4+1))
utc_pred_resp_mean_df.up5 <- as.data.frame(mapply(c,utc_pred_resp_mean.up5, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.up5 <- orderBy(~y + x, data=utc_pred_resp_mean_df.up5)
utc.pred.response.mean.grid.up5 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.up5$`log(utc.pred.response.mean.up5.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.up5 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.up5)
utc_pred_resp_sd_df.up5 <- as.data.frame(mapply(c,utc_pred_resp_sd.up5, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.up5 <- orderBy(~y + x, data=utc_pred_resp_sd_df.up5)
utc.pred.response.sd.grid.up5 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.up5$utc.pred.response.sd.up5)), nrow=utc.x.res*2)

utc_pred_lat_mean.up5 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.up5)
utc_pred_lat_mean_df.up5 <- as.data.frame(mapply(c,utc_pred_lat_mean.up5, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.up5 <- orderBy(~y + x, data=utc_pred_lat_mean_df.up5)
utc.pred.latent.mean.grid.up5 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.up5$utc.pred.latent.mean.up5)), nrow=utc.x.res*2)

utc_pred_lat_sd.up5 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.up5)
utc_pred_lat_sd_df.up5 <- as.data.frame(mapply(c,utc_pred_lat_sd.up5, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.up5 <- orderBy(~y + x, data=utc_pred_lat_sd_df.up5)
utc.pred.latent.sd.grid.up5 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.up5$utc.pred.latent.sd.up5)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.up5),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.up5), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.up5-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.up5-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.up5-utc.pred.latent.mean.grid.up5)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.up5-utc.pred.latent.mean.grid.up5)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.up5,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.up5, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.up5,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.up5, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.up5.mod <- exp(utc.pred.response.mean.up5)-1
utc.pred.response.mean.up5.mod[utc.pred.response.mean.up5.mod<0]<-0
res.15as.up5<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.up5.mod/4
RMSE.15as.up5 <- sqrt(mean(res.15as.up5^2))
RMSE.15as.up5 
correl.15as.up5 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.up5.mod/4) 
correl.15as.up5 

####################################
# real case sampled scenario + worst prior
# full model w/ log transformation + pc prior est at 30 arcsec full data
# est at 3 arcsec sampled 50% data (not weighted) + pred at 3 arcsec
####################################
utc.stack.est.15as.up6 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_u6$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_u6, 1),
                                     effects = list(c(utc.30as.s.index.pc, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_u6$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_u6$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.up6 <- inla.stack(utc.stack.est.15as.up6, utc.stack.latent.wp.p, utc.stack.pred.15as.wp.p)

utc.output.pred.15as.up6 <- inla(utc.formula.log.30as.pc, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.up6, spde=utc.30as.spde.pc), family="gaussian",
                                 control.fixed=utc.30as.fixed.normal, 
                                 control.family=list(hyper=list(prec=utc.30as.prec.pc)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.up6), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.up6)

utc.output.pred.15as.up6$dic$dic
utc.output.pred.15as.up6$waic$waic
slcpo(utc.output.pred.15as.up6)

### extract parameter estimates
utc.fixed.out.up6 <- round(utc.output.pred.15as.up6$summary.fixed[,1:5],3)
utc.fixed1_marg.up6  <- utc.output.pred.15as.up6$marginals.fixed[[1]]
utc.fixed2_marg.up6  <- utc.output.pred.15as.up6$marginals.fixed[[2]]
utc.fixed3_marg.up6  <- utc.output.pred.15as.up6$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.up6 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.up6$marginals.hyperpar[[1]])
utc.sigma2e_m1.up6 <- inla.emarginal(function(x) x, utc.sigma2e_marg.up6)
utc.sigma2e_m2.up6 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.up6)
utc.sigma2e_stdev.up6 <- sqrt(utc.sigma2e_m2.up6 - utc.sigma2e_m1.up6^2)
utc.sigma2e_quantiles.up6 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.up6)

utc.spde.result.up6 <- inla.spde2.result(utc.output.pred.15as.up6, name="utc.30as.spatial.field.pc", utc.30as.spde.pc)

# variance of spatial effects
utc.var.nom.marg.up6 <- utc.spde.result.up6$marginals.variance.nominal[[1]]
utc.var.nom.m1.up6 <- inla.emarginal(function(x) x, utc.var.nom.marg.up6)
utc.var.nom.m2.up6 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.up6)
utc.var.nom.stdev.up6 <- sqrt(utc.var.nom.m2.up6 - utc.var.nom.m1.up6^2)
utc.var.nom.quantiles.up6 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.up6)

utc.range.nom.marg.up6 <- utc.spde.result.up6$marginals.range.nominal[[1]]
utc.range.nom.m1.up6 <- inla.emarginal(function(x) x, utc.range.nom.marg.up6)
utc.range.nom.m2.up6 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.up6)
utc.range.nom.stdev.up6 <- sqrt(utc.range.nom.m2.up6 - utc.range.nom.m1.up6^2)
utc.range.nom.quantiles.up6 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.up6)

### extract predicted values
utc.index.pred.response.up6 <- inla.stack.index(utc.join.stack.pred.15as.up6, "pred.response")$data 
utc.index.pred.latent.up6 <- inla.stack.index(utc.join.stack.pred.15as.up6, "pred.latent")$data 

utc.pred.response.mean.up6 <- utc.output.pred.15as.up6$summary.linear.predictor[utc.index.pred.response.up6, "mean"]
utc.pred.response.sd.up6 <- utc.output.pred.15as.up6$summary.linear.predictor[utc.index.pred.response.up6, "sd"]

utc.pred.latent.mean.up6 <- utc.output.pred.15as.up6$summary.linear.predictor[utc.index.pred.latent.up6, "mean"]
utc.pred.latent.sd.up6 <- utc.output.pred.15as.up6$summary.linear.predictor[utc.index.pred.latent.up6, "sd"]

utc_pred_resp_mean.up6 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.up6.mod/4+1))
utc_pred_resp_mean_df.up6 <- as.data.frame(mapply(c,utc_pred_resp_mean.up6, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.up6 <- orderBy(~y + x, data=utc_pred_resp_mean_df.up6)
utc.pred.response.mean.grid.up6 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.up6$`log(utc.pred.response.mean.up6.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.up6 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.up6)
utc_pred_resp_sd_df.up6 <- as.data.frame(mapply(c,utc_pred_resp_sd.up6, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.up6 <- orderBy(~y + x, data=utc_pred_resp_sd_df.up6)
utc.pred.response.sd.grid.up6 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.up6$utc.pred.response.sd.up6)), nrow=utc.x.res*2)

utc_pred_lat_mean.up6 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.up6)
utc_pred_lat_mean_df.up6 <- as.data.frame(mapply(c,utc_pred_lat_mean.up6, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.up6 <- orderBy(~y + x, data=utc_pred_lat_mean_df.up6)
utc.pred.latent.mean.grid.up6 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.up6$utc.pred.latent.mean.up6)), nrow=utc.x.res*2)

utc_pred_lat_sd.up6 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.up6)
utc_pred_lat_sd_df.up6 <- as.data.frame(mapply(c,utc_pred_lat_sd.up6, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.up6 <- orderBy(~y + x, data=utc_pred_lat_sd_df.up6)
utc.pred.latent.sd.grid.up6 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.up6$utc.pred.latent.sd.up6)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.up6),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.up6), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.up6-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.up6-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.up6-utc.pred.latent.mean.grid.up6)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.up6-utc.pred.latent.mean.grid.up6)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.up6,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.up6, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.up6,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.up6, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.up6.mod <- exp(utc.pred.response.mean.up6)-1
utc.pred.response.mean.up6.mod[utc.pred.response.mean.up6.mod<0]<-0
res.15as.up6<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.up6.mod/4
RMSE.15as.up6 <- sqrt(mean(res.15as.up6^2))
RMSE.15as.up6 
correl.15as.up6 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.up6.mod/4) 
correl.15as.up6 

####################################
# real case sampled scenario + vague prior
# full model w/ log transformation + vague prior
# est at 3 arcsec sampled 1% data (not weighted) + pred at 3 arcsec
####################################
utc.stack.est.15as.uv1 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_u1$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_u1, 1),
                                     effects = list(c(utc.s.index, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_u1$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_u1$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.uv1 <- inla.stack(utc.stack.est.15as.uv1, utc.stack.latent, utc.stack.pred.15as.log)

utc.output.pred.15as.uv1 <- inla(utc.formula.log, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.uv1, spde=utc.spde), family="gaussian",
                                 control.fixed=utc.fixed.vague, 
                                 control.family=list(hyper=list(prec=utc.prec.vague)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.uv1), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.uv1)

utc.output.pred.15as.uv1$dic$dic
utc.output.pred.15as.uv1$waic$waic
slcpo(utc.output.pred.15as.uv1)

### extract parameter estimates
utc.fixed.out.uv1 <- round(utc.output.pred.15as.uv1$summary.fixed[,1:5],3) 
utc.fixed1_marg.uv1  <- utc.output.pred.15as.uv1$marginals.fixed[[1]]
utc.fixed2_marg.uv1  <- utc.output.pred.15as.uv1$marginals.fixed[[2]]
utc.fixed3_marg.uv1  <- utc.output.pred.15as.uv1$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.uv1 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.uv1$marginals.hyperpar[[1]])
utc.sigma2e_m1.uv1 <- inla.emarginal(function(x) x, utc.sigma2e_marg.uv1)
utc.sigma2e_m2.uv1 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.uv1)
utc.sigma2e_stdev.uv1 <- sqrt(utc.sigma2e_m2.uv1 - utc.sigma2e_m1.uv1^2)
utc.sigma2e_quantiles.uv1 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.uv1)

utc.spde.result.uv1 <- inla.spde2.result(utc.output.pred.15as.uv1, name="utc.spatial.field", utc.spde)

# variance of spatial effects
utc.var.nom.marg.uv1 <- utc.spde.result.uv1$marginals.variance.nominal[[1]]
utc.var.nom.m1.uv1 <- inla.emarginal(function(x) x, utc.var.nom.marg.uv1)
utc.var.nom.m2.uv1 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.uv1)
utc.var.nom.stdev.uv1 <- sqrt(utc.var.nom.m2.uv1 - utc.var.nom.m1.uv1^2)
utc.var.nom.quantiles.uv1 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.uv1)

utc.range.nom.marg.uv1 <- utc.spde.result.uv1$marginals.range.nominal[[1]]
utc.range.nom.m1.uv1 <- inla.emarginal(function(x) x, utc.range.nom.marg.uv1)
utc.range.nom.m2.uv1 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.uv1)
utc.range.nom.stdev.uv1 <- sqrt(utc.range.nom.m2.uv1 - utc.range.nom.m1.uv1^2)
utc.range.nom.quantiles.uv1 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.uv1)

### extract predicted values
utc.index.pred.response.uv1 <- inla.stack.index(utc.join.stack.pred.15as.uv1, "pred.response")$data 
utc.index.pred.latent.uv1 <- inla.stack.index(utc.join.stack.pred.15as.uv1, "pred.latent")$data 

utc.pred.response.mean.uv1 <- utc.output.pred.15as.uv1$summary.linear.predictor[utc.index.pred.response.uv1, "mean"]
utc.pred.response.sd.uv1 <- utc.output.pred.15as.uv1$summary.linear.predictor[utc.index.pred.response.uv1, "sd"]

utc.pred.latent.mean.uv1 <- utc.output.pred.15as.uv1$summary.linear.predictor[utc.index.pred.latent.uv1, "mean"]
utc.pred.latent.sd.uv1 <- utc.output.pred.15as.uv1$summary.linear.predictor[utc.index.pred.latent.uv1, "sd"]

utc_pred_resp_mean.uv1 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.uv1.mod/4+1))
utc_pred_resp_mean_df.uv1 <- as.data.frame(mapply(c,utc_pred_resp_mean.uv1, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.uv1 <- orderBy(~y + x, data=utc_pred_resp_mean_df.uv1)
utc.pred.response.mean.grid.uv1 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.uv1$`log(utc.pred.response.mean.uv1.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.uv1 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.uv1)
utc_pred_resp_sd_df.uv1 <- as.data.frame(mapply(c,utc_pred_resp_sd.uv1, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.uv1 <- orderBy(~y + x, data=utc_pred_resp_sd_df.uv1)
utc.pred.response.sd.grid.uv1 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.uv1$utc.pred.response.sd.uv1)), nrow=utc.x.res*2)

utc_pred_lat_mean.uv1 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.uv1)
utc_pred_lat_mean_df.uv1 <- as.data.frame(mapply(c,utc_pred_lat_mean.uv1, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.uv1 <- orderBy(~y + x, data=utc_pred_lat_mean_df.uv1)
utc.pred.latent.mean.grid.uv1 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.uv1$utc.pred.latent.mean.uv1)), nrow=utc.x.res*2)

utc_pred_lat_sd.uv1 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.uv1)
utc_pred_lat_sd_df.uv1 <- as.data.frame(mapply(c,utc_pred_lat_sd.uv1, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.uv1 <- orderBy(~y + x, data=utc_pred_lat_sd_df.uv1)
utc.pred.latent.sd.grid.uv1 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.uv1$utc.pred.latent.sd.uv1)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.uv1),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.uv1), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.uv1-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.uv1-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.uv1-utc.pred.latent.mean.grid.uv1)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.uv1-utc.pred.latent.mean.grid.uv1)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.uv1,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.uv1, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.uv1,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.uv1, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.uv1.mod <- exp(utc.pred.response.mean.uv1)-1
utc.pred.response.mean.uv1.mod[utc.pred.response.mean.uv1.mod<0]<-0
res.15as.uv1<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.uv1.mod/4
RMSE.15as.uv1 <- sqrt(mean(res.15as.uv1^2))
RMSE.15as.uv1 
correl.15as.uv1 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.uv1.mod/4) 
correl.15as.uv1 

####################################
# real case sampled scenario + vague prior
# full model w/ log transformation + vague prior
# est at 3 arcsec sampled 2% data (not weighted) + pred at 3 arcsec
####################################
utc.stack.est.15as.uv2 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_u2$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_u2, 1),
                                     effects = list(c(utc.s.index, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_u2$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_u2$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.uv2 <- inla.stack(utc.stack.est.15as.uv2, utc.stack.latent, utc.stack.pred.15as.log)

utc.output.pred.15as.uv2 <- inla(utc.formula.log, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.uv2, spde=utc.spde), family="gaussian",
                                 control.fixed=utc.fixed.vague, 
                                 control.family=list(hyper=list(prec=utc.prec.vague)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.uv2), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.uv2)

utc.output.pred.15as.uv2$dic$dic
utc.output.pred.15as.uv2$waic$waic
slcpo(utc.output.pred.15as.uv2)

### extract parameter estimates
utc.fixed.out.uv2 <- round(utc.output.pred.15as.uv2$summary.fixed[,1:5],3) 
utc.fixed1_marg.uv2  <- utc.output.pred.15as.uv2$marginals.fixed[[1]]
utc.fixed2_marg.uv2  <- utc.output.pred.15as.uv2$marginals.fixed[[2]]
utc.fixed3_marg.uv2  <- utc.output.pred.15as.uv2$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.uv2 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.uv2$marginals.hyperpar[[1]])
utc.sigma2e_m1.uv2 <- inla.emarginal(function(x) x, utc.sigma2e_marg.uv2)
utc.sigma2e_m2.uv2 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.uv2)
utc.sigma2e_stdev.uv2 <- sqrt(utc.sigma2e_m2.uv2 - utc.sigma2e_m1.uv2^2)
utc.sigma2e_quantiles.uv2 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.uv2)

utc.spde.result.uv2 <- inla.spde2.result(utc.output.pred.15as.uv2, name="utc.spatial.field", utc.spde)

# variance of spatial effects
utc.var.nom.marg.uv2 <- utc.spde.result.uv2$marginals.variance.nominal[[1]]
utc.var.nom.m1.uv2 <- inla.emarginal(function(x) x, utc.var.nom.marg.uv2)
utc.var.nom.m2.uv2 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.uv2)
utc.var.nom.stdev.uv2 <- sqrt(utc.var.nom.m2.uv2 - utc.var.nom.m1.uv2^2)
utc.var.nom.quantiles.uv2 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.uv2)

utc.range.nom.marg.uv2 <- utc.spde.result.uv2$marginals.range.nominal[[1]]
utc.range.nom.m1.uv2 <- inla.emarginal(function(x) x, utc.range.nom.marg.uv2)
utc.range.nom.m2.uv2 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.uv2)
utc.range.nom.stdev.uv2 <- sqrt(utc.range.nom.m2.uv2 - utc.range.nom.m1.uv2^2)
utc.range.nom.quantiles.uv2 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.uv2)

### extract predicted values
utc.index.pred.response.uv2 <- inla.stack.index(utc.join.stack.pred.15as.uv2, "pred.response")$data 
utc.index.pred.latent.uv2 <- inla.stack.index(utc.join.stack.pred.15as.uv2, "pred.latent")$data 

utc.pred.response.mean.uv2 <- utc.output.pred.15as.uv2$summary.linear.predictor[utc.index.pred.response.uv2, "mean"]
utc.pred.response.sd.uv2 <- utc.output.pred.15as.uv2$summary.linear.predictor[utc.index.pred.response.uv2, "sd"]

utc.pred.latent.mean.uv2 <- utc.output.pred.15as.uv2$summary.linear.predictor[utc.index.pred.latent.uv2, "mean"]
utc.pred.latent.sd.uv2 <- utc.output.pred.15as.uv2$summary.linear.predictor[utc.index.pred.latent.uv2, "sd"]

utc_pred_resp_mean.uv2 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.uv2.mod/4+1))
utc_pred_resp_mean_df.uv2 <- as.data.frame(mapply(c,utc_pred_resp_mean.uv2, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.uv2 <- orderBy(~y + x, data=utc_pred_resp_mean_df.uv2)
utc.pred.response.mean.grid.uv2 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.uv2$`log(utc.pred.response.mean.uv2.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.uv2 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.uv2)
utc_pred_resp_sd_df.uv2 <- as.data.frame(mapply(c,utc_pred_resp_sd.uv2, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.uv2 <- orderBy(~y + x, data=utc_pred_resp_sd_df.uv2)
utc.pred.response.sd.grid.uv2 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.uv2$utc.pred.response.sd.uv2)), nrow=utc.x.res*2)

utc_pred_lat_mean.uv2 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.uv2)
utc_pred_lat_mean_df.uv2 <- as.data.frame(mapply(c,utc_pred_lat_mean.uv2, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.uv2 <- orderBy(~y + x, data=utc_pred_lat_mean_df.uv2)
utc.pred.latent.mean.grid.uv2 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.uv2$utc.pred.latent.mean.uv2)), nrow=utc.x.res*2)

utc_pred_lat_sd.uv2 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.uv2)
utc_pred_lat_sd_df.uv2 <- as.data.frame(mapply(c,utc_pred_lat_sd.uv2, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.uv2 <- orderBy(~y + x, data=utc_pred_lat_sd_df.uv2)
utc.pred.latent.sd.grid.uv2 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.uv2$utc.pred.latent.sd.uv2)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.uv2),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.uv2), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.uv2-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.uv2-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.uv2-utc.pred.latent.mean.grid.uv2)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.uv2-utc.pred.latent.mean.grid.uv2)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.uv2,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.uv2, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.uv2,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.uv2, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.uv2.mod <- exp(utc.pred.response.mean.uv2)-1
utc.pred.response.mean.uv2.mod[utc.pred.response.mean.uv2.mod<0]<-0
res.15as.uv2<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.uv2.mod/4
RMSE.15as.uv2 <- sqrt(mean(res.15as.uv2^2))
RMSE.15as.uv2 
correl.15as.uv2 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.uv2.mod/4) 
correl.15as.uv2 

####################################
# real case sampled scenario + vague prior
# full model w/ log transformation + vague prior
# est at 3 arcsec sampled 5% data (not weighted) + pred at 3 arcsec
####################################
utc.stack.est.15as.uv3 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_u3$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_u3, 1),
                                     effects = list(c(utc.s.index, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_u3$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_u3$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.uv3 <- inla.stack(utc.stack.est.15as.uv3, utc.stack.latent, utc.stack.pred.15as.log)

utc.output.pred.15as.uv3 <- inla(utc.formula.log, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.uv3, spde=utc.spde), family="gaussian",
                                 control.fixed=utc.fixed.vague, 
                                 control.family=list(hyper=list(prec=utc.prec.vague)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.uv3), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.uv3)

utc.output.pred.15as.uv3$dic$dic
utc.output.pred.15as.uv3$waic$waic
slcpo(utc.output.pred.15as.uv3)

### extract parameter estimates
utc.fixed.out.uv3 <- round(utc.output.pred.15as.uv3$summary.fixed[,1:5],3) 
utc.fixed1_marg.uv3  <- utc.output.pred.15as.uv3$marginals.fixed[[1]]
utc.fixed2_marg.uv3  <- utc.output.pred.15as.uv3$marginals.fixed[[2]]
utc.fixed3_marg.uv3  <- utc.output.pred.15as.uv3$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.uv3 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.uv3$marginals.hyperpar[[1]])
utc.sigma2e_m1.uv3 <- inla.emarginal(function(x) x, utc.sigma2e_marg.uv3)
utc.sigma2e_m2.uv3 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.uv3)
utc.sigma2e_stdev.uv3 <- sqrt(utc.sigma2e_m2.uv3 - utc.sigma2e_m1.uv3^2)
utc.sigma2e_quantiles.uv3 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.uv3)

utc.spde.result.uv3 <- inla.spde2.result(utc.output.pred.15as.uv3, name="utc.spatial.field", utc.spde)

# variance of spatial effects
utc.var.nom.marg.uv3 <- utc.spde.result.uv3$marginals.variance.nominal[[1]]
utc.var.nom.m1.uv3 <- inla.emarginal(function(x) x, utc.var.nom.marg.uv3)
utc.var.nom.m2.uv3 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.uv3)
utc.var.nom.stdev.uv3 <- sqrt(utc.var.nom.m2.uv3 - utc.var.nom.m1.uv3^2)
utc.var.nom.quantiles.uv3 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.uv3)

utc.range.nom.marg.uv3 <- utc.spde.result.uv3$marginals.range.nominal[[1]]
utc.range.nom.m1.uv3 <- inla.emarginal(function(x) x, utc.range.nom.marg.uv3)
utc.range.nom.m2.uv3 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.uv3)
utc.range.nom.stdev.uv3 <- sqrt(utc.range.nom.m2.uv3 - utc.range.nom.m1.uv3^2)
utc.range.nom.quantiles.uv3 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.uv3)

### extract predicted values
utc.index.pred.response.uv3 <- inla.stack.index(utc.join.stack.pred.15as.uv3, "pred.response")$data 
utc.index.pred.latent.uv3 <- inla.stack.index(utc.join.stack.pred.15as.uv3, "pred.latent")$data 

utc.pred.response.mean.uv3 <- utc.output.pred.15as.uv3$summary.linear.predictor[utc.index.pred.response.uv3, "mean"]
utc.pred.response.sd.uv3 <- utc.output.pred.15as.uv3$summary.linear.predictor[utc.index.pred.response.uv3, "sd"]

utc.pred.latent.mean.uv3 <- utc.output.pred.15as.uv3$summary.linear.predictor[utc.index.pred.latent.uv3, "mean"]
utc.pred.latent.sd.uv3 <- utc.output.pred.15as.uv3$summary.linear.predictor[utc.index.pred.latent.uv3, "sd"]

utc_pred_resp_mean.uv3 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.uv3.mod/4+1))
utc_pred_resp_mean_df.uv3 <- as.data.frame(mapply(c,utc_pred_resp_mean.uv3, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.uv3 <- orderBy(~y + x, data=utc_pred_resp_mean_df.uv3)
utc.pred.response.mean.grid.uv3 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.uv3$`log(utc.pred.response.mean.uv3.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.uv3 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.uv3)
utc_pred_resp_sd_df.uv3 <- as.data.frame(mapply(c,utc_pred_resp_sd.uv3, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.uv3 <- orderBy(~y + x, data=utc_pred_resp_sd_df.uv3)
utc.pred.response.sd.grid.uv3 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.uv3$utc.pred.response.sd.uv3)), nrow=utc.x.res*2)

utc_pred_lat_mean.uv3 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.uv3)
utc_pred_lat_mean_df.uv3 <- as.data.frame(mapply(c,utc_pred_lat_mean.uv3, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.uv3 <- orderBy(~y + x, data=utc_pred_lat_mean_df.uv3)
utc.pred.latent.mean.grid.uv3 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.uv3$utc.pred.latent.mean.uv3)), nrow=utc.x.res*2)

utc_pred_lat_sd.uv3 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.uv3)
utc_pred_lat_sd_df.uv3 <- as.data.frame(mapply(c,utc_pred_lat_sd.uv3, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.uv3 <- orderBy(~y + x, data=utc_pred_lat_sd_df.uv3)
utc.pred.latent.sd.grid.uv3 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.uv3$utc.pred.latent.sd.uv3)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.uv3),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.uv3), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.uv3-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.uv3-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.uv3-utc.pred.latent.mean.grid.uv3)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.uv3-utc.pred.latent.mean.grid.uv3)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.uv3,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.uv3, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.uv3,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.uv3, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.uv3.mod <- exp(utc.pred.response.mean.uv3)-1
utc.pred.response.mean.uv3.mod[utc.pred.response.mean.uv3.mod<0]<-0
res.15as.uv3<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.uv3.mod/4
RMSE.15as.uv3 <- sqrt(mean(res.15as.uv3^2))
RMSE.15as.uv3 
correl.15as.uv3 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.uv3.mod/4) 
correl.15as.uv3 

####################################
# real case sampled scenario + vague prior
# full model w/ log transformation + vague prior
# est at 3 arcsec sampled 10% data (not weighted) + pred at 3 arcsec
####################################
utc.stack.est.15as.uv4 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_u4$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_u4, 1),
                                     effects = list(c(utc.s.index, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_u4$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_u4$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.uv4 <- inla.stack(utc.stack.est.15as.uv4, utc.stack.latent, utc.stack.pred.15as.log)

utc.output.pred.15as.uv4 <- inla(utc.formula.log, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.uv4, spde=utc.spde), family="gaussian",
                                 control.fixed=utc.fixed.vague, 
                                 control.family=list(hyper=list(prec=utc.prec.vague)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.uv4), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.uv4)

utc.output.pred.15as.uv4$dic$dic
utc.output.pred.15as.uv4$waic$waic
slcpo(utc.output.pred.15as.uv4)

### extract parameter estimates
utc.fixed.out.uv4 <- round(utc.output.pred.15as.uv4$summary.fixed[,1:5],3)
utc.fixed1_marg.uv4  <- utc.output.pred.15as.uv4$marginals.fixed[[1]]
utc.fixed2_marg.uv4  <- utc.output.pred.15as.uv4$marginals.fixed[[2]]
utc.fixed3_marg.uv4  <- utc.output.pred.15as.uv4$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.uv4 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.uv4$marginals.hyperpar[[1]])
utc.sigma2e_m1.uv4 <- inla.emarginal(function(x) x, utc.sigma2e_marg.uv4)
utc.sigma2e_m2.uv4 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.uv4)
utc.sigma2e_stdev.uv4 <- sqrt(utc.sigma2e_m2.uv4 - utc.sigma2e_m1.uv4^2)
utc.sigma2e_quantiles.uv4 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.uv4)

utc.spde.result.uv4 <- inla.spde2.result(utc.output.pred.15as.uv4, name="utc.spatial.field", utc.spde)

# variance of spatial effects
utc.var.nom.marg.uv4 <- utc.spde.result.uv4$marginals.variance.nominal[[1]]
utc.var.nom.m1.uv4 <- inla.emarginal(function(x) x, utc.var.nom.marg.uv4)
utc.var.nom.m2.uv4 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.uv4)
utc.var.nom.stdev.uv4 <- sqrt(utc.var.nom.m2.uv4 - utc.var.nom.m1.uv4^2)
utc.var.nom.quantiles.uv4 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.uv4)

utc.range.nom.marg.uv4 <- utc.spde.result.uv4$marginals.range.nominal[[1]]
utc.range.nom.m1.uv4 <- inla.emarginal(function(x) x, utc.range.nom.marg.uv4)
utc.range.nom.m2.uv4 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.uv4)
utc.range.nom.stdev.uv4 <- sqrt(utc.range.nom.m2.uv4 - utc.range.nom.m1.uv4^2)
utc.range.nom.quantiles.uv4 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.uv4)

### extract predicted values
utc.index.pred.response.uv4 <- inla.stack.index(utc.join.stack.pred.15as.uv4, "pred.response")$data 
utc.index.pred.latent.uv4 <- inla.stack.index(utc.join.stack.pred.15as.uv4, "pred.latent")$data 

utc.pred.response.mean.uv4 <- utc.output.pred.15as.uv4$summary.linear.predictor[utc.index.pred.response.uv4, "mean"]
utc.pred.response.sd.uv4 <- utc.output.pred.15as.uv4$summary.linear.predictor[utc.index.pred.response.uv4, "sd"]

utc.pred.latent.mean.uv4 <- utc.output.pred.15as.uv4$summary.linear.predictor[utc.index.pred.latent.uv4, "mean"]
utc.pred.latent.sd.uv4 <- utc.output.pred.15as.uv4$summary.linear.predictor[utc.index.pred.latent.uv4, "sd"]

utc_pred_resp_mean.uv4 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.uv4.mod/4+1))
utc_pred_resp_mean_df.uv4 <- as.data.frame(mapply(c,utc_pred_resp_mean.uv4, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.uv4 <- orderBy(~y + x, data=utc_pred_resp_mean_df.uv4)
utc.pred.response.mean.grid.uv4 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.uv4$`log(utc.pred.response.mean.uv4.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.uv4 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.uv4)
utc_pred_resp_sd_df.uv4 <- as.data.frame(mapply(c,utc_pred_resp_sd.uv4, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.uv4 <- orderBy(~y + x, data=utc_pred_resp_sd_df.uv4)
utc.pred.response.sd.grid.uv4 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.uv4$utc.pred.response.sd.uv4)), nrow=utc.x.res*2)

utc_pred_lat_mean.uv4 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.uv4)
utc_pred_lat_mean_df.uv4 <- as.data.frame(mapply(c,utc_pred_lat_mean.uv4, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.uv4 <- orderBy(~y + x, data=utc_pred_lat_mean_df.uv4)
utc.pred.latent.mean.grid.uv4 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.uv4$utc.pred.latent.mean.uv4)), nrow=utc.x.res*2)

utc_pred_lat_sd.uv4 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.uv4)
utc_pred_lat_sd_df.uv4 <- as.data.frame(mapply(c,utc_pred_lat_sd.uv4, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.uv4 <- orderBy(~y + x, data=utc_pred_lat_sd_df.uv4)
utc.pred.latent.sd.grid.uv4 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.uv4$utc.pred.latent.sd.uv4)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.uv4),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.uv4), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.uv4-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.uv4-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.uv4-utc.pred.latent.mean.grid.uv4)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.uv4-utc.pred.latent.mean.grid.uv4)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.uv4,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.uv4, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.uv4,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.uv4, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.uv4.mod <- exp(utc.pred.response.mean.uv4)-1
utc.pred.response.mean.uv4.mod[utc.pred.response.mean.uv4.mod<0]<-0
res.15as.uv4<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.uv4.mod/4
RMSE.15as.uv4 <- sqrt(mean(res.15as.uv4^2))
RMSE.15as.uv4 
correl.15as.uv4 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.uv4.mod/4) 
correl.15as.uv4 

####################################
# real case sampled scenario + vague prior
# full model w/ log transformation + vague prior
# est at 3 arcsec sampled 20% data (not weighted) + pred at 3 arcsec
####################################
utc.stack.est.15as.uv5 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_u5$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_u5, 1),
                                     effects = list(c(utc.s.index, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_u5$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_u5$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.uv5 <- inla.stack(utc.stack.est.15as.uv5, utc.stack.latent, utc.stack.pred.15as.log)

utc.output.pred.15as.uv5 <- inla(utc.formula.log, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.uv5, spde=utc.spde), family="gaussian",
                                 control.fixed=utc.fixed.vague, 
                                 control.family=list(hyper=list(prec=utc.prec.vague)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.uv5), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.uv5)

utc.output.pred.15as.uv5$dic$dic
utc.output.pred.15as.uv5$waic$waic
slcpo(utc.output.pred.15as.uv5)

### extract parameter estimates
utc.fixed.out.uv5 <- round(utc.output.pred.15as.uv5$summary.fixed[,1:5],3) 
utc.fixed1_marg.uv5  <- utc.output.pred.15as.uv5$marginals.fixed[[1]]
utc.fixed2_marg.uv5  <- utc.output.pred.15as.uv5$marginals.fixed[[2]]
utc.fixed3_marg.uv5  <- utc.output.pred.15as.uv5$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.uv5 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.uv5$marginals.hyperpar[[1]])
utc.sigma2e_m1.uv5 <- inla.emarginal(function(x) x, utc.sigma2e_marg.uv5)
utc.sigma2e_m2.uv5 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.uv5)
utc.sigma2e_stdev.uv5 <- sqrt(utc.sigma2e_m2.uv5 - utc.sigma2e_m1.uv5^2)
utc.sigma2e_quantiles.uv5 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.uv5)

utc.spde.result.uv5 <- inla.spde2.result(utc.output.pred.15as.uv5, name="utc.spatial.field", utc.spde)

# variance of spatial effects
utc.var.nom.marg.uv5 <- utc.spde.result.uv5$marginals.variance.nominal[[1]]
utc.var.nom.m1.uv5 <- inla.emarginal(function(x) x, utc.var.nom.marg.uv5)
utc.var.nom.m2.uv5 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.uv5)
utc.var.nom.stdev.uv5 <- sqrt(utc.var.nom.m2.uv5 - utc.var.nom.m1.uv5^2)
utc.var.nom.quantiles.uv5 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.uv5)

utc.range.nom.marg.uv5 <- utc.spde.result.uv5$marginals.range.nominal[[1]]
utc.range.nom.m1.uv5 <- inla.emarginal(function(x) x, utc.range.nom.marg.uv5)
utc.range.nom.m2.uv5 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.uv5)
utc.range.nom.stdev.uv5 <- sqrt(utc.range.nom.m2.uv5 - utc.range.nom.m1.uv5^2)
utc.range.nom.quantiles.uv5 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.uv5)

### extract predicted values
utc.index.pred.response.uv5 <- inla.stack.index(utc.join.stack.pred.15as.uv5, "pred.response")$data 
utc.index.pred.latent.uv5 <- inla.stack.index(utc.join.stack.pred.15as.uv5, "pred.latent")$data 

utc.pred.response.mean.uv5 <- utc.output.pred.15as.uv5$summary.linear.predictor[utc.index.pred.response.uv5, "mean"]
utc.pred.response.sd.uv5 <- utc.output.pred.15as.uv5$summary.linear.predictor[utc.index.pred.response.uv5, "sd"]

utc.pred.latent.mean.uv5 <- utc.output.pred.15as.uv5$summary.linear.predictor[utc.index.pred.latent.uv5, "mean"]
utc.pred.latent.sd.uv5 <- utc.output.pred.15as.uv5$summary.linear.predictor[utc.index.pred.latent.uv5, "sd"]

utc_pred_resp_mean.uv5 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.uv5.mod/4+1))
utc_pred_resp_mean_df.uv5 <- as.data.frame(mapply(c,utc_pred_resp_mean.uv5, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.uv5 <- orderBy(~y + x, data=utc_pred_resp_mean_df.uv5)
utc.pred.response.mean.grid.uv5 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.uv5$`log(utc.pred.response.mean.uv5.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.uv5 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.uv5)
utc_pred_resp_sd_df.uv5 <- as.data.frame(mapply(c,utc_pred_resp_sd.uv5, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.uv5 <- orderBy(~y + x, data=utc_pred_resp_sd_df.uv5)
utc.pred.response.sd.grid.uv5 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.uv5$utc.pred.response.sd.uv5)), nrow=utc.x.res*2)

utc_pred_lat_mean.uv5 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.uv5)
utc_pred_lat_mean_df.uv5 <- as.data.frame(mapply(c,utc_pred_lat_mean.uv5, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.uv5 <- orderBy(~y + x, data=utc_pred_lat_mean_df.uv5)
utc.pred.latent.mean.grid.uv5 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.uv5$utc.pred.latent.mean.uv5)), nrow=utc.x.res*2)

utc_pred_lat_sd.uv5 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.uv5)
utc_pred_lat_sd_df.uv5 <- as.data.frame(mapply(c,utc_pred_lat_sd.uv5, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.uv5 <- orderBy(~y + x, data=utc_pred_lat_sd_df.uv5)
utc.pred.latent.sd.grid.uv5 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.uv5$utc.pred.latent.sd.uv5)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.uv5),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.uv5), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.uv5-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.uv5-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.uv5-utc.pred.latent.mean.grid.uv5)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.uv5-utc.pred.latent.mean.grid.uv5)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.uv5,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.uv5, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.uv5,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.uv5, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.uv5.mod <- exp(utc.pred.response.mean.uv5)-1
utc.pred.response.mean.uv5.mod[utc.pred.response.mean.uv5.mod<0]<-0
res.15as.uv5<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.uv5.mod/4
RMSE.15as.uv5 <- sqrt(mean(res.15as.uv5^2))
RMSE.15as.uv5 
correl.15as.uv5 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.uv5.mod/4) 
correl.15as.uv5 

####################################
# real case sampled scenario + vague prior
# full model w/ log transformation + vague prior
# est at 3 arcsec sampled 50% data (not weighted) + pred at 3 arcsec
####################################
utc.stack.est.15as.uv6 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_u6$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_u6, 1),
                                     effects = list(c(utc.s.index, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_u6$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_u6$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.uv6 <- inla.stack(utc.stack.est.15as.uv6, utc.stack.latent, utc.stack.pred.15as.log)

utc.output.pred.15as.uv6 <- inla(utc.formula.log, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.uv6, spde=utc.spde), family="gaussian",
                                 control.fixed=utc.fixed.vague, 
                                 control.family=list(hyper=list(prec=utc.prec.vague)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.uv6), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.uv6)

utc.output.pred.15as.uv6$dic$dic
utc.output.pred.15as.uv6$waic$waic
slcpo(utc.output.pred.15as.uv6)

### extract parameter estimates
utc.fixed.out.uv6 <- round(utc.output.pred.15as.uv6$summary.fixed[,1:5],3) 
utc.fixed1_marg.uv6  <- utc.output.pred.15as.uv6$marginals.fixed[[1]]
utc.fixed2_marg.uv6  <- utc.output.pred.15as.uv6$marginals.fixed[[2]]
utc.fixed3_marg.uv6  <- utc.output.pred.15as.uv6$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.uv6 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.uv6$marginals.hyperpar[[1]])
utc.sigma2e_m1.uv6 <- inla.emarginal(function(x) x, utc.sigma2e_marg.uv6)
utc.sigma2e_m2.uv6 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.uv6)
utc.sigma2e_stdev.uv6 <- sqrt(utc.sigma2e_m2.uv6 - utc.sigma2e_m1.uv6^2)
utc.sigma2e_quantiles.uv6 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.uv6)

utc.spde.result.uv6 <- inla.spde2.result(utc.output.pred.15as.uv6, name="utc.spatial.field", utc.spde)

# variance of spatial effects
utc.var.nom.marg.uv6 <- utc.spde.result.uv6$marginals.variance.nominal[[1]]
utc.var.nom.m1.uv6 <- inla.emarginal(function(x) x, utc.var.nom.marg.uv6)
utc.var.nom.m2.uv6 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.uv6)
utc.var.nom.stdev.uv6 <- sqrt(utc.var.nom.m2.uv6 - utc.var.nom.m1.uv6^2)
utc.var.nom.quantiles.uv6 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.uv6)

utc.range.nom.marg.uv6 <- utc.spde.result.uv6$marginals.range.nominal[[1]]
utc.range.nom.m1.uv6 <- inla.emarginal(function(x) x, utc.range.nom.marg.uv6)
utc.range.nom.m2.uv6 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.uv6)
utc.range.nom.stdev.uv6 <- sqrt(utc.range.nom.m2.uv6 - utc.range.nom.m1.uv6^2)
utc.range.nom.quantiles.uv6 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.uv6)

### extract predicted values
utc.index.pred.response.uv6 <- inla.stack.index(utc.join.stack.pred.15as.uv6, "pred.response")$data 
utc.index.pred.latent.uv6 <- inla.stack.index(utc.join.stack.pred.15as.uv6, "pred.latent")$data 

utc.pred.response.mean.uv6 <- utc.output.pred.15as.uv6$summary.linear.predictor[utc.index.pred.response.uv6, "mean"]
utc.pred.response.sd.uv6 <- utc.output.pred.15as.uv6$summary.linear.predictor[utc.index.pred.response.uv6, "sd"]

utc.pred.latent.mean.uv6 <- utc.output.pred.15as.uv6$summary.linear.predictor[utc.index.pred.latent.uv6, "mean"]
utc.pred.latent.sd.uv6 <- utc.output.pred.15as.uv6$summary.linear.predictor[utc.index.pred.latent.uv6, "sd"]

utc_pred_resp_mean.uv6 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.uv6.mod/4+1))
utc_pred_resp_mean_df.uv6 <- as.data.frame(mapply(c,utc_pred_resp_mean.uv6, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.uv6 <- orderBy(~y + x, data=utc_pred_resp_mean_df.uv6)
utc.pred.response.mean.grid.uv6 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.uv6$`log(utc.pred.response.mean.uv6.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.uv6 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.uv6)
utc_pred_resp_sd_df.uv6 <- as.data.frame(mapply(c,utc_pred_resp_sd.uv6, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.uv6 <- orderBy(~y + x, data=utc_pred_resp_sd_df.uv6)
utc.pred.response.sd.grid.uv6 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.uv6$utc.pred.response.sd.uv6)), nrow=utc.x.res*2)

utc_pred_lat_mean.uv6 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.uv6)
utc_pred_lat_mean_df.uv6 <- as.data.frame(mapply(c,utc_pred_lat_mean.uv6, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.uv6 <- orderBy(~y + x, data=utc_pred_lat_mean_df.uv6)
utc.pred.latent.mean.grid.uv6 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.uv6$utc.pred.latent.mean.uv6)), nrow=utc.x.res*2)

utc_pred_lat_sd.uv6 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.uv6)
utc_pred_lat_sd_df.uv6 <- as.data.frame(mapply(c,utc_pred_lat_sd.uv6, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.uv6 <- orderBy(~y + x, data=utc_pred_lat_sd_df.uv6)
utc.pred.latent.sd.grid.uv6 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.uv6$utc.pred.latent.sd.uv6)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.uv6),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.uv6), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.uv6-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.uv6-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.uv6-utc.pred.latent.mean.grid.uv6)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.uv6-utc.pred.latent.mean.grid.uv6)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.uv6,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.uv6, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.uv6,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.uv6, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.uv6.mod <- exp(utc.pred.response.mean.uv6)-1
utc.pred.response.mean.uv6.mod[utc.pred.response.mean.uv6.mod<0]<-0
res.15as.uv6<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.uv6.mod/4
RMSE.15as.uv6 <- sqrt(mean(res.15as.uv6^2))
RMSE.15as.uv6 
correl.15as.uv6 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.uv6.mod/4) 
correl.15as.uv6 

####################################
# real case sampled scenario + worst prior
# full model w/ log transformation + normal prior est at 30 arcsec full data
# est at 3 arcsec sampled 1% data (weighted) + pred at 3 arcsec
####################################
utc.loc_sampled_15as_w1 <- cbind(stacked_utc_sampled_w1$x, stacked_utc_sampled_w1$y)
utc.A.sampled.15as_w1 <- inla.spde.make.A(mesh=utc.mesh, loc=utc.loc_sampled_15as_w1)

utc.stack.est.15as.wn1 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_w1$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_w1, 1),
                                     effects = list(c(utc.30as.s.index.normal, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_w1$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_w1$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.wn1 <- inla.stack(utc.stack.est.15as.wn1, utc.stack.latent.wp.n, utc.stack.pred.15as.wp.n)

utc.output.pred.15as.wn1 <- inla(utc.formula.log.30as.normal, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.wn1, spde=utc.30as.spde.normal), family="gaussian",
                                 control.fixed=utc.30as.fixed.normal, 
                                 control.family=list(hyper=list(prec=utc.30as.prec.normal)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.wn1), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.wn1)

utc.output.pred.15as.wn1$dic$dic
utc.output.pred.15as.wn1$waic$waic
slcpo(utc.output.pred.15as.wn1)

### extract parameter estimates
utc.fixed.out.wn1 <- round(utc.output.pred.15as.wn1$summary.fixed[,1:5],3) 
utc.fixed1_marg.wn1  <- utc.output.pred.15as.wn1$marginals.fixed[[1]]
utc.fixed2_marg.wn1  <- utc.output.pred.15as.wn1$marginals.fixed[[2]]
utc.fixed3_marg.wn1  <- utc.output.pred.15as.wn1$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.wn1 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.wn1$marginals.hyperpar[[1]])
utc.sigma2e_m1.wn1 <- inla.emarginal(function(x) x, utc.sigma2e_marg.wn1)
utc.sigma2e_m2.wn1 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.wn1)
utc.sigma2e_stdev.wn1 <- sqrt(utc.sigma2e_m2.wn1 - utc.sigma2e_m1.wn1^2)
utc.sigma2e_quantiles.wn1 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.wn1)

utc.spde.result.wn1 <- inla.spde2.result(utc.output.pred.15as.wn1, name="utc.30as.spatial.field.normal", utc.30as.spde.normal)

# variance of spatial effects
utc.var.nom.marg.wn1 <- utc.spde.result.wn1$marginals.variance.nominal[[1]]
utc.var.nom.m1.wn1 <- inla.emarginal(function(x) x, utc.var.nom.marg.wn1)
utc.var.nom.m2.wn1 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.wn1)
utc.var.nom.stdev.wn1 <- sqrt(utc.var.nom.m2.wn1 - utc.var.nom.m1.wn1^2)
utc.var.nom.quantiles.wn1 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.wn1)

utc.range.nom.marg.wn1 <- utc.spde.result.wn1$marginals.range.nominal[[1]]
utc.range.nom.m1.wn1 <- inla.emarginal(function(x) x, utc.range.nom.marg.wn1)
utc.range.nom.m2.wn1 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.wn1)
utc.range.nom.stdev.wn1 <- sqrt(utc.range.nom.m2.wn1 - utc.range.nom.m1.wn1^2)
utc.range.nom.quantiles.wn1 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.wn1)

### extract predicted values
utc.index.pred.response.wn1 <- inla.stack.index(utc.join.stack.pred.15as.wn1, "pred.response")$data 
utc.index.pred.latent.wn1 <- inla.stack.index(utc.join.stack.pred.15as.wn1, "pred.latent")$data 

utc.pred.response.mean.wn1 <- utc.output.pred.15as.wn1$summary.linear.predictor[utc.index.pred.response.wn1, "mean"]
utc.pred.response.sd.wn1 <- utc.output.pred.15as.wn1$summary.linear.predictor[utc.index.pred.response.wn1, "sd"]

utc.pred.latent.mean.wn1 <- utc.output.pred.15as.wn1$summary.linear.predictor[utc.index.pred.latent.wn1, "mean"]
utc.pred.latent.sd.wn1 <- utc.output.pred.15as.wn1$summary.linear.predictor[utc.index.pred.latent.wn1, "sd"]

utc_pred_resp_mean.wn1 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.wn1.mod/4+1))
utc_pred_resp_mean_df.wn1 <- as.data.frame(mapply(c,utc_pred_resp_mean.wn1, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.wn1 <- orderBy(~y + x, data=utc_pred_resp_mean_df.wn1)
utc.pred.response.mean.grid.wn1 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.wn1$`log(utc.pred.response.mean.wn1.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.wn1 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.wn1)
utc_pred_resp_sd_df.wn1 <- as.data.frame(mapply(c,utc_pred_resp_sd.wn1, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.wn1 <- orderBy(~y + x, data=utc_pred_resp_sd_df.wn1)
utc.pred.response.sd.grid.wn1 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.wn1$utc.pred.response.sd.wn1)), nrow=utc.x.res*2)

utc_pred_lat_mean.wn1 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.wn1)
utc_pred_lat_mean_df.wn1 <- as.data.frame(mapply(c,utc_pred_lat_mean.wn1, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.wn1 <- orderBy(~y + x, data=utc_pred_lat_mean_df.wn1)
utc.pred.latent.mean.grid.wn1 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.wn1$utc.pred.latent.mean.wn1)), nrow=utc.x.res*2)

utc_pred_lat_sd.wn1 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.wn1)
utc_pred_lat_sd_df.wn1 <- as.data.frame(mapply(c,utc_pred_lat_sd.wn1, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.wn1 <- orderBy(~y + x, data=utc_pred_lat_sd_df.wn1)
utc.pred.latent.sd.grid.wn1 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.wn1$utc.pred.latent.sd.wn1)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wn1),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wn1), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wn1-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wn1-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wn1-utc.pred.latent.mean.grid.wn1)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wn1-utc.pred.latent.mean.grid.wn1)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wn1,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wn1, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wn1,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wn1, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.wn1.mod <- exp(utc.pred.response.mean.wn1)-1
utc.pred.response.mean.wn1.mod[utc.pred.response.mean.wn1.mod<0]<-0
res.15as.wn1<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.wn1.mod/4
RMSE.15as.wn1 <- sqrt(mean(res.15as.wn1^2))
RMSE.15as.wn1 
correl.15as.wn1 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.wn1.mod/4) 
correl.15as.wn1 


####################################
# real case sampled scenario + worst prior
# full model w/ log transformation + normal prior est at 30 arcsec full data
# est at 3 arcsec sampled 2% data (weighted) + pred at 3 arcsec
####################################
utc.loc_sampled_15as_w2 <- cbind(stacked_utc_sampled_w2$x, stacked_utc_sampled_w2$y)
utc.A.sampled.15as_w2 <- inla.spde.make.A(mesh=utc.mesh, loc=utc.loc_sampled_15as_w2)

utc.stack.est.15as.wn2 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_w2$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_w2, 1),
                                     effects = list(c(utc.30as.s.index.normal, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_w2$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_w2$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.wn2 <- inla.stack(utc.stack.est.15as.wn2, utc.stack.latent.wp.n, utc.stack.pred.15as.wp.n)

utc.output.pred.15as.wn2 <- inla(utc.formula.log.30as.normal, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.wn2, spde=utc.30as.spde.normal), family="gaussian",
                                 control.fixed=utc.30as.fixed.normal, 
                                 control.family=list(hyper=list(prec=utc.30as.prec.normal)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.wn2), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.wn2)

utc.output.pred.15as.wn2$dic$dic
utc.output.pred.15as.wn2$waic$waic
slcpo(utc.output.pred.15as.wn2)

### extract parameter estimates
utc.fixed.out.wn2 <- round(utc.output.pred.15as.wn2$summary.fixed[,1:5],3) 
utc.fixed1_marg.wn2  <- utc.output.pred.15as.wn2$marginals.fixed[[1]]
utc.fixed2_marg.wn2  <- utc.output.pred.15as.wn2$marginals.fixed[[2]]
utc.fixed3_marg.wn2  <- utc.output.pred.15as.wn2$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.wn2 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.wn2$marginals.hyperpar[[1]])
utc.sigma2e_m1.wn2 <- inla.emarginal(function(x) x, utc.sigma2e_marg.wn2)
utc.sigma2e_m2.wn2 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.wn2)
utc.sigma2e_stdev.wn2 <- sqrt(utc.sigma2e_m2.wn2 - utc.sigma2e_m1.wn2^2)
utc.sigma2e_quantiles.wn2 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.wn2)

utc.spde.result.wn2 <- inla.spde2.result(utc.output.pred.15as.wn2, name="utc.30as.spatial.field.normal", utc.30as.spde.normal)

# variance of spatial effects
utc.var.nom.marg.wn2 <- utc.spde.result.wn2$marginals.variance.nominal[[1]]
utc.var.nom.m1.wn2 <- inla.emarginal(function(x) x, utc.var.nom.marg.wn2)
utc.var.nom.m2.wn2 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.wn2)
utc.var.nom.stdev.wn2 <- sqrt(utc.var.nom.m2.wn2 - utc.var.nom.m1.wn2^2)
utc.var.nom.quantiles.wn2 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.wn2)

utc.range.nom.marg.wn2 <- utc.spde.result.wn2$marginals.range.nominal[[1]]
utc.range.nom.m1.wn2 <- inla.emarginal(function(x) x, utc.range.nom.marg.wn2)
utc.range.nom.m2.wn2 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.wn2)
utc.range.nom.stdev.wn2 <- sqrt(utc.range.nom.m2.wn2 - utc.range.nom.m1.wn2^2)
utc.range.nom.quantiles.wn2 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.wn2)

### extract predicted values
utc.index.pred.response.wn2 <- inla.stack.index(utc.join.stack.pred.15as.wn2, "pred.response")$data 
utc.index.pred.latent.wn2 <- inla.stack.index(utc.join.stack.pred.15as.wn2, "pred.latent")$data 

utc.pred.response.mean.wn2 <- utc.output.pred.15as.wn2$summary.linear.predictor[utc.index.pred.response.wn2, "mean"]
utc.pred.response.sd.wn2 <- utc.output.pred.15as.wn2$summary.linear.predictor[utc.index.pred.response.wn2, "sd"]

utc.pred.latent.mean.wn2 <- utc.output.pred.15as.wn2$summary.linear.predictor[utc.index.pred.latent.wn2, "mean"]
utc.pred.latent.sd.wn2 <- utc.output.pred.15as.wn2$summary.linear.predictor[utc.index.pred.latent.wn2, "sd"]

utc_pred_resp_mean.wn2 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.wn2.mod/4+1))
utc_pred_resp_mean_df.wn2 <- as.data.frame(mapply(c,utc_pred_resp_mean.wn2, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.wn2 <- orderBy(~y + x, data=utc_pred_resp_mean_df.wn2)
utc.pred.response.mean.grid.wn2 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.wn2$`log(utc.pred.response.mean.wn2.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.wn2 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.wn2)
utc_pred_resp_sd_df.wn2 <- as.data.frame(mapply(c,utc_pred_resp_sd.wn2, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.wn2 <- orderBy(~y + x, data=utc_pred_resp_sd_df.wn2)
utc.pred.response.sd.grid.wn2 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.wn2$utc.pred.response.sd.wn2)), nrow=utc.x.res*2)

utc_pred_lat_mean.wn2 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.wn2)
utc_pred_lat_mean_df.wn2 <- as.data.frame(mapply(c,utc_pred_lat_mean.wn2, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.wn2 <- orderBy(~y + x, data=utc_pred_lat_mean_df.wn2)
utc.pred.latent.mean.grid.wn2 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.wn2$utc.pred.latent.mean.wn2)), nrow=utc.x.res*2)

utc_pred_lat_sd.wn2 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.wn2)
utc_pred_lat_sd_df.wn2 <- as.data.frame(mapply(c,utc_pred_lat_sd.wn2, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.wn2 <- orderBy(~y + x, data=utc_pred_lat_sd_df.wn2)
utc.pred.latent.sd.grid.wn2 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.wn2$utc.pred.latent.sd.wn2)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wn2),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wn2), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wn2-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wn2-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wn2-utc.pred.latent.mean.grid.wn2)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wn2-utc.pred.latent.mean.grid.wn2)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wn2,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wn2, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wn2,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wn2, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.wn2.mod <- exp(utc.pred.response.mean.wn2)-1
utc.pred.response.mean.wn2.mod[utc.pred.response.mean.wn2.mod<0]<-0
res.15as.wn2<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.wn2.mod/4
RMSE.15as.wn2 <- sqrt(mean(res.15as.wn2^2))
RMSE.15as.wn2 
correl.15as.wn2 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.wn2.mod/4) 
correl.15as.wn2 

####################################
# real case sampled scenario + worst prior
# full model w/ log transformation + normal prior est at 30 arcsec full data
# est at 3 arcsec sampled 5% data (weighted) + pred at 3 arcsec
####################################
utc.loc_sampled_15as_w3 <- cbind(stacked_utc_sampled_w3$x, stacked_utc_sampled_w3$y)
utc.A.sampled.15as_w3 <- inla.spde.make.A(mesh=utc.mesh, loc=utc.loc_sampled_15as_w3)

utc.stack.est.15as.wn3 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_w3$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_w3, 1),
                                     effects = list(c(utc.30as.s.index.normal, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_w3$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_w3$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.wn3 <- inla.stack(utc.stack.est.15as.wn3, utc.stack.latent.wp.n, utc.stack.pred.15as.wp.n)

utc.output.pred.15as.wn3 <- inla(utc.formula.log.30as.normal, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.wn3, spde=utc.30as.spde.normal), family="gaussian",
                                 control.fixed=utc.30as.fixed.normal, 
                                 control.family=list(hyper=list(prec=utc.30as.prec.normal)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.wn3), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.wn3)

utc.output.pred.15as.wn3$dic$dic
utc.output.pred.15as.wn3$waic$waic
slcpo(utc.output.pred.15as.wn3)

### extract parameter estimates
utc.fixed.out.wn3 <- round(utc.output.pred.15as.wn3$summary.fixed[,1:5],3)
utc.fixed1_marg.wn3  <- utc.output.pred.15as.wn3$marginals.fixed[[1]]
utc.fixed2_marg.wn3  <- utc.output.pred.15as.wn3$marginals.fixed[[2]]
utc.fixed3_marg.wn3  <- utc.output.pred.15as.wn3$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.wn3 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.wn3$marginals.hyperpar[[1]])
utc.sigma2e_m1.wn3 <- inla.emarginal(function(x) x, utc.sigma2e_marg.wn3)
utc.sigma2e_m2.wn3 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.wn3)
utc.sigma2e_stdev.wn3 <- sqrt(utc.sigma2e_m2.wn3 - utc.sigma2e_m1.wn3^2)
utc.sigma2e_quantiles.wn3 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.wn3)

utc.spde.result.wn3 <- inla.spde2.result(utc.output.pred.15as.wn3, name="utc.30as.spatial.field.normal", utc.30as.spde.normal)

# variance of spatial effects
utc.var.nom.marg.wn3 <- utc.spde.result.wn3$marginals.variance.nominal[[1]]
utc.var.nom.m1.wn3 <- inla.emarginal(function(x) x, utc.var.nom.marg.wn3)
utc.var.nom.m2.wn3 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.wn3)
utc.var.nom.stdev.wn3 <- sqrt(utc.var.nom.m2.wn3 - utc.var.nom.m1.wn3^2)
utc.var.nom.quantiles.wn3 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.wn3)

utc.range.nom.marg.wn3 <- utc.spde.result.wn3$marginals.range.nominal[[1]]
utc.range.nom.m1.wn3 <- inla.emarginal(function(x) x, utc.range.nom.marg.wn3)
utc.range.nom.m2.wn3 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.wn3)
utc.range.nom.stdev.wn3 <- sqrt(utc.range.nom.m2.wn3 - utc.range.nom.m1.wn3^2)
utc.range.nom.quantiles.wn3 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.wn3)

### extract predicted values
utc.index.pred.response.wn3 <- inla.stack.index(utc.join.stack.pred.15as.wn3, "pred.response")$data 
utc.index.pred.latent.wn3 <- inla.stack.index(utc.join.stack.pred.15as.wn3, "pred.latent")$data 

utc.pred.response.mean.wn3 <- utc.output.pred.15as.wn3$summary.linear.predictor[utc.index.pred.response.wn3, "mean"]
utc.pred.response.sd.wn3 <- utc.output.pred.15as.wn3$summary.linear.predictor[utc.index.pred.response.wn3, "sd"]

utc.pred.latent.mean.wn3 <- utc.output.pred.15as.wn3$summary.linear.predictor[utc.index.pred.latent.wn3, "mean"]
utc.pred.latent.sd.wn3 <- utc.output.pred.15as.wn3$summary.linear.predictor[utc.index.pred.latent.wn3, "sd"]

utc_pred_resp_mean.wn3 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.wn3.mod/4+1))
utc_pred_resp_mean_df.wn3 <- as.data.frame(mapply(c,utc_pred_resp_mean.wn3, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.wn3 <- orderBy(~y + x, data=utc_pred_resp_mean_df.wn3)
utc.pred.response.mean.grid.wn3 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.wn3$`log(utc.pred.response.mean.wn3.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.wn3 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.wn3)
utc_pred_resp_sd_df.wn3 <- as.data.frame(mapply(c,utc_pred_resp_sd.wn3, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.wn3 <- orderBy(~y + x, data=utc_pred_resp_sd_df.wn3)
utc.pred.response.sd.grid.wn3 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.wn3$utc.pred.response.sd.wn3)), nrow=utc.x.res*2)

utc_pred_lat_mean.wn3 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.wn3)
utc_pred_lat_mean_df.wn3 <- as.data.frame(mapply(c,utc_pred_lat_mean.wn3, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.wn3 <- orderBy(~y + x, data=utc_pred_lat_mean_df.wn3)
utc.pred.latent.mean.grid.wn3 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.wn3$utc.pred.latent.mean.wn3)), nrow=utc.x.res*2)

utc_pred_lat_sd.wn3 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.wn3)
utc_pred_lat_sd_df.wn3 <- as.data.frame(mapply(c,utc_pred_lat_sd.wn3, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.wn3 <- orderBy(~y + x, data=utc_pred_lat_sd_df.wn3)
utc.pred.latent.sd.grid.wn3 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.wn3$utc.pred.latent.sd.wn3)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wn3),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wn3), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wn3-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wn3-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wn3-utc.pred.latent.mean.grid.wn3)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wn3-utc.pred.latent.mean.grid.wn3)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wn3,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wn3, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wn3,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wn3, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.wn3.mod <- exp(utc.pred.response.mean.wn3)-1
utc.pred.response.mean.wn3.mod[utc.pred.response.mean.wn3.mod<0]<-0
res.15as.wn3<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.wn3.mod/4
RMSE.15as.wn3 <- sqrt(mean(res.15as.wn3^2))
RMSE.15as.wn3 
correl.15as.wn3 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.wn3.mod/4) 
correl.15as.wn3 


####################################
# real case sampled scenario + worst prior
# full model w/ log transformation + normal prior est at 30 arcsec full data
# est at 3 arcsec sampled 10% data (weighted) + pred at 3 arcsec
####################################
utc.loc_sampled_15as_w4 <- cbind(stacked_utc_sampled_w4$x, stacked_utc_sampled_w4$y)
utc.A.sampled.15as_w4 <- inla.spde.make.A(mesh=utc.mesh, loc=utc.loc_sampled_15as_w4)

utc.stack.est.15as.wn4 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_w4$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_w4, 1),
                                     effects = list(c(utc.30as.s.index.normal, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_w4$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_w4$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.wn4 <- inla.stack(utc.stack.est.15as.wn4, utc.stack.latent.wp.n, utc.stack.pred.15as.wp.n)

utc.output.pred.15as.wn4 <- inla(utc.formula.log.30as.normal, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.wn4, spde=utc.30as.spde.normal), family="gaussian",
                                 control.fixed=utc.30as.fixed.normal, 
                                 control.family=list(hyper=list(prec=utc.30as.prec.normal)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.wn4), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.wn4)

utc.output.pred.15as.wn4$dic$dic
utc.output.pred.15as.wn4$waic$waic
slcpo(utc.output.pred.15as.wn4)

### extract parameter estimates
utc.fixed.out.wn4 <- round(utc.output.pred.15as.wn4$summary.fixed[,1:5],3) 
utc.fixed1_marg.wn4  <- utc.output.pred.15as.wn4$marginals.fixed[[1]]
utc.fixed2_marg.wn4  <- utc.output.pred.15as.wn4$marginals.fixed[[2]]
utc.fixed3_marg.wn4  <- utc.output.pred.15as.wn4$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.wn4 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.wn4$marginals.hyperpar[[1]])
utc.sigma2e_m1.wn4 <- inla.emarginal(function(x) x, utc.sigma2e_marg.wn4)
utc.sigma2e_m2.wn4 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.wn4)
utc.sigma2e_stdev.wn4 <- sqrt(utc.sigma2e_m2.wn4 - utc.sigma2e_m1.wn4^2)
utc.sigma2e_quantiles.wn4 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.wn4)

utc.spde.result.wn4 <- inla.spde2.result(utc.output.pred.15as.wn4, name="utc.30as.spatial.field.normal", utc.30as.spde.normal)

# variance of spatial effects
utc.var.nom.marg.wn4 <- utc.spde.result.wn4$marginals.variance.nominal[[1]]
utc.var.nom.m1.wn4 <- inla.emarginal(function(x) x, utc.var.nom.marg.wn4)
utc.var.nom.m2.wn4 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.wn4)
utc.var.nom.stdev.wn4 <- sqrt(utc.var.nom.m2.wn4 - utc.var.nom.m1.wn4^2)
utc.var.nom.quantiles.wn4 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.wn4)

utc.range.nom.marg.wn4 <- utc.spde.result.wn4$marginals.range.nominal[[1]]
utc.range.nom.m1.wn4 <- inla.emarginal(function(x) x, utc.range.nom.marg.wn4)
utc.range.nom.m2.wn4 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.wn4)
utc.range.nom.stdev.wn4 <- sqrt(utc.range.nom.m2.wn4 - utc.range.nom.m1.wn4^2)
utc.range.nom.quantiles.wn4 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.wn4)

### extract predicted values
utc.index.pred.response.wn4 <- inla.stack.index(utc.join.stack.pred.15as.wn4, "pred.response")$data 
utc.index.pred.latent.wn4 <- inla.stack.index(utc.join.stack.pred.15as.wn4, "pred.latent")$data 

utc.pred.response.mean.wn4 <- utc.output.pred.15as.wn4$summary.linear.predictor[utc.index.pred.response.wn4, "mean"]
utc.pred.response.sd.wn4 <- utc.output.pred.15as.wn4$summary.linear.predictor[utc.index.pred.response.wn4, "sd"]

utc.pred.latent.mean.wn4 <- utc.output.pred.15as.wn4$summary.linear.predictor[utc.index.pred.latent.wn4, "mean"]
utc.pred.latent.sd.wn4 <- utc.output.pred.15as.wn4$summary.linear.predictor[utc.index.pred.latent.wn4, "sd"]

utc_pred_resp_mean.wn4 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.wn4.mod/4+1))
utc_pred_resp_mean_df.wn4 <- as.data.frame(mapply(c,utc_pred_resp_mean.wn4, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.wn4 <- orderBy(~y + x, data=utc_pred_resp_mean_df.wn4)
utc.pred.response.mean.grid.wn4 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.wn4$`log(utc.pred.response.mean.wn4.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.wn4 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.wn4)
utc_pred_resp_sd_df.wn4 <- as.data.frame(mapply(c,utc_pred_resp_sd.wn4, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.wn4 <- orderBy(~y + x, data=utc_pred_resp_sd_df.wn4)
utc.pred.response.sd.grid.wn4 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.wn4$utc.pred.response.sd.wn4)), nrow=utc.x.res*2)

utc_pred_lat_mean.wn4 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.wn4)
utc_pred_lat_mean_df.wn4 <- as.data.frame(mapply(c,utc_pred_lat_mean.wn4, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.wn4 <- orderBy(~y + x, data=utc_pred_lat_mean_df.wn4)
utc.pred.latent.mean.grid.wn4 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.wn4$utc.pred.latent.mean.wn4)), nrow=utc.x.res*2)

utc_pred_lat_sd.wn4 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.wn4)
utc_pred_lat_sd_df.wn4 <- as.data.frame(mapply(c,utc_pred_lat_sd.wn4, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.wn4 <- orderBy(~y + x, data=utc_pred_lat_sd_df.wn4)
utc.pred.latent.sd.grid.wn4 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.wn4$utc.pred.latent.sd.wn4)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wn4),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wn4), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wn4-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wn4-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wn4-utc.pred.latent.mean.grid.wn4)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wn4-utc.pred.latent.mean.grid.wn4)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wn4,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wn4, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wn4,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wn4, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.wn4.mod <- exp(utc.pred.response.mean.wn4)-1
utc.pred.response.mean.wn4.mod[utc.pred.response.mean.wn4.mod<0]<-0
res.15as.wn4<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.wn4.mod/4
RMSE.15as.wn4 <- sqrt(mean(res.15as.wn4^2))
RMSE.15as.wn4 
correl.15as.wn4 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.wn4.mod/4) 
correl.15as.wn4 


####################################
# real case sampled scenario + worst prior
# full model w/ log transformation + normal prior est at 30 arcsec full data
# est at 3 arcsec sampled 20% data (weighted) + pred at 3 arcsec
####################################
utc.loc_sampled_15as_w5 <- cbind(stacked_utc_sampled_w5$x, stacked_utc_sampled_w5$y)
utc.A.sampled.15as_w5 <- inla.spde.make.A(mesh=utc.mesh, loc=utc.loc_sampled_15as_w5)

utc.stack.est.15as.wn5 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_w5$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_w5, 1),
                                     effects = list(c(utc.30as.s.index.normal, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_w5$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_w5$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.wn5 <- inla.stack(utc.stack.est.15as.wn5, utc.stack.latent.wp.n, utc.stack.pred.15as.wp.n)

utc.output.pred.15as.wn5 <- inla(utc.formula.log.30as.normal, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.wn5, spde=utc.30as.spde.normal), family="gaussian",
                                 control.fixed=utc.30as.fixed.normal, 
                                 control.family=list(hyper=list(prec=utc.30as.prec.normal)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.wn5), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.wn5)

utc.output.pred.15as.wn5$dic$dic
utc.output.pred.15as.wn5$waic$waic
slcpo(utc.output.pred.15as.wn5)

### extract parameter estimates
utc.fixed.out.wn5 <- round(utc.output.pred.15as.wn5$summary.fixed[,1:5],3) 
utc.fixed1_marg.wn5  <- utc.output.pred.15as.wn5$marginals.fixed[[1]]
utc.fixed2_marg.wn5  <- utc.output.pred.15as.wn5$marginals.fixed[[2]]
utc.fixed3_marg.wn5  <- utc.output.pred.15as.wn5$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.wn5 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.wn5$marginals.hyperpar[[1]])
utc.sigma2e_m1.wn5 <- inla.emarginal(function(x) x, utc.sigma2e_marg.wn5)
utc.sigma2e_m2.wn5 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.wn5)
utc.sigma2e_stdev.wn5 <- sqrt(utc.sigma2e_m2.wn5 - utc.sigma2e_m1.wn5^2)
utc.sigma2e_quantiles.wn5 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.wn5)

utc.spde.result.wn5 <- inla.spde2.result(utc.output.pred.15as.wn5, name="utc.30as.spatial.field.normal", utc.30as.spde.normal)

# variance of spatial effects
utc.var.nom.marg.wn5 <- utc.spde.result.wn5$marginals.variance.nominal[[1]]
utc.var.nom.m1.wn5 <- inla.emarginal(function(x) x, utc.var.nom.marg.wn5)
utc.var.nom.m2.wn5 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.wn5)
utc.var.nom.stdev.wn5 <- sqrt(utc.var.nom.m2.wn5 - utc.var.nom.m1.wn5^2)
utc.var.nom.quantiles.wn5 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.wn5)

utc.range.nom.marg.wn5 <- utc.spde.result.wn5$marginals.range.nominal[[1]]
utc.range.nom.m1.wn5 <- inla.emarginal(function(x) x, utc.range.nom.marg.wn5)
utc.range.nom.m2.wn5 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.wn5)
utc.range.nom.stdev.wn5 <- sqrt(utc.range.nom.m2.wn5 - utc.range.nom.m1.wn5^2)
utc.range.nom.quantiles.wn5 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.wn5)

### extract predicted values
utc.index.pred.response.wn5 <- inla.stack.index(utc.join.stack.pred.15as.wn5, "pred.response")$data 
utc.index.pred.latent.wn5 <- inla.stack.index(utc.join.stack.pred.15as.wn5, "pred.latent")$data 

utc.pred.response.mean.wn5 <- utc.output.pred.15as.wn5$summary.linear.predictor[utc.index.pred.response.wn5, "mean"]
utc.pred.response.sd.wn5 <- utc.output.pred.15as.wn5$summary.linear.predictor[utc.index.pred.response.wn5, "sd"]

utc.pred.latent.mean.wn5 <- utc.output.pred.15as.wn5$summary.linear.predictor[utc.index.pred.latent.wn5, "mean"]
utc.pred.latent.sd.wn5 <- utc.output.pred.15as.wn5$summary.linear.predictor[utc.index.pred.latent.wn5, "sd"]

utc_pred_resp_mean.wn5 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.wn5.mod/4+1))
utc_pred_resp_mean_df.wn5 <- as.data.frame(mapply(c,utc_pred_resp_mean.wn5, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.wn5 <- orderBy(~y + x, data=utc_pred_resp_mean_df.wn5)
utc.pred.response.mean.grid.wn5 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.wn5$`log(utc.pred.response.mean.wn5.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.wn5 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.wn5)
utc_pred_resp_sd_df.wn5 <- as.data.frame(mapply(c,utc_pred_resp_sd.wn5, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.wn5 <- orderBy(~y + x, data=utc_pred_resp_sd_df.wn5)
utc.pred.response.sd.grid.wn5 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.wn5$utc.pred.response.sd.wn5)), nrow=utc.x.res*2)

utc_pred_lat_mean.wn5 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.wn5)
utc_pred_lat_mean_df.wn5 <- as.data.frame(mapply(c,utc_pred_lat_mean.wn5, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.wn5 <- orderBy(~y + x, data=utc_pred_lat_mean_df.wn5)
utc.pred.latent.mean.grid.wn5 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.wn5$utc.pred.latent.mean.wn5)), nrow=utc.x.res*2)

utc_pred_lat_sd.wn5 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.wn5)
utc_pred_lat_sd_df.wn5 <- as.data.frame(mapply(c,utc_pred_lat_sd.wn5, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.wn5 <- orderBy(~y + x, data=utc_pred_lat_sd_df.wn5)
utc.pred.latent.sd.grid.wn5 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.wn5$utc.pred.latent.sd.wn5)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wn5),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wn5), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wn5-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wn5-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wn5-utc.pred.latent.mean.grid.wn5)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wn5-utc.pred.latent.mean.grid.wn5)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wn5,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wn5, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wn5,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wn5, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.wn5.mod <- exp(utc.pred.response.mean.wn5)-1
utc.pred.response.mean.wn5.mod[utc.pred.response.mean.wn5.mod<0]<-0
res.15as.wn5<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.wn5.mod/4
RMSE.15as.wn5 <- sqrt(mean(res.15as.wn5^2))
RMSE.15as.wn5 
correl.15as.wn5 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.wn5.mod/4) 
correl.15as.wn5 

####################################
# real case sampled scenario + worst prior
# full model w/ log transformation + normal prior est at 30 arcsec full data
# est at 3 arcsec sampled 50% data (weighted) + pred at 3 arcsec
####################################
utc.loc_sampled_15as_w6 <- cbind(stacked_utc_sampled_w6$x, stacked_utc_sampled_w6$y)
utc.A.sampled.15as_w6 <- inla.spde.make.A(mesh=utc.mesh, loc=utc.loc_sampled_15as_w6)

utc.stack.est.15as.wn6 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_w6$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_w6, 1),
                                     effects = list(c(utc.30as.s.index.normal, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_w6$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_w6$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.wn6 <- inla.stack(utc.stack.est.15as.wn6, utc.stack.latent.wp.n, utc.stack.pred.15as.wp.n)

utc.output.pred.15as.wn6 <- inla(utc.formula.log.30as.normal, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.wn6, spde=utc.30as.spde.normal), family="gaussian",
                                 control.fixed=utc.30as.fixed.normal, 
                                 control.family=list(hyper=list(prec=utc.30as.prec.normal)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.wn6), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.wn6)

utc.output.pred.15as.wn6$dic$dic
utc.output.pred.15as.wn6$waic$waic
slcpo(utc.output.pred.15as.wn6)

### extract parameter estimates
utc.fixed.out.wn6 <- round(utc.output.pred.15as.wn6$summary.fixed[,1:5],3) 
utc.fixed1_marg.wn6  <- utc.output.pred.15as.wn6$marginals.fixed[[1]]
utc.fixed2_marg.wn6  <- utc.output.pred.15as.wn6$marginals.fixed[[2]]
utc.fixed3_marg.wn6  <- utc.output.pred.15as.wn6$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.wn6 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.wn6$marginals.hyperpar[[1]])
utc.sigma2e_m1.wn6 <- inla.emarginal(function(x) x, utc.sigma2e_marg.wn6)
utc.sigma2e_m2.wn6 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.wn6)
utc.sigma2e_stdev.wn6 <- sqrt(utc.sigma2e_m2.wn6 - utc.sigma2e_m1.wn6^2)
utc.sigma2e_quantiles.wn6 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.wn6)

utc.spde.result.wn6 <- inla.spde2.result(utc.output.pred.15as.wn6, name="utc.30as.spatial.field.normal", utc.30as.spde.normal)

# variance of spatial effects
utc.var.nom.marg.wn6 <- utc.spde.result.wn6$marginals.variance.nominal[[1]]
utc.var.nom.m1.wn6 <- inla.emarginal(function(x) x, utc.var.nom.marg.wn6)
utc.var.nom.m2.wn6 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.wn6)
utc.var.nom.stdev.wn6 <- sqrt(utc.var.nom.m2.wn6 - utc.var.nom.m1.wn6^2)
utc.var.nom.quantiles.wn6 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.wn6)

utc.range.nom.marg.wn6 <- utc.spde.result.wn6$marginals.range.nominal[[1]]
utc.range.nom.m1.wn6 <- inla.emarginal(function(x) x, utc.range.nom.marg.wn6)
utc.range.nom.m2.wn6 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.wn6)
utc.range.nom.stdev.wn6 <- sqrt(utc.range.nom.m2.wn6 - utc.range.nom.m1.wn6^2)
utc.range.nom.quantiles.wn6 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.wn6)

### extract predicted values
utc.index.pred.response.wn6 <- inla.stack.index(utc.join.stack.pred.15as.wn6, "pred.response")$data 
utc.index.pred.latent.wn6 <- inla.stack.index(utc.join.stack.pred.15as.wn6, "pred.latent")$data 

utc.pred.response.mean.wn6 <- utc.output.pred.15as.wn6$summary.linear.predictor[utc.index.pred.response.wn6, "mean"]
utc.pred.response.sd.wn6 <- utc.output.pred.15as.wn6$summary.linear.predictor[utc.index.pred.response.wn6, "sd"]

utc.pred.latent.mean.wn6 <- utc.output.pred.15as.wn6$summary.linear.predictor[utc.index.pred.latent.wn6, "mean"]
utc.pred.latent.sd.wn6 <- utc.output.pred.15as.wn6$summary.linear.predictor[utc.index.pred.latent.wn6, "sd"]

utc_pred_resp_mean.wn6 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.wn6.mod/4+1))
utc_pred_resp_mean_df.wn6 <- as.data.frame(mapply(c,utc_pred_resp_mean.wn6, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.wn6 <- orderBy(~y + x, data=utc_pred_resp_mean_df.wn6)
utc.pred.response.mean.grid.wn6 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.wn6$`log(utc.pred.response.mean.wn6.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.wn6 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.wn6)
utc_pred_resp_sd_df.wn6 <- as.data.frame(mapply(c,utc_pred_resp_sd.wn6, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.wn6 <- orderBy(~y + x, data=utc_pred_resp_sd_df.wn6)
utc.pred.response.sd.grid.wn6 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.wn6$utc.pred.response.sd.wn6)), nrow=utc.x.res*2)

utc_pred_lat_mean.wn6 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.wn6)
utc_pred_lat_mean_df.wn6 <- as.data.frame(mapply(c,utc_pred_lat_mean.wn6, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.wn6 <- orderBy(~y + x, data=utc_pred_lat_mean_df.wn6)
utc.pred.latent.mean.grid.wn6 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.wn6$utc.pred.latent.mean.wn6)), nrow=utc.x.res*2)

utc_pred_lat_sd.wn6 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.wn6)
utc_pred_lat_sd_df.wn6 <- as.data.frame(mapply(c,utc_pred_lat_sd.wn6, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.wn6 <- orderBy(~y + x, data=utc_pred_lat_sd_df.wn6)
utc.pred.latent.sd.grid.wn6 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.wn6$utc.pred.latent.sd.wn6)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wn6),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wn6), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wn6-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wn6-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wn6-utc.pred.latent.mean.grid.wn6)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wn6-utc.pred.latent.mean.grid.wn6)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wn6,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wn6, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wn6,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wn6, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.wn6.mod <- exp(utc.pred.response.mean.wn6)-1
utc.pred.response.mean.wn6.mod[utc.pred.response.mean.wn6.mod<0]<-0
res.15as.wn6<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.wn6.mod/4
RMSE.15as.wn6 <- sqrt(mean(res.15as.wn6^2))
RMSE.15as.wn6 
correl.15as.wn6 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.wn6.mod/4) 
correl.15as.wn6 

####################################
# real case sampled scenario + worst prior
# full model w/ log transformation + pc prior est at 30 arcsec full data
# est at 3 arcsec sampled 1% data (weighted) + pred at 3 arcsec
####################################
utc.stack.est.15as.wp1 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_w1$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_w1, 1),
                                     effects = list(c(utc.30as.s.index.pc, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_w1$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_w1$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.wp1 <- inla.stack(utc.stack.est.15as.wp1, utc.stack.latent.wp.p, utc.stack.pred.15as.wp.p)

utc.output.pred.15as.wp1 <- inla(utc.formula.log.30as.pc, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.wp1, spde=utc.30as.spde.pc), family="gaussian",
                                 control.fixed=utc.30as.fixed.normal, 
                                 control.family=list(hyper=list(prec=utc.30as.prec.pc)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.wp1), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.wp1)

utc.output.pred.15as.wp1$dic$dic
utc.output.pred.15as.wp1$waic$waic
slcpo(utc.output.pred.15as.wp1)

### extract parameter estimates
utc.fixed.out.wp1 <- round(utc.output.pred.15as.wp1$summary.fixed[,1:5],3) 
utc.fixed1_marg.wp1  <- utc.output.pred.15as.wp1$marginals.fixed[[1]]
utc.fixed2_marg.wp1  <- utc.output.pred.15as.wp1$marginals.fixed[[2]]
utc.fixed3_marg.wp1  <- utc.output.pred.15as.wp1$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.wp1 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.wp1$marginals.hyperpar[[1]])
utc.sigma2e_m1.wp1 <- inla.emarginal(function(x) x, utc.sigma2e_marg.wp1)
utc.sigma2e_m2.wp1 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.wp1)
utc.sigma2e_stdev.wp1 <- sqrt(utc.sigma2e_m2.wp1 - utc.sigma2e_m1.wp1^2)
utc.sigma2e_quantiles.wp1 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.wp1)

utc.spde.result.wp1 <- inla.spde2.result(utc.output.pred.15as.wp1, name="utc.30as.spatial.field.pc", utc.30as.spde.pc)

# variance of spatial effects
utc.var.nom.marg.wp1 <- utc.spde.result.wp1$marginals.variance.nominal[[1]]
utc.var.nom.m1.wp1 <- inla.emarginal(function(x) x, utc.var.nom.marg.wp1)
utc.var.nom.m2.wp1 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.wp1)
utc.var.nom.stdev.wp1 <- sqrt(utc.var.nom.m2.wp1 - utc.var.nom.m1.wp1^2)
utc.var.nom.quantiles.wp1 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.wp1)

utc.range.nom.marg.wp1 <- utc.spde.result.wp1$marginals.range.nominal[[1]]
utc.range.nom.m1.wp1 <- inla.emarginal(function(x) x, utc.range.nom.marg.wp1)
utc.range.nom.m2.wp1 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.wp1)
utc.range.nom.stdev.wp1 <- sqrt(utc.range.nom.m2.wp1 - utc.range.nom.m1.wp1^2)
utc.range.nom.quantiles.wp1 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.wp1)

### extract predicted values
utc.index.pred.response.wp1 <- inla.stack.index(utc.join.stack.pred.15as.wp1, "pred.response")$data 
utc.index.pred.latent.wp1 <- inla.stack.index(utc.join.stack.pred.15as.wp1, "pred.latent")$data 

utc.pred.response.mean.wp1 <- utc.output.pred.15as.wp1$summary.linear.predictor[utc.index.pred.response.wp1, "mean"]
utc.pred.response.sd.wp1 <- utc.output.pred.15as.wp1$summary.linear.predictor[utc.index.pred.response.wp1, "sd"]

utc.pred.latent.mean.wp1 <- utc.output.pred.15as.wp1$summary.linear.predictor[utc.index.pred.latent.wp1, "mean"]
utc.pred.latent.sd.wp1 <- utc.output.pred.15as.wp1$summary.linear.predictor[utc.index.pred.latent.wp1, "sd"]

utc_pred_resp_mean.wp1 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.wp1.mod/4+1))
utc_pred_resp_mean_df.wp1 <- as.data.frame(mapply(c,utc_pred_resp_mean.wp1, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.wp1 <- orderBy(~y + x, data=utc_pred_resp_mean_df.wp1)
utc.pred.response.mean.grid.wp1 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.wp1$`log(utc.pred.response.mean.wp1.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.wp1 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.wp1)
utc_pred_resp_sd_df.wp1 <- as.data.frame(mapply(c,utc_pred_resp_sd.wp1, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.wp1 <- orderBy(~y + x, data=utc_pred_resp_sd_df.wp1)
utc.pred.response.sd.grid.wp1 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.wp1$utc.pred.response.sd.wp1)), nrow=utc.x.res*2)

utc_pred_lat_mean.wp1 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.wp1)
utc_pred_lat_mean_df.wp1 <- as.data.frame(mapply(c,utc_pred_lat_mean.wp1, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.wp1 <- orderBy(~y + x, data=utc_pred_lat_mean_df.wp1)
utc.pred.latent.mean.grid.wp1 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.wp1$utc.pred.latent.mean.wp1)), nrow=utc.x.res*2)

utc_pred_lat_sd.wp1 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.wp1)
utc_pred_lat_sd_df.wp1 <- as.data.frame(mapply(c,utc_pred_lat_sd.wp1, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.wp1 <- orderBy(~y + x, data=utc_pred_lat_sd_df.wp1)
utc.pred.latent.sd.grid.wp1 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.wp1$utc.pred.latent.sd.wp1)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wp1),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wp1), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wp1-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wp1-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wp1-utc.pred.latent.mean.grid.wp1)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wp1-utc.pred.latent.mean.grid.wp1)-log(100)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wp1,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wp1, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wp1,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wp1, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.wp1.mod <- exp(utc.pred.response.mean.wp1)-1
utc.pred.response.mean.wp1.mod[utc.pred.response.mean.wp1.mod<0]<-0
res.15as.wp1<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.wp1.mod/4
RMSE.15as.wp1 <- sqrt(mean(res.15as.wp1^2))
RMSE.15as.wp1 
correl.15as.wp1 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.wp1.mod/4) 
correl.15as.wp1 


####################################
# real case sampled scenario + worst prior
# full model w/ log transformation + pc prior est at 30 arcsec full data
# est at 3 arcsec sampled 2% data (weighted) + pred at 3 arcsec
####################################
utc.stack.est.15as.wp2 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_w2$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_w2, 1),
                                     effects = list(c(utc.30as.s.index.pc, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_w2$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_w2$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.wp2 <- inla.stack(utc.stack.est.15as.wp2, utc.stack.latent.wp.p, utc.stack.pred.15as.wp.p)

utc.output.pred.15as.wp2 <- inla(utc.formula.log.30as.pc, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.wp2, spde=utc.30as.spde.pc), family="gaussian",
                                 control.fixed=utc.30as.fixed.normal, 
                                 control.family=list(hyper=list(prec=utc.30as.prec.pc)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.wp2), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.wp2)

utc.output.pred.15as.wp2$dic$dic
utc.output.pred.15as.wp2$waic$waic
slcpo(utc.output.pred.15as.wp2)

### extract parameter estimates
utc.fixed.out.wp2 <- round(utc.output.pred.15as.wp2$summary.fixed[,1:5],3) 
utc.fixed1_marg.wp2  <- utc.output.pred.15as.wp2$marginals.fixed[[1]]
utc.fixed2_marg.wp2  <- utc.output.pred.15as.wp2$marginals.fixed[[2]]
utc.fixed3_marg.wp2  <- utc.output.pred.15as.wp2$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.wp2 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.wp2$marginals.hyperpar[[1]])
utc.sigma2e_m1.wp2 <- inla.emarginal(function(x) x, utc.sigma2e_marg.wp2)
utc.sigma2e_m2.wp2 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.wp2)
utc.sigma2e_stdev.wp2 <- sqrt(utc.sigma2e_m2.wp2 - utc.sigma2e_m1.wp2^2)
utc.sigma2e_quantiles.wp2 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.wp2)

utc.spde.result.wp2 <- inla.spde2.result(utc.output.pred.15as.wp2, name="utc.30as.spatial.field.pc", utc.30as.spde.pc)

# variance of spatial effects
utc.var.nom.marg.wp2 <- utc.spde.result.wp2$marginals.variance.nominal[[1]]
utc.var.nom.m1.wp2 <- inla.emarginal(function(x) x, utc.var.nom.marg.wp2)
utc.var.nom.m2.wp2 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.wp2)
utc.var.nom.stdev.wp2 <- sqrt(utc.var.nom.m2.wp2 - utc.var.nom.m1.wp2^2)
utc.var.nom.quantiles.wp2 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.wp2)

utc.range.nom.marg.wp2 <- utc.spde.result.wp2$marginals.range.nominal[[1]]
utc.range.nom.m1.wp2 <- inla.emarginal(function(x) x, utc.range.nom.marg.wp2)
utc.range.nom.m2.wp2 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.wp2)
utc.range.nom.stdev.wp2 <- sqrt(utc.range.nom.m2.wp2 - utc.range.nom.m1.wp2^2)
utc.range.nom.quantiles.wp2 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.wp2)

### extract predicted values
utc.index.pred.response.wp2 <- inla.stack.index(utc.join.stack.pred.15as.wp2, "pred.response")$data 
utc.index.pred.latent.wp2 <- inla.stack.index(utc.join.stack.pred.15as.wp2, "pred.latent")$data 

utc.pred.response.mean.wp2 <- utc.output.pred.15as.wp2$summary.linear.predictor[utc.index.pred.response.wp2, "mean"]
utc.pred.response.sd.wp2 <- utc.output.pred.15as.wp2$summary.linear.predictor[utc.index.pred.response.wp2, "sd"]

utc.pred.latent.mean.wp2 <- utc.output.pred.15as.wp2$summary.linear.predictor[utc.index.pred.latent.wp2, "mean"]
utc.pred.latent.sd.wp2 <- utc.output.pred.15as.wp2$summary.linear.predictor[utc.index.pred.latent.wp2, "sd"]

utc_pred_resp_mean.wp2 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.wp2.mod/4+1))
utc_pred_resp_mean_df.wp2 <- as.data.frame(mapply(c,utc_pred_resp_mean.wp2, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.wp2 <- orderBy(~y + x, data=utc_pred_resp_mean_df.wp2)
utc.pred.response.mean.grid.wp2 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.wp2$`log(utc.pred.response.mean.wp2.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.wp2 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.wp2)
utc_pred_resp_sd_df.wp2 <- as.data.frame(mapply(c,utc_pred_resp_sd.wp2, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.wp2 <- orderBy(~y + x, data=utc_pred_resp_sd_df.wp2)
utc.pred.response.sd.grid.wp2 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.wp2$utc.pred.response.sd.wp2)), nrow=utc.x.res*2)

utc_pred_lat_mean.wp2 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.wp2)
utc_pred_lat_mean_df.wp2 <- as.data.frame(mapply(c,utc_pred_lat_mean.wp2, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.wp2 <- orderBy(~y + x, data=utc_pred_lat_mean_df.wp2)
utc.pred.latent.mean.grid.wp2 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.wp2$utc.pred.latent.mean.wp2)), nrow=utc.x.res*2)

utc_pred_lat_sd.wp2 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.wp2)
utc_pred_lat_sd_df.wp2 <- as.data.frame(mapply(c,utc_pred_lat_sd.wp2, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.wp2 <- orderBy(~y + x, data=utc_pred_lat_sd_df.wp2)
utc.pred.latent.sd.grid.wp2 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.wp2$utc.pred.latent.sd.wp2)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wp2),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wp2), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wp2-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wp2-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wp2-utc.pred.latent.mean.grid.wp2)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wp2-utc.pred.latent.mean.grid.wp2)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wp2,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wp2, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wp2,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wp2, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.wp2.mod <- exp(utc.pred.response.mean.wp2)-1
utc.pred.response.mean.wp2.mod[utc.pred.response.mean.wp2.mod<0]<-0
res.15as.wp2<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.wp2.mod/4
RMSE.15as.wp2 <- sqrt(mean(res.15as.wp2^2))
RMSE.15as.wp2 
correl.15as.wp2 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.wp2.mod/4) 
correl.15as.wp2 

####################################
# real case sampled scenario + worst prior
# full model w/ log transformation + pc prior est at 30 arcsec full data
# est at 3 arcsec sampled 5% data (weighted) + pred at 3 arcsec
####################################
utc.stack.est.15as.wp3 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_w3$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_w3, 1),
                                     effects = list(c(utc.30as.s.index.pc, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_w3$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_w3$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.wp3 <- inla.stack(utc.stack.est.15as.wp3, utc.stack.latent.wp.p, utc.stack.pred.15as.wp.p)

utc.output.pred.15as.wp3 <- inla(utc.formula.log.30as.pc, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.wp3, spde=utc.30as.spde.pc), family="gaussian",
                                 control.fixed=utc.30as.fixed.normal, 
                                 control.family=list(hyper=list(prec=utc.30as.prec.pc)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.wp3), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.wp3)

utc.output.pred.15as.wp3$dic$dic
utc.output.pred.15as.wp3$waic$waic
slcpo(utc.output.pred.15as.wp3)

### extract parameter estimates
utc.fixed.out.wp3 <- round(utc.output.pred.15as.wp3$summary.fixed[,1:5],3) 
utc.fixed1_marg.wp3  <- utc.output.pred.15as.wp3$marginals.fixed[[1]]
utc.fixed2_marg.wp3  <- utc.output.pred.15as.wp3$marginals.fixed[[2]]
utc.fixed3_marg.wp3  <- utc.output.pred.15as.wp3$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.wp3 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.wp3$marginals.hyperpar[[1]])
utc.sigma2e_m1.wp3 <- inla.emarginal(function(x) x, utc.sigma2e_marg.wp3)
utc.sigma2e_m2.wp3 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.wp3)
utc.sigma2e_stdev.wp3 <- sqrt(utc.sigma2e_m2.wp3 - utc.sigma2e_m1.wp3^2)
utc.sigma2e_quantiles.wp3 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.wp3)

utc.spde.result.wp3 <- inla.spde2.result(utc.output.pred.15as.wp3, name="utc.30as.spatial.field.pc", utc.30as.spde.pc)

# variance of spatial effects
utc.var.nom.marg.wp3 <- utc.spde.result.wp3$marginals.variance.nominal[[1]]
utc.var.nom.m1.wp3 <- inla.emarginal(function(x) x, utc.var.nom.marg.wp3)
utc.var.nom.m2.wp3 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.wp3)
utc.var.nom.stdev.wp3 <- sqrt(utc.var.nom.m2.wp3 - utc.var.nom.m1.wp3^2)
utc.var.nom.quantiles.wp3 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.wp3)

utc.range.nom.marg.wp3 <- utc.spde.result.wp3$marginals.range.nominal[[1]]
utc.range.nom.m1.wp3 <- inla.emarginal(function(x) x, utc.range.nom.marg.wp3)
utc.range.nom.m2.wp3 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.wp3)
utc.range.nom.stdev.wp3 <- sqrt(utc.range.nom.m2.wp3 - utc.range.nom.m1.wp3^2)
utc.range.nom.quantiles.wp3 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.wp3)

### extract predicted values
utc.index.pred.response.wp3 <- inla.stack.index(utc.join.stack.pred.15as.wp3, "pred.response")$data 
utc.index.pred.latent.wp3 <- inla.stack.index(utc.join.stack.pred.15as.wp3, "pred.latent")$data 

utc.pred.response.mean.wp3 <- utc.output.pred.15as.wp3$summary.linear.predictor[utc.index.pred.response.wp3, "mean"]
utc.pred.response.sd.wp3 <- utc.output.pred.15as.wp3$summary.linear.predictor[utc.index.pred.response.wp3, "sd"]

utc.pred.latent.mean.wp3 <- utc.output.pred.15as.wp3$summary.linear.predictor[utc.index.pred.latent.wp3, "mean"]
utc.pred.latent.sd.wp3 <- utc.output.pred.15as.wp3$summary.linear.predictor[utc.index.pred.latent.wp3, "sd"]

utc_pred_resp_mean.wp3 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.wp3.mod/4+1))
utc_pred_resp_mean_df.wp3 <- as.data.frame(mapply(c,utc_pred_resp_mean.wp3, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.wp3 <- orderBy(~y + x, data=utc_pred_resp_mean_df.wp3)
utc.pred.response.mean.grid.wp3 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.wp3$`log(utc.pred.response.mean.wp3.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.wp3 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.wp3)
utc_pred_resp_sd_df.wp3 <- as.data.frame(mapply(c,utc_pred_resp_sd.wp3, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.wp3 <- orderBy(~y + x, data=utc_pred_resp_sd_df.wp3)
utc.pred.response.sd.grid.wp3 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.wp3$utc.pred.response.sd.wp3)), nrow=utc.x.res*2)

utc_pred_lat_mean.wp3 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.wp3)
utc_pred_lat_mean_df.wp3 <- as.data.frame(mapply(c,utc_pred_lat_mean.wp3, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.wp3 <- orderBy(~y + x, data=utc_pred_lat_mean_df.wp3)
utc.pred.latent.mean.grid.wp3 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.wp3$utc.pred.latent.mean.wp3)), nrow=utc.x.res*2)

utc_pred_lat_sd.wp3 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.wp3)
utc_pred_lat_sd_df.wp3 <- as.data.frame(mapply(c,utc_pred_lat_sd.wp3, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.wp3 <- orderBy(~y + x, data=utc_pred_lat_sd_df.wp3)
utc.pred.latent.sd.grid.wp3 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.wp3$utc.pred.latent.sd.wp3)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wp3),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wp3), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wp3-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wp3-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wp3-utc.pred.latent.mean.grid.wp3)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wp3-utc.pred.latent.mean.grid.wp3)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wp3,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wp3, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wp3,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wp3, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.wp3.mod <- exp(utc.pred.response.mean.wp3)-1
utc.pred.response.mean.wp3.mod[utc.pred.response.mean.wp3.mod<0]<-0
res.15as.wp3<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.wp3.mod/4
RMSE.15as.wp3 <- sqrt(mean(res.15as.wp3^2))
RMSE.15as.wp3 
correl.15as.wp3 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.wp3.mod/4) 
correl.15as.wp3 

####################################
# real case sampled scenario + worst prior
# full model w/ log transformation + pc prior est at 30 arcsec full data
# est at 3 arcsec sampled 10% data (weighted) + pred at 3 arcsec
####################################
utc.stack.est.15as.wp4 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_w4$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_w4, 1),
                                     effects = list(c(utc.30as.s.index.pc, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_w4$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_w4$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.wp4 <- inla.stack(utc.stack.est.15as.wp4, utc.stack.latent.wp.p, utc.stack.pred.15as.wp.p)

utc.output.pred.15as.wp4 <- inla(utc.formula.log.30as.pc, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.wp4, spde=utc.30as.spde.pc), family="gaussian",
                                 control.fixed=utc.30as.fixed.normal, 
                                 control.family=list(hyper=list(prec=utc.30as.prec.pc)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.wp4), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.wp4)

utc.output.pred.15as.wp4$dic$dic
utc.output.pred.15as.wp4$waic$waic
slcpo(utc.output.pred.15as.wp4)

### extract parameter estimates
utc.fixed.out.wp4 <- round(utc.output.pred.15as.wp4$summary.fixed[,1:5],3) 
utc.fixed1_marg.wp4  <- utc.output.pred.15as.wp4$marginals.fixed[[1]]
utc.fixed2_marg.wp4  <- utc.output.pred.15as.wp4$marginals.fixed[[2]]
utc.fixed3_marg.wp4  <- utc.output.pred.15as.wp4$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.wp4 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.wp4$marginals.hyperpar[[1]])
utc.sigma2e_m1.wp4 <- inla.emarginal(function(x) x, utc.sigma2e_marg.wp4)
utc.sigma2e_m2.wp4 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.wp4)
utc.sigma2e_stdev.wp4 <- sqrt(utc.sigma2e_m2.wp4 - utc.sigma2e_m1.wp4^2)
utc.sigma2e_quantiles.wp4 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.wp4)

utc.spde.result.wp4 <- inla.spde2.result(utc.output.pred.15as.wp4, name="utc.30as.spatial.field.pc", utc.30as.spde.pc)

# variance of spatial effects
utc.var.nom.marg.wp4 <- utc.spde.result.wp4$marginals.variance.nominal[[1]]
utc.var.nom.m1.wp4 <- inla.emarginal(function(x) x, utc.var.nom.marg.wp4)
utc.var.nom.m2.wp4 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.wp4)
utc.var.nom.stdev.wp4 <- sqrt(utc.var.nom.m2.wp4 - utc.var.nom.m1.wp4^2)
utc.var.nom.quantiles.wp4 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.wp4)

utc.range.nom.marg.wp4 <- utc.spde.result.wp4$marginals.range.nominal[[1]]
utc.range.nom.m1.wp4 <- inla.emarginal(function(x) x, utc.range.nom.marg.wp4)
utc.range.nom.m2.wp4 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.wp4)
utc.range.nom.stdev.wp4 <- sqrt(utc.range.nom.m2.wp4 - utc.range.nom.m1.wp4^2)
utc.range.nom.quantiles.wp4 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.wp4)

### extract predicted values
utc.index.pred.response.wp4 <- inla.stack.index(utc.join.stack.pred.15as.wp4, "pred.response")$data 
utc.index.pred.latent.wp4 <- inla.stack.index(utc.join.stack.pred.15as.wp4, "pred.latent")$data 

utc.pred.response.mean.wp4 <- utc.output.pred.15as.wp4$summary.linear.predictor[utc.index.pred.response.wp4, "mean"]
utc.pred.response.sd.wp4 <- utc.output.pred.15as.wp4$summary.linear.predictor[utc.index.pred.response.wp4, "sd"]

utc.pred.latent.mean.wp4 <- utc.output.pred.15as.wp4$summary.linear.predictor[utc.index.pred.latent.wp4, "mean"]
utc.pred.latent.sd.wp4 <- utc.output.pred.15as.wp4$summary.linear.predictor[utc.index.pred.latent.wp4, "sd"]

utc_pred_resp_mean.wp4 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.wp4.mod/4+1))
utc_pred_resp_mean_df.wp4 <- as.data.frame(mapply(c,utc_pred_resp_mean.wp4, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.wp4 <- orderBy(~y + x, data=utc_pred_resp_mean_df.wp4)
utc.pred.response.mean.grid.wp4 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.wp4$`log(utc.pred.response.mean.wp4.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.wp4 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.wp4)
utc_pred_resp_sd_df.wp4 <- as.data.frame(mapply(c,utc_pred_resp_sd.wp4, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.wp4 <- orderBy(~y + x, data=utc_pred_resp_sd_df.wp4)
utc.pred.response.sd.grid.wp4 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.wp4$utc.pred.response.sd.wp4)), nrow=utc.x.res*2)

utc_pred_lat_mean.wp4 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.wp4)
utc_pred_lat_mean_df.wp4 <- as.data.frame(mapply(c,utc_pred_lat_mean.wp4, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.wp4 <- orderBy(~y + x, data=utc_pred_lat_mean_df.wp4)
utc.pred.latent.mean.grid.wp4 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.wp4$utc.pred.latent.mean.wp4)), nrow=utc.x.res*2)

utc_pred_lat_sd.wp4 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.wp4)
utc_pred_lat_sd_df.wp4 <- as.data.frame(mapply(c,utc_pred_lat_sd.wp4, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.wp4 <- orderBy(~y + x, data=utc_pred_lat_sd_df.wp4)
utc.pred.latent.sd.grid.wp4 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.wp4$utc.pred.latent.sd.wp4)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wp4),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wp4), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wp4-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wp4-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wp4-utc.pred.latent.mean.grid.wp4)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wp4-utc.pred.latent.mean.grid.wp4)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wp4,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wp4, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wp4,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wp4, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.wp4.mod <- exp(utc.pred.response.mean.wp4)-1
utc.pred.response.mean.wp4.mod[utc.pred.response.mean.wp4.mod<0]<-0
res.15as.wp4<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.wp4.mod/4
RMSE.15as.wp4 <- sqrt(mean(res.15as.wp4^2))
RMSE.15as.wp4 
correl.15as.wp4 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.wp4.mod/4) 
correl.15as.wp4 


####################################
# real case sampled scenario + worst prior
# full model w/ log transformation + pc prior est at 30 arcsec full data
# est at 3 arcsec sampled 20% data (weighted) + pred at 3 arcsec
####################################
utc.stack.est.15as.wp5 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_w5$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_w5, 1),
                                     effects = list(c(utc.30as.s.index.pc, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_w5$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_w5$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.wp5 <- inla.stack(utc.stack.est.15as.wp5, utc.stack.latent.wp.p, utc.stack.pred.15as.wp.p)

utc.output.pred.15as.wp5 <- inla(utc.formula.log.30as.pc, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.wp5, spde=utc.30as.spde.pc), family="gaussian",
                                 control.fixed=utc.30as.fixed.normal, 
                                 control.family=list(hyper=list(prec=utc.30as.prec.pc)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.wp5), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.wp5)

utc.output.pred.15as.wp5$dic$dic
utc.output.pred.15as.wp5$waic$waic
slcpo(utc.output.pred.15as.wp5)

### extract parameter estimates
utc.fixed.out.wp5 <- round(utc.output.pred.15as.wp5$summary.fixed[,1:5],3) 
utc.fixed1_marg.wp5  <- utc.output.pred.15as.wp5$marginals.fixed[[1]]
utc.fixed2_marg.wp5  <- utc.output.pred.15as.wp5$marginals.fixed[[2]]
utc.fixed3_marg.wp5  <- utc.output.pred.15as.wp5$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.wp5 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.wp5$marginals.hyperpar[[1]])
utc.sigma2e_m1.wp5 <- inla.emarginal(function(x) x, utc.sigma2e_marg.wp5)
utc.sigma2e_m2.wp5 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.wp5)
utc.sigma2e_stdev.wp5 <- sqrt(utc.sigma2e_m2.wp5 - utc.sigma2e_m1.wp5^2)
utc.sigma2e_quantiles.wp5 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.wp5)

utc.spde.result.wp5 <- inla.spde2.result(utc.output.pred.15as.wp5, name="utc.30as.spatial.field.pc", utc.30as.spde.pc)

# variance of spatial effects
utc.var.nom.marg.wp5 <- utc.spde.result.wp5$marginals.variance.nominal[[1]]
utc.var.nom.m1.wp5 <- inla.emarginal(function(x) x, utc.var.nom.marg.wp5)
utc.var.nom.m2.wp5 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.wp5)
utc.var.nom.stdev.wp5 <- sqrt(utc.var.nom.m2.wp5 - utc.var.nom.m1.wp5^2)
utc.var.nom.quantiles.wp5 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.wp5)

utc.range.nom.marg.wp5 <- utc.spde.result.wp5$marginals.range.nominal[[1]]
utc.range.nom.m1.wp5 <- inla.emarginal(function(x) x, utc.range.nom.marg.wp5)
utc.range.nom.m2.wp5 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.wp5)
utc.range.nom.stdev.wp5 <- sqrt(utc.range.nom.m2.wp5 - utc.range.nom.m1.wp5^2)
utc.range.nom.quantiles.wp5 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.wp5)

### extract predicted values
utc.index.pred.response.wp5 <- inla.stack.index(utc.join.stack.pred.15as.wp5, "pred.response")$data 
utc.index.pred.latent.wp5 <- inla.stack.index(utc.join.stack.pred.15as.wp5, "pred.latent")$data 

utc.pred.response.mean.wp5 <- utc.output.pred.15as.wp5$summary.linear.predictor[utc.index.pred.response.wp5, "mean"]
utc.pred.response.sd.wp5 <- utc.output.pred.15as.wp5$summary.linear.predictor[utc.index.pred.response.wp5, "sd"]

utc.pred.latent.mean.wp5 <- utc.output.pred.15as.wp5$summary.linear.predictor[utc.index.pred.latent.wp5, "mean"]
utc.pred.latent.sd.wp5 <- utc.output.pred.15as.wp5$summary.linear.predictor[utc.index.pred.latent.wp5, "sd"]

utc_pred_resp_mean.wp5 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.wp5.mod/4+1))
utc_pred_resp_mean_df.wp5 <- as.data.frame(mapply(c,utc_pred_resp_mean.wp5, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.wp5 <- orderBy(~y + x, data=utc_pred_resp_mean_df.wp5)
utc.pred.response.mean.grid.wp5 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.wp5$`log(utc.pred.response.mean.wp5.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.wp5 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.wp5)
utc_pred_resp_sd_df.wp5 <- as.data.frame(mapply(c,utc_pred_resp_sd.wp5, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.wp5 <- orderBy(~y + x, data=utc_pred_resp_sd_df.wp5)
utc.pred.response.sd.grid.wp5 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.wp5$utc.pred.response.sd.wp5)), nrow=utc.x.res*2)

utc_pred_lat_mean.wp5 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.wp5)
utc_pred_lat_mean_df.wp5 <- as.data.frame(mapply(c,utc_pred_lat_mean.wp5, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.wp5 <- orderBy(~y + x, data=utc_pred_lat_mean_df.wp5)
utc.pred.latent.mean.grid.wp5 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.wp5$utc.pred.latent.mean.wp5)), nrow=utc.x.res*2)

utc_pred_lat_sd.wp5 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.wp5)
utc_pred_lat_sd_df.wp5 <- as.data.frame(mapply(c,utc_pred_lat_sd.wp5, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.wp5 <- orderBy(~y + x, data=utc_pred_lat_sd_df.wp5)
utc.pred.latent.sd.grid.wp5 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.wp5$utc.pred.latent.sd.wp5)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wp5),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wp5), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wp5-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wp5-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wp5-utc.pred.latent.mean.grid.wp5)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wp5-utc.pred.latent.mean.grid.wp5)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wp5,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wp5, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wp5,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wp5, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.wp5.mod <- exp(utc.pred.response.mean.wp5)-1
utc.pred.response.mean.wp5.mod[utc.pred.response.mean.wp5.mod<0]<-0
res.15as.wp5<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.wp5.mod/4
RMSE.15as.wp5 <- sqrt(mean(res.15as.wp5^2))
RMSE.15as.wp5 
correl.15as.wp5 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.wp5.mod/4) 
correl.15as.wp5 


####################################
# real case sampled scenario + worst prior
# full model w/ log transformation + pc prior est at 30 arcsec full data
# est at 3 arcsec sampled 50% data (weighted) + pred at 3 arcsec
####################################
utc.stack.est.15as.wp6 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_w6$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_w6, 1),
                                     effects = list(c(utc.30as.s.index.pc, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_w6$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_w6$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.wp6 <- inla.stack(utc.stack.est.15as.wp6, utc.stack.latent.wp.p, utc.stack.pred.15as.wp.p)

utc.output.pred.15as.wp6 <- inla(utc.formula.log.30as.pc, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.wp6, spde=utc.30as.spde.pc), family="gaussian",
                                 control.fixed=utc.30as.fixed.normal, 
                                 control.family=list(hyper=list(prec=utc.30as.prec.pc)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.wp6), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.wp6)

utc.output.pred.15as.wp6$dic$dic
utc.output.pred.15as.wp6$waic$waic
slcpo(utc.output.pred.15as.wp6)

### extract parameter estimates
utc.fixed.out.wp6 <- round(utc.output.pred.15as.wp6$summary.fixed[,1:5],3) 
utc.fixed1_marg.wp6  <- utc.output.pred.15as.wp6$marginals.fixed[[1]]
utc.fixed2_marg.wp6  <- utc.output.pred.15as.wp6$marginals.fixed[[2]]
utc.fixed3_marg.wp6  <- utc.output.pred.15as.wp6$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.wp6 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.wp6$marginals.hyperpar[[1]])
utc.sigma2e_m1.wp6 <- inla.emarginal(function(x) x, utc.sigma2e_marg.wp6)
utc.sigma2e_m2.wp6 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.wp6)
utc.sigma2e_stdev.wp6 <- sqrt(utc.sigma2e_m2.wp6 - utc.sigma2e_m1.wp6^2)
utc.sigma2e_quantiles.wp6 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.wp6)

utc.spde.result.wp6 <- inla.spde2.result(utc.output.pred.15as.wp6, name="utc.30as.spatial.field.pc", utc.30as.spde.pc)

# variance of spatial effects
utc.var.nom.marg.wp6 <- utc.spde.result.wp6$marginals.variance.nominal[[1]]
utc.var.nom.m1.wp6 <- inla.emarginal(function(x) x, utc.var.nom.marg.wp6)
utc.var.nom.m2.wp6 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.wp6)
utc.var.nom.stdev.wp6 <- sqrt(utc.var.nom.m2.wp6 - utc.var.nom.m1.wp6^2)
utc.var.nom.quantiles.wp6 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.wp6)

utc.range.nom.marg.wp6 <- utc.spde.result.wp6$marginals.range.nominal[[1]]
utc.range.nom.m1.wp6 <- inla.emarginal(function(x) x, utc.range.nom.marg.wp6)
utc.range.nom.m2.wp6 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.wp6)
utc.range.nom.stdev.wp6 <- sqrt(utc.range.nom.m2.wp6 - utc.range.nom.m1.wp6^2)
utc.range.nom.quantiles.wp6 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.wp6)

### extract predicted values
utc.index.pred.response.wp6 <- inla.stack.index(utc.join.stack.pred.15as.wp6, "pred.response")$data 
utc.index.pred.latent.wp6 <- inla.stack.index(utc.join.stack.pred.15as.wp6, "pred.latent")$data 

utc.pred.response.mean.wp6 <- utc.output.pred.15as.wp6$summary.linear.predictor[utc.index.pred.response.wp6, "mean"]
utc.pred.response.sd.wp6 <- utc.output.pred.15as.wp6$summary.linear.predictor[utc.index.pred.response.wp6, "sd"]

utc.pred.latent.mean.wp6 <- utc.output.pred.15as.wp6$summary.linear.predictor[utc.index.pred.latent.wp6, "mean"]
utc.pred.latent.sd.wp6 <- utc.output.pred.15as.wp6$summary.linear.predictor[utc.index.pred.latent.wp6, "sd"]

utc_pred_resp_mean.wp6 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.wp6.mod/4+1))
utc_pred_resp_mean_df.wp6 <- as.data.frame(mapply(c,utc_pred_resp_mean.wp6, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.wp6 <- orderBy(~y + x, data=utc_pred_resp_mean_df.wp6)
utc.pred.response.mean.grid.wp6 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.wp6$`log(utc.pred.response.mean.wp6.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.wp6 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.wp6)
utc_pred_resp_sd_df.wp6 <- as.data.frame(mapply(c,utc_pred_resp_sd.wp6, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.wp6 <- orderBy(~y + x, data=utc_pred_resp_sd_df.wp6)
utc.pred.response.sd.grid.wp6 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.wp6$utc.pred.response.sd.wp6)), nrow=utc.x.res*2)

utc_pred_lat_mean.wp6 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.wp6)
utc_pred_lat_mean_df.wp6 <- as.data.frame(mapply(c,utc_pred_lat_mean.wp6, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.wp6 <- orderBy(~y + x, data=utc_pred_lat_mean_df.wp6)
utc.pred.latent.mean.grid.wp6 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.wp6$utc.pred.latent.mean.wp6)), nrow=utc.x.res*2)

utc_pred_lat_sd.wp6 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.wp6)
utc_pred_lat_sd_df.wp6 <- as.data.frame(mapply(c,utc_pred_lat_sd.wp6, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.wp6 <- orderBy(~y + x, data=utc_pred_lat_sd_df.wp6)
utc.pred.latent.sd.grid.wp6 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.wp6$utc.pred.latent.sd.wp6)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wp6),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wp6), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wp6-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wp6-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wp6-utc.pred.latent.mean.grid.wp6)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wp6-utc.pred.latent.mean.grid.wp6)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wp6,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wp6, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wp6,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wp6, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.wp6.mod <- exp(utc.pred.response.mean.wp6)-1
utc.pred.response.mean.wp6.mod[utc.pred.response.mean.wp6.mod<0]<-0
res.15as.wp6<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.wp6.mod/4
RMSE.15as.wp6 <- sqrt(mean(res.15as.wp6^2))
RMSE.15as.wp6 
correl.15as.wp6 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.wp6.mod/4) 
correl.15as.wp6 

####################################
# real case sampled scenario + vague prior
# full model w/ log transformation + vague prior
# est at 3 arcsec sampled 1% data (weighted) + pred at 3 arcsec
####################################
utc.stack.est.15as.wv1 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_w1$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_w1, 1),
                                     effects = list(c(utc.s.index, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_w1$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_w1$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.wv1 <- inla.stack(utc.stack.est.15as.wv1, utc.stack.latent, utc.stack.pred.15as.log)

utc.output.pred.15as.wv1 <- inla(utc.formula.log, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.wv1, spde=utc.spde), family="gaussian",
                                 control.fixed=utc.fixed.vague, 
                                 control.family=list(hyper=list(prec=utc.prec.vague)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.wv1), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.wv1)

utc.output.pred.15as.wv1$dic$dic
utc.output.pred.15as.wv1$waic$waic
slcpo(utc.output.pred.15as.wv1)

### extract parameter estimates
utc.fixed.out.wv1 <- round(utc.output.pred.15as.wv1$summary.fixed[,1:5],3) 
utc.fixed1_marg.wv1  <- utc.output.pred.15as.wv1$marginals.fixed[[1]]
utc.fixed2_marg.wv1  <- utc.output.pred.15as.wv1$marginals.fixed[[2]]
utc.fixed3_marg.wv1  <- utc.output.pred.15as.wv1$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.wv1 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.wv1$marginals.hyperpar[[1]])
utc.sigma2e_m1.wv1 <- inla.emarginal(function(x) x, utc.sigma2e_marg.wv1)
utc.sigma2e_m2.wv1 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.wv1)
utc.sigma2e_stdev.wv1 <- sqrt(utc.sigma2e_m2.wv1 - utc.sigma2e_m1.wv1^2)
utc.sigma2e_quantiles.wv1 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.wv1)

utc.spde.result.wv1 <- inla.spde2.result(utc.output.pred.15as.wv1, name="utc.spatial.field", utc.spde)

# variance of spatial effects
utc.var.nom.marg.wv1 <- utc.spde.result.wv1$marginals.variance.nominal[[1]]
utc.var.nom.m1.wv1 <- inla.emarginal(function(x) x, utc.var.nom.marg.wv1)
utc.var.nom.m2.wv1 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.wv1)
utc.var.nom.stdev.wv1 <- sqrt(utc.var.nom.m2.wv1 - utc.var.nom.m1.wv1^2)
utc.var.nom.quantiles.wv1 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.wv1)

utc.range.nom.marg.wv1 <- utc.spde.result.wv1$marginals.range.nominal[[1]]
utc.range.nom.m1.wv1 <- inla.emarginal(function(x) x, utc.range.nom.marg.wv1)
utc.range.nom.m2.wv1 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.wv1)
utc.range.nom.stdev.wv1 <- sqrt(utc.range.nom.m2.wv1 - utc.range.nom.m1.wv1^2)
utc.range.nom.quantiles.wv1 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.wv1)

### extract predicted values
utc.index.pred.response.wv1 <- inla.stack.index(utc.join.stack.pred.15as.wv1, "pred.response")$data 
utc.index.pred.latent.wv1 <- inla.stack.index(utc.join.stack.pred.15as.wv1, "pred.latent")$data 

utc.pred.response.mean.wv1 <- utc.output.pred.15as.wv1$summary.linear.predictor[utc.index.pred.response.wv1, "mean"]
utc.pred.response.sd.wv1 <- utc.output.pred.15as.wv1$summary.linear.predictor[utc.index.pred.response.wv1, "sd"]

utc.pred.latent.mean.wv1 <- utc.output.pred.15as.wv1$summary.linear.predictor[utc.index.pred.latent.wv1, "mean"]
utc.pred.latent.sd.wv1 <- utc.output.pred.15as.wv1$summary.linear.predictor[utc.index.pred.latent.wv1, "sd"]

utc_pred_resp_mean.wv1 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.wv1.mod/4+1))
utc_pred_resp_mean_df.wv1 <- as.data.frame(mapply(c,utc_pred_resp_mean.wv1, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.wv1 <- orderBy(~y + x, data=utc_pred_resp_mean_df.wv1)
utc.pred.response.mean.grid.wv1 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.wv1$`log(utc.pred.response.mean.wv1.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.wv1 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.wv1)
utc_pred_resp_sd_df.wv1 <- as.data.frame(mapply(c,utc_pred_resp_sd.wv1, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.wv1 <- orderBy(~y + x, data=utc_pred_resp_sd_df.wv1)
utc.pred.response.sd.grid.wv1 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.wv1$utc.pred.response.sd.wv1)), nrow=utc.x.res*2)

utc_pred_lat_mean.wv1 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.wv1)
utc_pred_lat_mean_df.wv1 <- as.data.frame(mapply(c,utc_pred_lat_mean.wv1, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.wv1 <- orderBy(~y + x, data=utc_pred_lat_mean_df.wv1)
utc.pred.latent.mean.grid.wv1 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.wv1$utc.pred.latent.mean.wv1)), nrow=utc.x.res*2)

utc_pred_lat_sd.wv1 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.wv1)
utc_pred_lat_sd_df.wv1 <- as.data.frame(mapply(c,utc_pred_lat_sd.wv1, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.wv1 <- orderBy(~y + x, data=utc_pred_lat_sd_df.wv1)
utc.pred.latent.sd.grid.wv1 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.wv1$utc.pred.latent.sd.wv1)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wv1),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wv1), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wv1-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wv1-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wv1-utc.pred.latent.mean.grid.wv1)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wv1-utc.pred.latent.mean.grid.wv1)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wv1,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wv1, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wv1,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wv1, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.wv1.mod <- exp(utc.pred.response.mean.wv1)-1
utc.pred.response.mean.wv1.mod[utc.pred.response.mean.wv1.mod<0]<-0
res.15as.wv1<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.wv1.mod/4
RMSE.15as.wv1 <- sqrt(mean(res.15as.wv1^2))
RMSE.15as.wv1 
correl.15as.wv1 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.wv1.mod/4) 
correl.15as.wv1 


####################################
# real case sampled scenario + vague prior
# full model w/ log transformation + vague prior
# est at 3 arcsec sampled 2% data (weighted) + pred at 3 arcsec
####################################
utc.stack.est.15as.wv2 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_w2$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_w2, 1),
                                     effects = list(c(utc.s.index, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_w2$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_w2$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.wv2 <- inla.stack(utc.stack.est.15as.wv2, utc.stack.latent, utc.stack.pred.15as.log)

utc.output.pred.15as.wv2 <- inla(utc.formula.log, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.wv2, spde=utc.spde), family="gaussian",
                                 control.fixed=utc.fixed.vague, 
                                 control.family=list(hyper=list(prec=utc.prec.vague)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.wv2), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.wv2)

utc.output.pred.15as.wv2$dic$dic
utc.output.pred.15as.wv2$waic$waic
slcpo(utc.output.pred.15as.wv2)

### extract parameter estimates
utc.fixed.out.wv2 <- round(utc.output.pred.15as.wv2$summary.fixed[,1:5],3) 
utc.fixed1_marg.wv2  <- utc.output.pred.15as.wv2$marginals.fixed[[1]]
utc.fixed2_marg.wv2  <- utc.output.pred.15as.wv2$marginals.fixed[[2]]
utc.fixed3_marg.wv2  <- utc.output.pred.15as.wv2$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.wv2 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.wv2$marginals.hyperpar[[1]])
utc.sigma2e_m1.wv2 <- inla.emarginal(function(x) x, utc.sigma2e_marg.wv2)
utc.sigma2e_m2.wv2 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.wv2)
utc.sigma2e_stdev.wv2 <- sqrt(utc.sigma2e_m2.wv2 - utc.sigma2e_m1.wv2^2)
utc.sigma2e_quantiles.wv2 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.wv2)

utc.spde.result.wv2 <- inla.spde2.result(utc.output.pred.15as.wv2, name="utc.spatial.field", utc.spde)

# variance of spatial effects
utc.var.nom.marg.wv2 <- utc.spde.result.wv2$marginals.variance.nominal[[1]]
utc.var.nom.m1.wv2 <- inla.emarginal(function(x) x, utc.var.nom.marg.wv2)
utc.var.nom.m2.wv2 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.wv2)
utc.var.nom.stdev.wv2 <- sqrt(utc.var.nom.m2.wv2 - utc.var.nom.m1.wv2^2)
utc.var.nom.quantiles.wv2 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.wv2)

utc.range.nom.marg.wv2 <- utc.spde.result.wv2$marginals.range.nominal[[1]]
utc.range.nom.m1.wv2 <- inla.emarginal(function(x) x, utc.range.nom.marg.wv2)
utc.range.nom.m2.wv2 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.wv2)
utc.range.nom.stdev.wv2 <- sqrt(utc.range.nom.m2.wv2 - utc.range.nom.m1.wv2^2)
utc.range.nom.quantiles.wv2 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.wv2)

### extract predicted values
utc.index.pred.response.wv2 <- inla.stack.index(utc.join.stack.pred.15as.wv2, "pred.response")$data 
utc.index.pred.latent.wv2 <- inla.stack.index(utc.join.stack.pred.15as.wv2, "pred.latent")$data 

utc.pred.response.mean.wv2 <- utc.output.pred.15as.wv2$summary.linear.predictor[utc.index.pred.response.wv2, "mean"]
utc.pred.response.sd.wv2 <- utc.output.pred.15as.wv2$summary.linear.predictor[utc.index.pred.response.wv2, "sd"]

utc.pred.latent.mean.wv2 <- utc.output.pred.15as.wv2$summary.linear.predictor[utc.index.pred.latent.wv2, "mean"]
utc.pred.latent.sd.wv2 <- utc.output.pred.15as.wv2$summary.linear.predictor[utc.index.pred.latent.wv2, "sd"]

utc_pred_resp_mean.wv2 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.wv2.mod/4+1))
utc_pred_resp_mean_df.wv2 <- as.data.frame(mapply(c,utc_pred_resp_mean.wv2, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.wv2 <- orderBy(~y + x, data=utc_pred_resp_mean_df.wv2)
utc.pred.response.mean.grid.wv2 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.wv2$`log(utc.pred.response.mean.wv2.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.wv2 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.wv2)
utc_pred_resp_sd_df.wv2 <- as.data.frame(mapply(c,utc_pred_resp_sd.wv2, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.wv2 <- orderBy(~y + x, data=utc_pred_resp_sd_df.wv2)
utc.pred.response.sd.grid.wv2 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.wv2$utc.pred.response.sd.wv2)), nrow=utc.x.res*2)

utc_pred_lat_mean.wv2 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.wv2)
utc_pred_lat_mean_df.wv2 <- as.data.frame(mapply(c,utc_pred_lat_mean.wv2, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.wv2 <- orderBy(~y + x, data=utc_pred_lat_mean_df.wv2)
utc.pred.latent.mean.grid.wv2 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.wv2$utc.pred.latent.mean.wv2)), nrow=utc.x.res*2)

utc_pred_lat_sd.wv2 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.wv2)
utc_pred_lat_sd_df.wv2 <- as.data.frame(mapply(c,utc_pred_lat_sd.wv2, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.wv2 <- orderBy(~y + x, data=utc_pred_lat_sd_df.wv2)
utc.pred.latent.sd.grid.wv2 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.wv2$utc.pred.latent.sd.wv2)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wv2),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wv2), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wv2-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wv2-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wv2-utc.pred.latent.mean.grid.wv2)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wv2-utc.pred.latent.mean.grid.wv2)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wv2,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wv2, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wv2,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wv2, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.wv2.mod <- exp(utc.pred.response.mean.wv2)-1
utc.pred.response.mean.wv2.mod[utc.pred.response.mean.wv2.mod<0]<-0
res.15as.wv2<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.wv2.mod/4
RMSE.15as.wv2 <- sqrt(mean(res.15as.wv2^2))
RMSE.15as.wv2 
correl.15as.wv2 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.wv2.mod/4) 
correl.15as.wv2 

####################################
# real case sampled scenario + vague prior
# full model w/ log transformation + vague prior
# est at 3 arcsec sampled 5% data (weighted) + pred at 3 arcsec
####################################
utc.stack.est.15as.wv3 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_w3$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_w3, 1),
                                     effects = list(c(utc.s.index, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_w3$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_w3$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.wv3 <- inla.stack(utc.stack.est.15as.wv3, utc.stack.latent, utc.stack.pred.15as.log)

utc.output.pred.15as.wv3 <- inla(utc.formula.log, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.wv3, spde=utc.spde), family="gaussian",
                                 control.fixed=utc.fixed.vague, 
                                 control.family=list(hyper=list(prec=utc.prec.vague)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.wv3), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.wv3)

utc.output.pred.15as.wv3$dic$dic
utc.output.pred.15as.wv3$waic$waic
slcpo(utc.output.pred.15as.wv3)

### extract parameter estimates
utc.fixed.out.wv3 <- round(utc.output.pred.15as.wv3$summary.fixed[,1:5],3) 
utc.fixed1_marg.wv3  <- utc.output.pred.15as.wv3$marginals.fixed[[1]]
utc.fixed2_marg.wv3  <- utc.output.pred.15as.wv3$marginals.fixed[[2]]
utc.fixed3_marg.wv3  <- utc.output.pred.15as.wv3$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.wv3 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.wv3$marginals.hyperpar[[1]])
utc.sigma2e_m1.wv3 <- inla.emarginal(function(x) x, utc.sigma2e_marg.wv3)
utc.sigma2e_m2.wv3 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.wv3)
utc.sigma2e_stdev.wv3 <- sqrt(utc.sigma2e_m2.wv3 - utc.sigma2e_m1.wv3^2)
utc.sigma2e_quantiles.wv3 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.wv3)

utc.spde.result.wv3 <- inla.spde2.result(utc.output.pred.15as.wv3, name="utc.spatial.field", utc.spde)

# variance of spatial effects
utc.var.nom.marg.wv3 <- utc.spde.result.wv3$marginals.variance.nominal[[1]]
utc.var.nom.m1.wv3 <- inla.emarginal(function(x) x, utc.var.nom.marg.wv3)
utc.var.nom.m2.wv3 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.wv3)
utc.var.nom.stdev.wv3 <- sqrt(utc.var.nom.m2.wv3 - utc.var.nom.m1.wv3^2)
utc.var.nom.quantiles.wv3 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.wv3)

utc.range.nom.marg.wv3 <- utc.spde.result.wv3$marginals.range.nominal[[1]]
utc.range.nom.m1.wv3 <- inla.emarginal(function(x) x, utc.range.nom.marg.wv3)
utc.range.nom.m2.wv3 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.wv3)
utc.range.nom.stdev.wv3 <- sqrt(utc.range.nom.m2.wv3 - utc.range.nom.m1.wv3^2)
utc.range.nom.quantiles.wv3 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.wv3)

### extract predicted values
utc.index.pred.response.wv3 <- inla.stack.index(utc.join.stack.pred.15as.wv3, "pred.response")$data 
utc.index.pred.latent.wv3 <- inla.stack.index(utc.join.stack.pred.15as.wv3, "pred.latent")$data 

utc.pred.response.mean.wv3 <- utc.output.pred.15as.wv3$summary.linear.predictor[utc.index.pred.response.wv3, "mean"]
utc.pred.response.sd.wv3 <- utc.output.pred.15as.wv3$summary.linear.predictor[utc.index.pred.response.wv3, "sd"]

utc.pred.latent.mean.wv3 <- utc.output.pred.15as.wv3$summary.linear.predictor[utc.index.pred.latent.wv3, "mean"]
utc.pred.latent.sd.wv3 <- utc.output.pred.15as.wv3$summary.linear.predictor[utc.index.pred.latent.wv3, "sd"]

utc_pred_resp_mean.wv3 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.wv3.mod/4+1))
utc_pred_resp_mean_df.wv3 <- as.data.frame(mapply(c,utc_pred_resp_mean.wv3, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.wv3 <- orderBy(~y + x, data=utc_pred_resp_mean_df.wv3)
utc.pred.response.mean.grid.wv3 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.wv3$`log(utc.pred.response.mean.wv3.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.wv3 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.wv3)
utc_pred_resp_sd_df.wv3 <- as.data.frame(mapply(c,utc_pred_resp_sd.wv3, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.wv3 <- orderBy(~y + x, data=utc_pred_resp_sd_df.wv3)
utc.pred.response.sd.grid.wv3 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.wv3$utc.pred.response.sd.wv3)), nrow=utc.x.res*2)

utc_pred_lat_mean.wv3 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.wv3)
utc_pred_lat_mean_df.wv3 <- as.data.frame(mapply(c,utc_pred_lat_mean.wv3, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.wv3 <- orderBy(~y + x, data=utc_pred_lat_mean_df.wv3)
utc.pred.latent.mean.grid.wv3 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.wv3$utc.pred.latent.mean.wv3)), nrow=utc.x.res*2)

utc_pred_lat_sd.wv3 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.wv3)
utc_pred_lat_sd_df.wv3 <- as.data.frame(mapply(c,utc_pred_lat_sd.wv3, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.wv3 <- orderBy(~y + x, data=utc_pred_lat_sd_df.wv3)
utc.pred.latent.sd.grid.wv3 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.wv3$utc.pred.latent.sd.wv3)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wv3),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wv3), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wv3-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wv3-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wv3-utc.pred.latent.mean.grid.wv3)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wv3-utc.pred.latent.mean.grid.wv3)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wv3,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wv3, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wv3,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wv3, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.wv3.mod <- exp(utc.pred.response.mean.wv3)-1
utc.pred.response.mean.wv3.mod[utc.pred.response.mean.wv3.mod<0]<-0
res.15as.wv3<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.wv3.mod/4
RMSE.15as.wv3 <- sqrt(mean(res.15as.wv3^2))
RMSE.15as.wv3 
correl.15as.wv3 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.wv3.mod/4) 
correl.15as.wv3 

####################################
# real case sampled scenario + vague prior
# full model w/ log transformation + vague prior
# est at 3 arcsec sampled 10% data (weighted) + pred at 3 arcsec
####################################
utc.stack.est.15as.wv4 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_w4$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_w4, 1),
                                     effects = list(c(utc.s.index, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_w4$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_w4$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.wv4 <- inla.stack(utc.stack.est.15as.wv4, utc.stack.latent, utc.stack.pred.15as.log)

utc.output.pred.15as.wv4 <- inla(utc.formula.log, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.wv4, spde=utc.spde), family="gaussian",
                                 control.fixed=utc.fixed.vague, 
                                 control.family=list(hyper=list(prec=utc.prec.vague)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.wv4), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.wv4)

utc.output.pred.15as.wv4$dic$dic
utc.output.pred.15as.wv4$waic$waic
slcpo(utc.output.pred.15as.wv4)

### extract parameter estimates
utc.fixed.out.wv4 <- round(utc.output.pred.15as.wv4$summary.fixed[,1:5],3) 
utc.fixed1_marg.wv4  <- utc.output.pred.15as.wv4$marginals.fixed[[1]]
utc.fixed2_marg.wv4  <- utc.output.pred.15as.wv4$marginals.fixed[[2]]
utc.fixed3_marg.wv4  <- utc.output.pred.15as.wv4$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.wv4 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.wv4$marginals.hyperpar[[1]])
utc.sigma2e_m1.wv4 <- inla.emarginal(function(x) x, utc.sigma2e_marg.wv4)
utc.sigma2e_m2.wv4 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.wv4)
utc.sigma2e_stdev.wv4 <- sqrt(utc.sigma2e_m2.wv4 - utc.sigma2e_m1.wv4^2)
utc.sigma2e_quantiles.wv4 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.wv4)

utc.spde.result.wv4 <- inla.spde2.result(utc.output.pred.15as.wv4, name="utc.spatial.field", utc.spde)

# variance of spatial effects
utc.var.nom.marg.wv4 <- utc.spde.result.wv4$marginals.variance.nominal[[1]]
utc.var.nom.m1.wv4 <- inla.emarginal(function(x) x, utc.var.nom.marg.wv4)
utc.var.nom.m2.wv4 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.wv4)
utc.var.nom.stdev.wv4 <- sqrt(utc.var.nom.m2.wv4 - utc.var.nom.m1.wv4^2)
utc.var.nom.quantiles.wv4 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.wv4)

utc.range.nom.marg.wv4 <- utc.spde.result.wv4$marginals.range.nominal[[1]]
utc.range.nom.m1.wv4 <- inla.emarginal(function(x) x, utc.range.nom.marg.wv4)
utc.range.nom.m2.wv4 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.wv4)
utc.range.nom.stdev.wv4 <- sqrt(utc.range.nom.m2.wv4 - utc.range.nom.m1.wv4^2)
utc.range.nom.quantiles.wv4 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.wv4)

### extract predicted values
utc.index.pred.response.wv4 <- inla.stack.index(utc.join.stack.pred.15as.wv4, "pred.response")$data 
utc.index.pred.latent.wv4 <- inla.stack.index(utc.join.stack.pred.15as.wv4, "pred.latent")$data 

utc.pred.response.mean.wv4 <- utc.output.pred.15as.wv4$summary.linear.predictor[utc.index.pred.response.wv4, "mean"]
utc.pred.response.sd.wv4 <- utc.output.pred.15as.wv4$summary.linear.predictor[utc.index.pred.response.wv4, "sd"]

utc.pred.latent.mean.wv4 <- utc.output.pred.15as.wv4$summary.linear.predictor[utc.index.pred.latent.wv4, "mean"]
utc.pred.latent.sd.wv4 <- utc.output.pred.15as.wv4$summary.linear.predictor[utc.index.pred.latent.wv4, "sd"]

utc_pred_resp_mean.wv4 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.wv4.mod/4+1))
utc_pred_resp_mean_df.wv4 <- as.data.frame(mapply(c,utc_pred_resp_mean.wv4, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.wv4 <- orderBy(~y + x, data=utc_pred_resp_mean_df.wv4)
utc.pred.response.mean.grid.wv4 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.wv4$`log(utc.pred.response.mean.wv4.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.wv4 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.wv4)
utc_pred_resp_sd_df.wv4 <- as.data.frame(mapply(c,utc_pred_resp_sd.wv4, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.wv4 <- orderBy(~y + x, data=utc_pred_resp_sd_df.wv4)
utc.pred.response.sd.grid.wv4 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.wv4$utc.pred.response.sd.wv4)), nrow=utc.x.res*2)

utc_pred_lat_mean.wv4 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.wv4)
utc_pred_lat_mean_df.wv4 <- as.data.frame(mapply(c,utc_pred_lat_mean.wv4, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.wv4 <- orderBy(~y + x, data=utc_pred_lat_mean_df.wv4)
utc.pred.latent.mean.grid.wv4 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.wv4$utc.pred.latent.mean.wv4)), nrow=utc.x.res*2)

utc_pred_lat_sd.wv4 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.wv4)
utc_pred_lat_sd_df.wv4 <- as.data.frame(mapply(c,utc_pred_lat_sd.wv4, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.wv4 <- orderBy(~y + x, data=utc_pred_lat_sd_df.wv4)
utc.pred.latent.sd.grid.wv4 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.wv4$utc.pred.latent.sd.wv4)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wv4),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wv4), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wv4-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wv4-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wv4-utc.pred.latent.mean.grid.wv4)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wv4-utc.pred.latent.mean.grid.wv4)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wv4,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wv4, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wv4,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wv4, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.wv4.mod <- exp(utc.pred.response.mean.wv4)-1
utc.pred.response.mean.wv4.mod[utc.pred.response.mean.wv4.mod<0]<-0
res.15as.wv4<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.wv4.mod/4
RMSE.15as.wv4 <- sqrt(mean(res.15as.wv4^2))
RMSE.15as.wv4 
correl.15as.wv4 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.wv4.mod/4) 
correl.15as.wv4 

####################################
# real case sampled scenario + vague prior
# full model w/ log transformation + vague prior
# est at 3 arcsec sampled 20% data (weighted) + pred at 3 arcsec
####################################
utc.stack.est.15as.wv5 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_w5$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_w5, 1),
                                     effects = list(c(utc.s.index, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_w5$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_w5$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.wv5 <- inla.stack(utc.stack.est.15as.wv5, utc.stack.latent, utc.stack.pred.15as.log)

utc.output.pred.15as.wv5 <- inla(utc.formula.log, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.wv5, spde=utc.spde), family="gaussian",
                                 control.fixed=utc.fixed.vague, 
                                 control.family=list(hyper=list(prec=utc.prec.vague)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.wv5), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.wv5)

utc.output.pred.15as.wv5$dic$dic
utc.output.pred.15as.wv5$waic$waic
slcpo(utc.output.pred.15as.wv5)

### extract parameter estimates
utc.fixed.out.wv5 <- round(utc.output.pred.15as.wv5$summary.fixed[,1:5],3) 
utc.fixed1_marg.wv5  <- utc.output.pred.15as.wv5$marginals.fixed[[1]]
utc.fixed2_marg.wv5  <- utc.output.pred.15as.wv5$marginals.fixed[[2]]
utc.fixed3_marg.wv5  <- utc.output.pred.15as.wv5$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.wv5 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.wv5$marginals.hyperpar[[1]])
utc.sigma2e_m1.wv5 <- inla.emarginal(function(x) x, utc.sigma2e_marg.wv5)
utc.sigma2e_m2.wv5 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.wv5)
utc.sigma2e_stdev.wv5 <- sqrt(utc.sigma2e_m2.wv5 - utc.sigma2e_m1.wv5^2)
utc.sigma2e_quantiles.wv5 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.wv5)

utc.spde.result.wv5 <- inla.spde2.result(utc.output.pred.15as.wv5, name="utc.spatial.field", utc.spde)

# variance of spatial effects
utc.var.nom.marg.wv5 <- utc.spde.result.wv5$marginals.variance.nominal[[1]]
utc.var.nom.m1.wv5 <- inla.emarginal(function(x) x, utc.var.nom.marg.wv5)
utc.var.nom.m2.wv5 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.wv5)
utc.var.nom.stdev.wv5 <- sqrt(utc.var.nom.m2.wv5 - utc.var.nom.m1.wv5^2)
utc.var.nom.quantiles.wv5 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.wv5)

utc.range.nom.marg.wv5 <- utc.spde.result.wv5$marginals.range.nominal[[1]]
utc.range.nom.m1.wv5 <- inla.emarginal(function(x) x, utc.range.nom.marg.wv5)
utc.range.nom.m2.wv5 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.wv5)
utc.range.nom.stdev.wv5 <- sqrt(utc.range.nom.m2.wv5 - utc.range.nom.m1.wv5^2)
utc.range.nom.quantiles.wv5 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.wv5)

### extract predicted values
utc.index.pred.response.wv5 <- inla.stack.index(utc.join.stack.pred.15as.wv5, "pred.response")$data 
utc.index.pred.latent.wv5 <- inla.stack.index(utc.join.stack.pred.15as.wv5, "pred.latent")$data 

utc.pred.response.mean.wv5 <- utc.output.pred.15as.wv5$summary.linear.predictor[utc.index.pred.response.wv5, "mean"]
utc.pred.response.sd.wv5 <- utc.output.pred.15as.wv5$summary.linear.predictor[utc.index.pred.response.wv5, "sd"]

utc.pred.latent.mean.wv5 <- utc.output.pred.15as.wv5$summary.linear.predictor[utc.index.pred.latent.wv5, "mean"]
utc.pred.latent.sd.wv5 <- utc.output.pred.15as.wv5$summary.linear.predictor[utc.index.pred.latent.wv5, "sd"]

utc_pred_resp_mean.wv5 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.wv5.mod/4+1))
utc_pred_resp_mean_df.wv5 <- as.data.frame(mapply(c,utc_pred_resp_mean.wv5, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.wv5 <- orderBy(~y + x, data=utc_pred_resp_mean_df.wv5)
utc.pred.response.mean.grid.wv5 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.wv5$`log(utc.pred.response.mean.wv5.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.wv5 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.wv5)
utc_pred_resp_sd_df.wv5 <- as.data.frame(mapply(c,utc_pred_resp_sd.wv5, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.wv5 <- orderBy(~y + x, data=utc_pred_resp_sd_df.wv5)
utc.pred.response.sd.grid.wv5 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.wv5$utc.pred.response.sd.wv5)), nrow=utc.x.res*2)

utc_pred_lat_mean.wv5 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.wv5)
utc_pred_lat_mean_df.wv5 <- as.data.frame(mapply(c,utc_pred_lat_mean.wv5, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.wv5 <- orderBy(~y + x, data=utc_pred_lat_mean_df.wv5)
utc.pred.latent.mean.grid.wv5 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.wv5$utc.pred.latent.mean.wv5)), nrow=utc.x.res*2)

utc_pred_lat_sd.wv5 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.wv5)
utc_pred_lat_sd_df.wv5 <- as.data.frame(mapply(c,utc_pred_lat_sd.wv5, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.wv5 <- orderBy(~y + x, data=utc_pred_lat_sd_df.wv5)
utc.pred.latent.sd.grid.wv5 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.wv5$utc.pred.latent.sd.wv5)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wv5),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wv5), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wv5-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wv5-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wv5-utc.pred.latent.mean.grid.wv5)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wv5-utc.pred.latent.mean.grid.wv5)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wv5,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wv5, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wv5,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wv5, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.wv5.mod <- exp(utc.pred.response.mean.wv5)-1
utc.pred.response.mean.wv5.mod[utc.pred.response.mean.wv5.mod<0]<-0
res.15as.wv5<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.wv5.mod/4
RMSE.15as.wv5 <- sqrt(mean(res.15as.wv5^2))
RMSE.15as.wv5 
correl.15as.wv5 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.wv5.mod/4) 
correl.15as.wv5 

####################################
# real case sampled scenario + vague prior
# full model w/ log transformation + vague prior
# est at 3 arcsec sampled 50% data (weighted) + pred at 3 arcsec
####################################
utc.stack.est.15as.wv6 <- inla.stack(data = list(logpop=log(as.numeric(stacked_utc_sampled_w6$pop_utc_15_15as)*4+1)),
                                     A = list(utc.A.sampled.15as_w6, 1),
                                     effects = list(c(utc.s.index, list(Intercept=1)), 
                                                    list(ntl=as.numeric(stacked_utc_sampled_w6$ntl_utc_15_15as), 
                                                         slope=as.numeric(stacked_utc_sampled_w6$slope_utc_15as)
                                                    )),
                                     tag="est")

utc.join.stack.pred.15as.wv6 <- inla.stack(utc.stack.est.15as.wv6, utc.stack.latent, utc.stack.pred.15as.log)

utc.output.pred.15as.wv6 <- inla(utc.formula.log, 
                                 data=inla.stack.data(utc.join.stack.pred.15as.wv6, spde=utc.spde), family="gaussian",
                                 control.fixed=utc.fixed.vague, 
                                 control.family=list(hyper=list(prec=utc.prec.vague)), 
                                 control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.wv6), compute=TRUE),
                                 control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.wv6)

utc.output.pred.15as.wv6$dic$dic
utc.output.pred.15as.wv6$waic$waic
slcpo(utc.output.pred.15as.wv6)

### extract parameter estimates
utc.fixed.out.wv6 <- round(utc.output.pred.15as.wv6$summary.fixed[,1:5],3) 
utc.fixed1_marg.wv6  <- utc.output.pred.15as.wv6$marginals.fixed[[1]]
utc.fixed2_marg.wv6  <- utc.output.pred.15as.wv6$marginals.fixed[[2]]
utc.fixed3_marg.wv6  <- utc.output.pred.15as.wv6$marginals.fixed[[3]]

# variance of unstructured residual
utc.sigma2e_marg.wv6 <- inla.tmarginal(function(x) 1/x, utc.output.pred.15as.wv6$marginals.hyperpar[[1]])
utc.sigma2e_m1.wv6 <- inla.emarginal(function(x) x, utc.sigma2e_marg.wv6)
utc.sigma2e_m2.wv6 <- inla.emarginal(function(x) x^2, utc.sigma2e_marg.wv6)
utc.sigma2e_stdev.wv6 <- sqrt(utc.sigma2e_m2.wv6 - utc.sigma2e_m1.wv6^2)
utc.sigma2e_quantiles.wv6 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.sigma2e_marg.wv6)

utc.spde.result.wv6 <- inla.spde2.result(utc.output.pred.15as.wv6, name="utc.spatial.field", utc.spde)

# variance of spatial effects
utc.var.nom.marg.wv6 <- utc.spde.result.wv6$marginals.variance.nominal[[1]]
utc.var.nom.m1.wv6 <- inla.emarginal(function(x) x, utc.var.nom.marg.wv6)
utc.var.nom.m2.wv6 <- inla.emarginal(function(x) x^2, utc.var.nom.marg.wv6)
utc.var.nom.stdev.wv6 <- sqrt(utc.var.nom.m2.wv6 - utc.var.nom.m1.wv6^2)
utc.var.nom.quantiles.wv6 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.var.nom.marg.wv6)

utc.range.nom.marg.wv6 <- utc.spde.result.wv6$marginals.range.nominal[[1]]
utc.range.nom.m1.wv6 <- inla.emarginal(function(x) x, utc.range.nom.marg.wv6)
utc.range.nom.m2.wv6 <- inla.emarginal(function(x) x^2, utc.range.nom.marg.wv6)
utc.range.nom.stdev.wv6 <- sqrt(utc.range.nom.m2.wv6 - utc.range.nom.m1.wv6^2)
utc.range.nom.quantiles.wv6 <- inla.qmarginal(c(0.025, 0.5, 0.975), utc.range.nom.marg.wv6)

### extract predicted values
utc.index.pred.response.wv6 <- inla.stack.index(utc.join.stack.pred.15as.wv6, "pred.response")$data 
utc.index.pred.latent.wv6 <- inla.stack.index(utc.join.stack.pred.15as.wv6, "pred.latent")$data 

utc.pred.response.mean.wv6 <- utc.output.pred.15as.wv6$summary.linear.predictor[utc.index.pred.response.wv6, "mean"]
utc.pred.response.sd.wv6 <- utc.output.pred.15as.wv6$summary.linear.predictor[utc.index.pred.response.wv6, "sd"]

utc.pred.latent.mean.wv6 <- utc.output.pred.15as.wv6$summary.linear.predictor[utc.index.pred.latent.wv6, "mean"]
utc.pred.latent.sd.wv6 <- utc.output.pred.15as.wv6$summary.linear.predictor[utc.index.pred.latent.wv6, "sd"]

utc_pred_resp_mean.wv6 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],log(utc.pred.response.mean.wv6.mod/4+1))
utc_pred_resp_mean_df.wv6 <- as.data.frame(mapply(c,utc_pred_resp_mean.wv6, utc_pred_na_15as_sub))
utc_pred_resp_mean_order.wv6 <- orderBy(~y + x, data=utc_pred_resp_mean_df.wv6)
utc.pred.response.mean.grid.wv6 <- matrix(as.matrix(as.numeric(utc_pred_resp_mean_order.wv6$`log(utc.pred.response.mean.wv6.mod/4 + 1)`)), nrow=utc.x.res*2)

utc_pred_resp_sd.wv6 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.response.sd.wv6)
utc_pred_resp_sd_df.wv6 <- as.data.frame(mapply(c,utc_pred_resp_sd.wv6, utc_pred_na_15as_sub))
utc_pred_resp_sd_order.wv6 <- orderBy(~y + x, data=utc_pred_resp_sd_df.wv6)
utc.pred.response.sd.grid.wv6 <- matrix(as.matrix(as.numeric(utc_pred_resp_sd_order.wv6$utc.pred.response.sd.wv6)), nrow=utc.x.res*2)

utc_pred_lat_mean.wv6 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.mean.wv6)
utc_pred_lat_mean_df.wv6 <- as.data.frame(mapply(c,utc_pred_lat_mean.wv6, utc_pred_na_15as_sub))
utc_pred_lat_mean_order.wv6 <- orderBy(~y + x, data=utc_pred_lat_mean_df.wv6)
utc.pred.latent.mean.grid.wv6 <- matrix(as.matrix(as.numeric(utc_pred_lat_mean_order.wv6$utc.pred.latent.mean.wv6)), nrow=utc.x.res*2)

utc_pred_lat_sd.wv6 <- cbind(stacked_utc_pred_15as_fit[,c(1,2)],utc.pred.latent.sd.wv6)
utc_pred_lat_sd_df.wv6 <- as.data.frame(mapply(c,utc_pred_lat_sd.wv6, utc_pred_na_15as_sub))
utc_pred_lat_sd_order.wv6 <- orderBy(~y + x, data=utc_pred_lat_sd_df.wv6)
utc.pred.latent.sd.grid.wv6 <- matrix(as.matrix(as.numeric(utc_pred_lat_sd_order.wv6$utc.pred.latent.sd.wv6)), nrow=utc.x.res*2)

# on the response
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wv6),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0,12.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.response.mean.grid.wv6), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wv6-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,(utc.pred.latent.mean.grid.wv6-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the fixed effect
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wv6-utc.pred.latent.mean.grid.wv6)-log(4)),
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(-8,10)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as,((utc.pred.response.mean.grid.wv6-utc.pred.latent.mean.grid.wv6)-log(4)), add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# sd
# on the response  
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wv6,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(0, 10.1)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.response.sd.grid.wv6, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

# on the latent effects
image.plot(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wv6,
           xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T, legend.width=0.8, legend.mar=5,
           zlim = range(c(1, 8)))
contour(utc.seq.x.grid.15as,utc.seq.y.grid.15as, utc.pred.latent.sd.grid.wv6, add=T, lwd=2, labcex=1)
lines(utc_borderline,lwd=4)

utc.pred.response.mean.wv6.mod <- exp(utc.pred.response.mean.wv6)-1
utc.pred.response.mean.wv6.mod[utc.pred.response.mean.wv6.mod<0]<-0
res.15as.wv6<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.wv6.mod/4
RMSE.15as.wv6 <- sqrt(mean(res.15as.wv6^2))
RMSE.15as.wv6 
correl.15as.wv6 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.wv6.mod/4) 
correl.15as.wv6 

####################################
# LGCP worst case scenario + vague prior
# est at 3 arcsec sampled 100% data + pred at 3 arcsec
####################################
utc.stack.est.30as.lgcp.worst <- inla.stack(data = list(pop=round(as.numeric(stacked_utc_30as_fit$pop_utc_15_30as),digits=0)),
                                            A = list(utc.A.est.30as, 1),
                                            effects = list(c(utc.s.index, list(Intercept=1)), 
                                                           list(ntl=as.numeric(stacked_utc_30as_fit$ntl_utc_15_30as), 
                                                                slope=as.numeric(stacked_utc_30as_fit$slope_utc_30as)
                                                           )),
                                            tag="est")

utc.join.stack.pred.15as.lgcp.worst <- inla.stack(utc.stack.est.30as.lgcp.worst, utc.stack.latent, utc.stack.pred.15as)

utc.output.pred.15as.lgcp.worst <- inla(utc.formula, 
                                        data=inla.stack.data(utc.join.stack.pred.15as.lgcp.worst, spde=utc.spde), family="poisson",
                                        control.fixed=utc.fixed.vague, 
                                        control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.lgcp.worst), compute=TRUE, link=1),
                                        control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.lgcp.worst) 

utc.index.pred.response.lgcp.worst <- inla.stack.index(utc.join.stack.pred.15as.lgcp.worst, "pred.response")$data 

utc.pred.response.mean.lgcp.worst <- utc.output.pred.15as.lgcp.worst$summary.fitted.values[utc.index.pred.response.lgcp.worst, "mean"]

res.15as.lgcp.worst <- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.lgcp.worst/4
RMSE.15as.lgcp.worst <- sqrt(mean(res.15as.lgcp.worst^2))
RMSE.15as.lgcp.worst 
correl.15as.lgcp.worst <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.lgcp.worst/4) 
correl.15as.lgcp.worst 

####################################
# LGCP best case scenario + vague prior
# est at 3 arcsec sampled 100% data + pred at 3 arcsec
####################################
utc.stack.est.30as.lgcp.best <- inla.stack(data = list(pop=round(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as),digits=0)*4),
                                           A = list(utc.A.pred.15as, 1),
                                           effects = list(c(utc.s.index, list(Intercept=1)), 
                                                          list(ntl=as.numeric(stacked_utc_pred_15as_fit$ntl_utc_15_15as), 
                                                               slope=as.numeric(stacked_utc_pred_15as_fit$slope_utc_15as)
                                                          )),
                                           tag="est")

utc.join.stack.pred.15as.lgcp.best <- inla.stack(utc.stack.est.30as.lgcp.best, utc.stack.latent, utc.stack.pred.15as)

utc.output.pred.15as.lgcp.best <- inla(utc.formula, 
                                       data=inla.stack.data(utc.join.stack.pred.15as.lgcp.best, spde=utc.spde), family="poisson",
                                       control.fixed=utc.fixed.vague, 
                                       control.predictor=list(A=inla.stack.A(utc.join.stack.pred.15as.lgcp.best), compute=TRUE, link=1),
                                       control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(utc.output.pred.15as.lgcp.best) 

utc.index.pred.response.lgcp.best <- inla.stack.index(utc.join.stack.pred.15as.lgcp.best, "pred.response")$data 

utc.pred.response.mean.lgcp.best <- utc.output.pred.15as.lgcp.best$summary.fitted.values[utc.index.pred.response.lgcp.best, "mean"]

res.15as.lgcp.best <- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as) - utc.pred.response.mean.lgcp.best/4
RMSE.15as.lgcp.best <- sqrt(mean(res.15as.lgcp.best^2))
RMSE.15as.lgcp.best 
correl.15as.lgcp.best <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as), utc.pred.response.mean.lgcp.best/4) 
correl.15as.lgcp.best 

####################################
# bottom-up method best case scenario + vague prior
# est at 3 arcsec sampled 100% data + pred at 3 arcsec
####################################
library(rjags)
library(runjags)

bt.model <- "model{
  
  ## indexing
  # i = microcensus enumeration zone
  # t = residential type
  # k = regression coefficient
  
  for(i in 1:ni){
    
    # N = true number of people
    N[i] ~ dpois(D[i]*A[i])
    
    # D = true density of people
    D[i] ~ dlnorm(Dbar[i], pow(sigmaD0_region[type[i],region[i]],-2))
    
    # Dbar = expected density
    Dbar[i] <- alpha0_region[type[i],region[i]] + beta[1]*wpglobal[i] + beta[2]*schools[i] + beta[3]*hhsize[i] + beta[4]*settled[i]
    
    # posterior predictions
    Nhat[i] ~ dpois(Dhat[i]*A[i])
    Dhat[i] ~ dlnorm(Dbar[i], pow(sigmaD0_region[type[i],region[i]],-2))
  }
  
  # alpha = average population density (log-scale) for each settlement type per LGA (hierarchical by type/region/state/lga)
  mu_alpha_national ~ dnorm(0, 1e-3)
  sigma_alpha_national ~ dunif(0, 1e3)
  
  alpha0_national ~ dnorm(mu_alpha_national, pow(sigma_alpha_national, -2))
  
  for(t in 1:ntype){
    alpha0_type[t] ~ dnorm(mu_alpha_type[t], pow(sigma_alpha_type[t],-2))
    
    mu_alpha_type[t] ~ dnorm(mu_alpha_national, pow(sigma_alpha_national,-2))
    sigma_alpha_type[t] ~ dunif(0, sigma_alpha_national)
    
    for(r in 1:nregion){
      alpha0_region[t,r] ~ dnorm(mu_alpha_type[t], pow(sigma_alpha_type[t],-2))
    }
  }
  
  # beta = fixed effects of covariates
  for(b in 1:nbeta){
    beta[b] ~ dnorm(0, pow(5,-2))
  }
  
  # sigmaD = residual variance in population density (log-scale) for an LGA (hierarchical by type/region/state/lga)
  mu_sigmaD_national ~ dnorm(0, 1e-3) I(0,)
  sigma_sigmaD_national ~ dunif(0, 3)
  
  sigmaD0_national ~ dnorm(mu_sigmaD_national, pow(sigma_sigmaD_national, -2)) T(0,)
  
  for(t in 1:ntype){
    sigmaD0_type[t] ~ dnorm(mu_sigmaD_type[t], pow(sigma_sigmaD_type[t],-2)) T(0,)
    
    mu_sigmaD_type[t] ~ dnorm(mu_sigmaD_national, pow(sigma_sigmaD_national,-2)) T(0,)
    sigma_sigmaD_type[t] ~ dunif(0, sigma_sigmaD_national)
    
    for(r in 1:nregion){
      sigmaD0_region[t,r] ~ dnorm(mu_sigmaD_type[t], pow(sigma_sigmaD_type[t],-2)) T(0,)
    }
  }
}"

bt.data.full <- list(N=as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as_real), 
                     A=rep(1,dim(stacked_utc_pred_15as_fit)[1]),
                     type=as.factor(stacked_utc_pred_15as_fit$stt_utc_15_30as),
                     region=as.factor(stacked_utc_pred_15as_fit$muni_utc_15),
                     wpglobal=stacked_utc_pred_15as_fit$pop_utc_15_30as,
                     schools=stacked_utc_pred_15as_fit$psd_utc_15_15as,
                     hhsize=stacked_utc_pred_15as_fit$hhs_utc_15_15as,
                     settled=stacked_utc_pred_15as_fit$bsa_utc_15_15as/10000,
                     ni=dim(stacked_utc_pred_15as_fit)[1],
                     ntype=8,
                     nregion=26,
                     nbeta=4)

bt.out.full <- run.jags(bt.model, data=bt.data.full, n.chains=3, adapt = 1000, burnin = 1000, 
                        sample = 2000, monitor=c("mu_alpha_national", "sigma_alpha_national", 
                                                 "mu_sigmaD_national", "sigma_sigmaD_national", "beta"))
summary(bt.out.full)

utc.pred.response.mean.bt.full <- bt.out.full$summaries[1,4]+
  bt.out.full$summaries[5,4]*stacked_utc_pred_15as_fit$pop_utc_15_30as+
  bt.out.full$summaries[6,4]*stacked_utc_pred_15as_fit$psd_utc_15_15as+
  bt.out.full$summaries[7,4]*stacked_utc_pred_15as_fit$hhs_utc_15_15as+
  bt.out.full$summaries[8,4]*stacked_utc_pred_15as_fit$bsa_utc_15_15as/10000

res.15as.bt.full<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as_real) - exp(utc.pred.response.mean.bt.full)
RMSE.15as.bt.full <- sqrt(mean(res.15as.bt.full^2))
RMSE.15as.bt.full 
correl.15as.bt.full <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as_real), exp(utc.pred.response.mean.bt.full))
correl.15as.bt.full 

####################################
# bottom-up method + vague prior
# est at 3 arcsec sampled 1% data (not weighted) + pred at 3 arcsec
####################################
bt.data.u1 <- list(N=as.numeric(stacked_utc_sampled_u1$pop_utc_15_15as_real), 
                   A=rep(1,dim(stacked_utc_sampled_u1)[1]),
                   type=as.factor(stacked_utc_sampled_u1$stt_utc_15_30as),
                   region=as.factor(stacked_utc_sampled_u1$muni_utc_15),
                   wpglobal=stacked_utc_sampled_u1$pop_utc_15_30as,
                   schools=stacked_utc_sampled_u1$psd_utc_15_15as,
                   hhsize=stacked_utc_sampled_u1$hhs_utc_15_15as,
                   settled=stacked_utc_sampled_u1$bsa_utc_15_15as/10000,
                   ni=dim(stacked_utc_sampled_u1)[1],
                   ntype=8,
                   nregion=26,
                   nbeta=4)

bt.out.u1 <- run.jags(bt.model, data=bt.data.u1, n.chains=3, adapt = 1000, burnin = 1000, 
                      sample = 2000, monitor=c("mu_alpha_national", "sigma_alpha_national", 
                                               "mu_sigmaD_national", "sigma_sigmaD_national", "beta"))
summary(bt.out.u1)

par(mfrow=c(4,2), mar = c(5, 5, 2, 2))
traceplot(bt.out.u1$mcmc)

utc.pred.response.mean.bt.u1 <- bt.out.u1$summaries[1,4]+
  bt.out.u1$summaries[5,4]*stacked_utc_pred_15as_fit$pop_utc_15_30as+
  bt.out.u1$summaries[6,4]*stacked_utc_pred_15as_fit$psd_utc_15_15as+
  bt.out.u1$summaries[7,4]*stacked_utc_pred_15as_fit$hhs_utc_15_15as+
  bt.out.u1$summaries[8,4]*stacked_utc_pred_15as_fit$bsa_utc_15_15as/10000

res.15as.bt.u1<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as_real) - exp(utc.pred.response.mean.bt.u1)
RMSE.15as.bt.u1 <- sqrt(mean(res.15as.bt.u1^2))
RMSE.15as.bt.u1 
correl.15as.bt.u1 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as_real), exp(utc.pred.response.mean.bt.u1))
correl.15as.bt.u1 

####################################
# bottom-up method + vague prior
# est at 3 arcsec sampled 2% data (not weighted) + pred at 3 arcsec
####################################
bt.data.u2 <- list(N=as.numeric(stacked_utc_sampled_u2$pop_utc_15_15as_real), 
                   A=rep(1,dim(stacked_utc_sampled_u2)[1]),
                   type=as.factor(stacked_utc_sampled_u2$stt_utc_15_30as),
                   region=as.factor(stacked_utc_sampled_u2$muni_utc_15),
                   wpglobal=stacked_utc_sampled_u2$pop_utc_15_30as,
                   schools=stacked_utc_sampled_u2$psd_utc_15_15as,
                   hhsize=stacked_utc_sampled_u2$hhs_utc_15_15as,
                   settled=stacked_utc_sampled_u2$bsa_utc_15_15as/10000,
                   ni=dim(stacked_utc_sampled_u2)[1],
                   ntype=8,
                   nregion=26,
                   nbeta=4)

bt.out.u2 <- run.jags(bt.model, data=bt.data.u2, n.chains=3, adapt = 1000, burnin = 1000, 
                      sample = 2000, monitor=c("mu_alpha_national", "sigma_alpha_national", 
                                               "mu_sigmaD_national", "sigma_sigmaD_national", "beta"))
summary(bt.out.u2)

traceplot(bt.out.u2$mcmc)

utc.pred.response.mean.bt.u2 <- bt.out.u2$summaries[1,4]+
  bt.out.u2$summaries[5,4]*stacked_utc_pred_15as_fit$pop_utc_15_30as+
  bt.out.u2$summaries[6,4]*stacked_utc_pred_15as_fit$psd_utc_15_15as+
  bt.out.u2$summaries[7,4]*stacked_utc_pred_15as_fit$hhs_utc_15_15as+
  bt.out.u2$summaries[8,4]*stacked_utc_pred_15as_fit$bsa_utc_15_15as/10000

res.15as.bt.u2<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as_real) - exp(utc.pred.response.mean.bt.u2)
RMSE.15as.bt.u2 <- sqrt(mean(res.15as.bt.u2^2))
RMSE.15as.bt.u2 
correl.15as.bt.u2 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as_real), exp(utc.pred.response.mean.bt.u2))
correl.15as.bt.u2 

####################################
# bottom-up method + vague prior
# est at 3 arcsec sampled 5% data (not weighted) + pred at 3 arcsec
####################################
bt.data.u3 <- list(N=as.numeric(stacked_utc_sampled_u3$pop_utc_15_15as_real), 
                   A=rep(1,dim(stacked_utc_sampled_u3)[1]),
                   type=as.factor(stacked_utc_sampled_u3$stt_utc_15_30as),
                   region=as.factor(stacked_utc_sampled_u3$muni_utc_15),
                   wpglobal=stacked_utc_sampled_u3$pop_utc_15_30as,
                   schools=stacked_utc_sampled_u3$psd_utc_15_15as,
                   hhsize=stacked_utc_sampled_u3$hhs_utc_15_15as,
                   settled=stacked_utc_sampled_u3$bsa_utc_15_15as/10000,
                   ni=dim(stacked_utc_sampled_u3)[1],
                   ntype=8,
                   nregion=26,
                   nbeta=4)

bt.out.u3 <- run.jags(bt.model, data=bt.data.u3, n.chains=3, adapt = 1000, burnin = 1000, 
                      sample = 2000, monitor=c("mu_alpha_national", "sigma_alpha_national", 
                                               "mu_sigmaD_national", "sigma_sigmaD_national", "beta"))
summary(bt.out.u3)

traceplot(bt.out.u3$mcmc)

utc.pred.response.mean.bt.u3 <- bt.out.u3$summaries[1,4]+
  bt.out.u3$summaries[5,4]*stacked_utc_pred_15as_fit$pop_utc_15_30as+
  bt.out.u3$summaries[6,4]*stacked_utc_pred_15as_fit$psd_utc_15_15as+
  bt.out.u3$summaries[7,4]*stacked_utc_pred_15as_fit$hhs_utc_15_15as+
  bt.out.u3$summaries[8,4]*stacked_utc_pred_15as_fit$bsa_utc_15_15as/10000

res.15as.bt.u3<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as_real) - exp(utc.pred.response.mean.bt.u3)
RMSE.15as.bt.u3 <- sqrt(mean(res.15as.bt.u3^2))
RMSE.15as.bt.u3 
correl.15as.bt.u3 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as_real), exp(utc.pred.response.mean.bt.u3))
correl.15as.bt.u3 

####################################
# bottom-up method + vague prior
# est at 3 arcsec sampled 10% data (not weighted) + pred at 3 arcsec
####################################
bt.data.u4 <- list(N=as.numeric(stacked_utc_sampled_u4$pop_utc_15_15as_real), 
                   A=rep(1,dim(stacked_utc_sampled_u4)[1]),
                   type=as.factor(stacked_utc_sampled_u4$stt_utc_15_30as),
                   region=as.factor(stacked_utc_sampled_u4$muni_utc_15),
                   wpglobal=stacked_utc_sampled_u4$pop_utc_15_30as,
                   schools=stacked_utc_sampled_u4$psd_utc_15_15as,
                   hhsize=stacked_utc_sampled_u4$hhs_utc_15_15as,
                   settled=stacked_utc_sampled_u4$bsa_utc_15_15as/10000,
                   ni=dim(stacked_utc_sampled_u4)[1],
                   ntype=8,
                   nregion=26,
                   nbeta=4)

bt.out.u4 <- run.jags(bt.model, data=bt.data.u4, n.chains=3, adapt = 1000, burnin = 1000, 
                      sample = 2000, monitor=c("mu_alpha_national", "sigma_alpha_national", 
                                               "mu_sigmaD_national", "sigma_sigmaD_national", "beta"))
summary(bt.out.u4)

traceplot(bt.out.u4$mcmc)

utc.pred.response.mean.bt.u4 <- bt.out.u4$summaries[1,4]+
  bt.out.u4$summaries[5,4]*stacked_utc_pred_15as_fit$pop_utc_15_30as+
  bt.out.u4$summaries[6,4]*stacked_utc_pred_15as_fit$psd_utc_15_15as+
  bt.out.u4$summaries[7,4]*stacked_utc_pred_15as_fit$hhs_utc_15_15as+
  bt.out.u4$summaries[8,4]*stacked_utc_pred_15as_fit$bsa_utc_15_15as/10000

res.15as.bt.u4<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as_real) - exp(utc.pred.response.mean.bt.u4)
RMSE.15as.bt.u4 <- sqrt(mean(res.15as.bt.u4^2))
RMSE.15as.bt.u4 
correl.15as.bt.u4 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as_real), exp(utc.pred.response.mean.bt.u4))
correl.15as.bt.u4 

####################################
# bottom-up method + vague prior
# est at 3 arcsec sampled 20% data (not weighted) + pred at 3 arcsec
####################################
bt.data.u5 <- list(N=as.numeric(stacked_utc_sampled_u5$pop_utc_15_15as_real), 
                   A=rep(1,dim(stacked_utc_sampled_u5)[1]),
                   type=as.factor(stacked_utc_sampled_u5$stt_utc_15_30as),
                   region=as.factor(stacked_utc_sampled_u5$muni_utc_15),
                   wpglobal=stacked_utc_sampled_u5$pop_utc_15_30as,
                   schools=stacked_utc_sampled_u5$psd_utc_15_15as,
                   hhsize=stacked_utc_sampled_u5$hhs_utc_15_15as,
                   settled=stacked_utc_sampled_u5$bsa_utc_15_15as/10000,
                   ni=dim(stacked_utc_sampled_u5)[1],
                   ntype=8,
                   nregion=26,
                   nbeta=4)

bt.out.u5 <- run.jags(bt.model, data=bt.data.u5, n.chains=3, adapt = 1000, burnin = 1000, 
                      sample = 2000, monitor=c("mu_alpha_national", "sigma_alpha_national", 
                                               "mu_sigmaD_national", "sigma_sigmaD_national", "beta"))
summary(bt.out.u5)

traceplot(bt.out.u5$mcmc)

utc.pred.response.mean.bt.u5 <- bt.out.u5$summaries[1,4]+
  bt.out.u5$summaries[5,4]*stacked_utc_pred_15as_fit$pop_utc_15_30as+
  bt.out.u5$summaries[6,4]*stacked_utc_pred_15as_fit$psd_utc_15_15as+
  bt.out.u5$summaries[7,4]*stacked_utc_pred_15as_fit$hhs_utc_15_15as+
  bt.out.u5$summaries[8,4]*stacked_utc_pred_15as_fit$bsa_utc_15_15as/10000

res.15as.bt.u5<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as_real) - exp(utc.pred.response.mean.bt.u5)
RMSE.15as.bt.u5 <- sqrt(mean(res.15as.bt.u5^2))
RMSE.15as.bt.u5 
correl.15as.bt.u5 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as_real), exp(utc.pred.response.mean.bt.u5))
correl.15as.bt.u5 

####################################
# bottom-up method + vague prior
# est at 3 arcsec sampled 50% data (not weighted) + pred at 3 arcsec
####################################
bt.data.u6 <- list(N=as.numeric(stacked_utc_sampled_u6$pop_utc_15_15as_real), 
                   A=rep(1,dim(stacked_utc_sampled_u6)[1]),
                   type=as.factor(stacked_utc_sampled_u6$stt_utc_15_30as),
                   region=as.factor(stacked_utc_sampled_u6$muni_utc_15),
                   wpglobal=stacked_utc_sampled_u6$pop_utc_15_30as,
                   schools=stacked_utc_sampled_u6$psd_utc_15_15as,
                   hhsize=stacked_utc_sampled_u6$hhs_utc_15_15as,
                   settled=stacked_utc_sampled_u6$bsa_utc_15_15as/10000,
                   ni=dim(stacked_utc_sampled_u6)[1],
                   ntype=8,
                   nregion=26,
                   nbeta=4)

bt.out.u6 <- run.jags(bt.model, data=bt.data.u6, n.chains=3, adapt = 1000, burnin = 1000, 
                      sample = 2000, monitor=c("mu_alpha_national", "sigma_alpha_national", 
                                               "mu_sigmaD_national", "sigma_sigmaD_national", "beta"))
summary(bt.out.u6)

traceplot(bt.out.u6$mcmc)

utc.pred.response.mean.bt.u6 <- bt.out.u6$summaries[1,4]+
  bt.out.u6$summaries[5,4]*stacked_utc_pred_15as_fit$pop_utc_15_30as+
  bt.out.u6$summaries[6,4]*stacked_utc_pred_15as_fit$psd_utc_15_15as+
  bt.out.u6$summaries[7,4]*stacked_utc_pred_15as_fit$hhs_utc_15_15as+
  bt.out.u6$summaries[8,4]*stacked_utc_pred_15as_fit$bsa_utc_15_15as/10000

res.15as.bt.u6<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as_real) - exp(utc.pred.response.mean.bt.u6)
RMSE.15as.bt.u6 <- sqrt(mean(res.15as.bt.u6^2))
RMSE.15as.bt.u6 
correl.15as.bt.u6 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as_real), exp(utc.pred.response.mean.bt.u6))
correl.15as.bt.u6 

####################################
# bottom-up method + vague prior
# est at 3 arcsec sampled 1% data (weighted) + pred at 3 arcsec
####################################
bt.data.w1 <- list(N=as.numeric(stacked_utc_sampled_w1$pop_utc_15_15as_real), 
                   A=rep(1,dim(stacked_utc_sampled_w1)[1]),
                   type=as.factor(stacked_utc_sampled_w1$stt_utc_15_30as),
                   region=as.factor(stacked_utc_sampled_w1$muni_utc_15),
                   wpglobal=stacked_utc_sampled_w1$pop_utc_15_30as,
                   schools=stacked_utc_sampled_w1$psd_utc_15_15as,
                   hhsize=stacked_utc_sampled_w1$hhs_utc_15_15as,
                   settled=stacked_utc_sampled_w1$bsa_utc_15_15as/10000,
                   ni=dim(stacked_utc_sampled_w1)[1],
                   ntype=8,
                   nregion=26,
                   nbeta=4)

bt.out.w1 <- run.jags(bt.model, data=bt.data.w1, n.chains=3, adapt = 1000, burnin = 1000, 
                      sample = 2000, monitor=c("mu_alpha_national", "sigma_alpha_national", 
                                               "mu_sigmaD_national", "sigma_sigmaD_national", "beta"))
summary(bt.out.w1)

traceplot(bt.out.w1$mcmc)

utc.pred.response.mean.bt.w1 <- bt.out.w1$summaries[1,4]+
  bt.out.w1$summaries[5,4]*stacked_utc_pred_15as_fit$pop_utc_15_30as+
  bt.out.w1$summaries[6,4]*stacked_utc_pred_15as_fit$psd_utc_15_15as+
  bt.out.w1$summaries[7,4]*stacked_utc_pred_15as_fit$hhs_utc_15_15as+
  bt.out.w1$summaries[8,4]*stacked_utc_pred_15as_fit$bsa_utc_15_15as/10000

res.15as.bt.w1<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as_real) - exp(utc.pred.response.mean.bt.w1)
RMSE.15as.bt.w1 <- sqrt(mean(res.15as.bt.w1^2))
RMSE.15as.bt.w1 
correl.15as.bt.w1 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as_real), exp(utc.pred.response.mean.bt.w1))
correl.15as.bt.w1 

####################################
# bottom-up method + vague prior
# est at 3 arcsec sampled 2% data (weighted) + pred at 3 arcsec
####################################
bt.data.w2 <- list(N=as.numeric(stacked_utc_sampled_w2$pop_utc_15_15as_real), 
                   A=rep(1,dim(stacked_utc_sampled_w2)[1]),
                   type=as.factor(stacked_utc_sampled_w2$stt_utc_15_30as),
                   region=as.factor(stacked_utc_sampled_w2$muni_utc_15),
                   wpglobal=stacked_utc_sampled_w2$pop_utc_15_30as,
                   schools=stacked_utc_sampled_w2$psd_utc_15_15as,
                   hhsize=stacked_utc_sampled_w2$hhs_utc_15_15as,
                   settled=stacked_utc_sampled_w2$bsa_utc_15_15as/10000,
                   ni=dim(stacked_utc_sampled_w2)[1],
                   ntype=8,
                   nregion=26,
                   nbeta=4)

bt.out.w2 <- run.jags(bt.model, data=bt.data.w2, n.chains=3, adapt = 1000, burnin = 1000, 
                      sample = 2000, monitor=c("mu_alpha_national", "sigma_alpha_national", 
                                               "mu_sigmaD_national", "sigma_sigmaD_national", "beta"))
summary(bt.out.w2)

traceplot(bt.out.w2$mcmc)

utc.pred.response.mean.bt.w2 <- bt.out.w2$summaries[1,4]+
  bt.out.w2$summaries[5,4]*stacked_utc_pred_15as_fit$pop_utc_15_30as+
  bt.out.w2$summaries[6,4]*stacked_utc_pred_15as_fit$psd_utc_15_15as+
  bt.out.w2$summaries[7,4]*stacked_utc_pred_15as_fit$hhs_utc_15_15as+
  bt.out.w2$summaries[8,4]*stacked_utc_pred_15as_fit$bsa_utc_15_15as/10000

res.15as.bt.w2<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as_real) - exp(utc.pred.response.mean.bt.w2)
RMSE.15as.bt.w2 <- sqrt(mean(res.15as.bt.w2^2))
RMSE.15as.bt.w2 
correl.15as.bt.w2 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as_real), exp(utc.pred.response.mean.bt.w2))
correl.15as.bt.w2

####################################
# bottom-up method + vague prior
# est at 3 arcsec sampled 5% data (weighted) + pred at 3 arcsec
####################################
bt.data.w3 <- list(N=as.numeric(stacked_utc_sampled_w3$pop_utc_15_15as_real), 
                   A=rep(1,dim(stacked_utc_sampled_w3)[1]),
                   type=as.factor(stacked_utc_sampled_w3$stt_utc_15_30as),
                   region=as.factor(stacked_utc_sampled_w3$muni_utc_15),
                   wpglobal=stacked_utc_sampled_w3$pop_utc_15_30as,
                   schools=stacked_utc_sampled_w3$psd_utc_15_15as,
                   hhsize=stacked_utc_sampled_w3$hhs_utc_15_15as,
                   settled=stacked_utc_sampled_w3$bsa_utc_15_15as/10000,
                   ni=dim(stacked_utc_sampled_w3)[1],
                   ntype=8,
                   nregion=26,
                   nbeta=4)

bt.out.w3 <- run.jags(bt.model, data=bt.data.w3, n.chains=3, adapt = 1000, burnin = 1000, 
                      sample = 2000, monitor=c("mu_alpha_national", "sigma_alpha_national", 
                                               "mu_sigmaD_national", "sigma_sigmaD_national", "beta"))
summary(bt.out.w3)

traceplot(bt.out.w3$mcmc)

utc.pred.response.mean.bt.w3 <- bt.out.w3$summaries[1,4]+
  bt.out.w3$summaries[5,4]*stacked_utc_pred_15as_fit$pop_utc_15_30as+
  bt.out.w3$summaries[6,4]*stacked_utc_pred_15as_fit$psd_utc_15_15as+
  bt.out.w3$summaries[7,4]*stacked_utc_pred_15as_fit$hhs_utc_15_15as+
  bt.out.w3$summaries[8,4]*stacked_utc_pred_15as_fit$bsa_utc_15_15as/10000

res.15as.bt.w3<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as_real) - exp(utc.pred.response.mean.bt.w3)
RMSE.15as.bt.w3 <- sqrt(mean(res.15as.bt.w3^2))
RMSE.15as.bt.w3 
correl.15as.bt.w3 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as_real), exp(utc.pred.response.mean.bt.w3))
correl.15as.bt.w3

####################################
# bottom-up method + vague prior
# est at 3 arcsec sampled 10% data (weighted) + pred at 3 arcsec
####################################
bt.data.w4 <- list(N=as.numeric(stacked_utc_sampled_w4$pop_utc_15_15as_real), 
                   A=rep(1,dim(stacked_utc_sampled_w4)[1]),
                   type=as.factor(stacked_utc_sampled_w4$stt_utc_15_30as),
                   region=as.factor(stacked_utc_sampled_w4$muni_utc_15),
                   wpglobal=stacked_utc_sampled_w4$pop_utc_15_30as,
                   schools=stacked_utc_sampled_w4$psd_utc_15_15as,
                   hhsize=stacked_utc_sampled_w4$hhs_utc_15_15as,
                   settled=stacked_utc_sampled_w4$bsa_utc_15_15as/10000,
                   ni=dim(stacked_utc_sampled_w4)[1],
                   ntype=8,
                   nregion=26,
                   nbeta=4)

bt.out.w4 <- run.jags(bt.model, data=bt.data.w4, n.chains=3, adapt = 1000, burnin = 1000, 
                      sample = 2000, monitor=c("mu_alpha_national", "sigma_alpha_national", 
                                               "mu_sigmaD_national", "sigma_sigmaD_national", "beta"))
summary(bt.out.w4)

traceplot(bt.out.w4$mcmc)

utc.pred.response.mean.bt.w4 <- bt.out.w4$summaries[1,4]+
  bt.out.w4$summaries[5,4]*stacked_utc_pred_15as_fit$pop_utc_15_30as+
  bt.out.w4$summaries[6,4]*stacked_utc_pred_15as_fit$psd_utc_15_15as+
  bt.out.w4$summaries[7,4]*stacked_utc_pred_15as_fit$hhs_utc_15_15as+
  bt.out.w4$summaries[8,4]*stacked_utc_pred_15as_fit$bsa_utc_15_15as/10000

res.15as.bt.w4<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as_real) - exp(utc.pred.response.mean.bt.w4)
RMSE.15as.bt.w4 <- sqrt(mean(res.15as.bt.w4^2))
RMSE.15as.bt.w4 
correl.15as.bt.w4 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as_real), exp(utc.pred.response.mean.bt.w4))
correl.15as.bt.w4

####################################
# bottom-up method + vague prior
# est at 3 arcsec sampled 20% data (weighted) + pred at 3 arcsec
####################################
bt.data.w5 <- list(N=as.numeric(stacked_utc_sampled_w5$pop_utc_15_15as_real), 
                   A=rep(1,dim(stacked_utc_sampled_w5)[1]),
                   type=as.factor(stacked_utc_sampled_w5$stt_utc_15_30as),
                   region=as.factor(stacked_utc_sampled_w5$muni_utc_15),
                   wpglobal=stacked_utc_sampled_w5$pop_utc_15_30as,
                   schools=stacked_utc_sampled_w5$psd_utc_15_15as,
                   hhsize=stacked_utc_sampled_w5$hhs_utc_15_15as,
                   settled=stacked_utc_sampled_w5$bsa_utc_15_15as/10000,
                   ni=dim(stacked_utc_sampled_w5)[1],
                   ntype=8,
                   nregion=26,
                   nbeta=4)

bt.out.w5 <- run.jags(bt.model, data=bt.data.w5, n.chains=3, adapt = 1000, burnin = 1000, 
                      sample = 2000, monitor=c("mu_alpha_national", "sigma_alpha_national", 
                                               "mu_sigmaD_national", "sigma_sigmaD_national", "beta"))
summary(bt.out.w5)

traceplot(bt.out.w5$mcmc)

utc.pred.response.mean.bt.w5 <- bt.out.w5$summaries[1,4]+
  bt.out.w5$summaries[5,4]*stacked_utc_pred_15as_fit$pop_utc_15_30as+
  bt.out.w5$summaries[6,4]*stacked_utc_pred_15as_fit$psd_utc_15_15as+
  bt.out.w5$summaries[7,4]*stacked_utc_pred_15as_fit$hhs_utc_15_15as+
  bt.out.w5$summaries[8,4]*stacked_utc_pred_15as_fit$bsa_utc_15_15as/10000

res.15as.bt.w5<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as_real) - exp(utc.pred.response.mean.bt.w5)
RMSE.15as.bt.w5 <- sqrt(mean(res.15as.bt.w5^2))
RMSE.15as.bt.w5 
correl.15as.bt.w5 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as_real), exp(utc.pred.response.mean.bt.w5))
correl.15as.bt.w5

####################################
# bottom-up method + vague prior
# est at 3 arcsec sampled 50% data (weighted) + pred at 3 arcsec
####################################
bt.data.w6 <- list(N=as.numeric(stacked_utc_sampled_w6$pop_utc_15_15as_real), 
                   A=rep(1,dim(stacked_utc_sampled_w6)[1]),
                   type=as.factor(stacked_utc_sampled_w6$stt_utc_15_30as),
                   region=as.factor(stacked_utc_sampled_w6$muni_utc_15),
                   wpglobal=stacked_utc_sampled_w6$pop_utc_15_30as,
                   schools=stacked_utc_sampled_w6$psd_utc_15_15as,
                   hhsize=stacked_utc_sampled_w6$hhs_utc_15_15as,
                   settled=stacked_utc_sampled_w6$bsa_utc_15_15as/10000,
                   ni=dim(stacked_utc_sampled_w6)[1],
                   ntype=8,
                   nregion=26,
                   nbeta=4)

bt.out.w6 <- run.jags(bt.model, data=bt.data.w6, n.chains=3, adapt = 1000, burnin = 1000, 
                      sample = 2000, monitor=c("mu_alpha_national", "sigma_alpha_national", 
                                               "mu_sigmaD_national", "sigma_sigmaD_national", "beta"))
summary(bt.out.w6)

traceplot(bt.out.w6$mcmc)

utc.pred.response.mean.bt.w6 <- bt.out.w6$summaries[1,4]+
  bt.out.w6$summaries[5,4]*stacked_utc_pred_15as_fit$pop_utc_15_30as+
  bt.out.w6$summaries[6,4]*stacked_utc_pred_15as_fit$psd_utc_15_15as+
  bt.out.w6$summaries[7,4]*stacked_utc_pred_15as_fit$hhs_utc_15_15as+
  bt.out.w6$summaries[8,4]*stacked_utc_pred_15as_fit$bsa_utc_15_15as/10000

res.15as.bt.w6<- as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as_real) - exp(utc.pred.response.mean.bt.w6)
RMSE.15as.bt.w6 <- sqrt(mean(res.15as.bt.w6^2))
RMSE.15as.bt.w6 
correl.15as.bt.w6 <- cor(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as_real), exp(utc.pred.response.mean.bt.w6))
correl.15as.bt.w6

####################################
# check overdispersion and underdispersion 
####################################
par(mfrow=c(1,1), mar = c(5, 5, 2, 2))
hist.pit.lgcp.worst <- hist(utc.output.pred.15as.lgcp.worst$cpo$pit, breaks=20)
hist.pit.lgcp.best <- hist(utc.output.pred.15as.lgcp.best$cpo$pit, breaks=20)
hist.pit.log <- hist(utc.output.pred.15as.log$cpo$pit, breaks=20)
hist.pit.bsvp <- hist(utc.output.pred.15as.bsvp$cpo$pit, breaks=20)

# at 30 arcsec
plot(hist.pit.lgcp.worst, col=col_u, xlab="PIT", main=NULL)
plot(hist.pit.log, add=TRUE, col=col_w)
legend(x = "topright", legend=c("proposed","LGCP"), 
       col=c(col_w,col_u),pch=c(15,15),
       cex=0.8)

# at 15 arcsec
plot(hist.pit.lgcp.best, col=col_u, xlab="PIT", main=NULL)
plot(hist.pit.bsvp, add=TRUE, col=col_w)
legend(x = "topright", legend=c("proposed","LGCP"), 
       col=c(col_w,col_u),pch=c(15,15),
       cex=0.8)

####################################
# decide which transformation to be used
####################################
par(mfrow=c(3,2), mar = c(5, 5, 2, 2))
plot.a1.xmin <- min(c(utc.fixed1_marg[,1], utc.fixed1_marg.sqrt[,1], utc.fixed1_marg.log[,1]))
plot.a1.xmax <- max(c(utc.fixed1_marg[,1], utc.fixed1_marg.sqrt[,1], utc.fixed1_marg.log[,1]))
plot.a1.ymin <- min(c(utc.fixed1_marg[,2], utc.fixed1_marg.sqrt[,2], utc.fixed1_marg.log[,2]))
plot.a1.ymax <- max(c(utc.fixed1_marg[,2], utc.fixed1_marg.sqrt[,2], utc.fixed1_marg.log[,2]))
plot(utc.fixed1_marg,xlim=c(plot.a1.xmin,plot.a1.xmax),ylim=c(plot.a1.ymin,plot.a1.ymax),
     col="black",t="l",xlab=expression(beta[0]),ylab="")
lines(utc.fixed1_marg.sqrt,col="blue",t="l",xlab="",ylab="")
lines(utc.fixed1_marg.log,col="red",t="l",xlab="",ylab="")
legend(x = "topright", legend=c("no transformation","sqrt transformation","log transformation"), 
       col=c("black", "blue", "red"), lwd=c(1,1,1), cex=0.8)

plot.a2.xmin <- min(c(utc.fixed2_marg[,1], utc.fixed2_marg.sqrt[,1], utc.fixed2_marg.log[,1]))
plot.a2.xmax <- max(c(utc.fixed2_marg[,1], utc.fixed2_marg.sqrt[,1], utc.fixed2_marg.log[,1]))
plot.a2.ymin <- min(c(utc.fixed2_marg[,2], utc.fixed2_marg.sqrt[,2], utc.fixed2_marg.log[,2]))
plot.a2.ymax <- max(c(utc.fixed2_marg[,2], utc.fixed2_marg.sqrt[,2], utc.fixed2_marg.log[,2]))
plot(utc.fixed2_marg,xlim=c(plot.a2.xmin,plot.a2.xmax),ylim=c(plot.a2.ymin,plot.a2.ymax),
     col="black",t="l",xlab=expression(beta[1]),ylab="")
lines(utc.fixed2_marg.sqrt,col="blue",t="l",xlab="",ylab="")
lines(utc.fixed2_marg.log,col="red",t="l",xlab="",ylab="")
legend(x = "topright", legend=c("no transformation","sqrt transformation","log transformation"), 
       col=c("black", "blue", "red"), lwd=c(1,1,1), cex=0.8)

plot.a3.xmin <- min(c(utc.fixed3_marg[,1], utc.fixed3_marg.sqrt[,1], utc.fixed3_marg.log[,1]))
plot.a3.xmax <- max(c(utc.fixed3_marg[,1], utc.fixed3_marg.sqrt[,1], utc.fixed3_marg.log[,1]))
plot.a3.ymin <- min(c(utc.fixed3_marg[,2], utc.fixed3_marg.sqrt[,2], utc.fixed3_marg.log[,2]))
plot.a3.ymax <- max(c(utc.fixed3_marg[,2], utc.fixed3_marg.sqrt[,2], utc.fixed3_marg.log[,2]))
plot(utc.fixed3_marg,xlim=c(plot.a3.xmin,plot.a3.xmax),ylim=c(plot.a3.ymin,plot.a3.ymax),
     col="black",t="l",xlab=expression(beta[2]),ylab="")
lines(utc.fixed3_marg.sqrt,col="blue",t="l",xlab="",ylab="")
lines(utc.fixed3_marg.log,col="red",t="l",xlab="",ylab="")
legend(x = "topright", legend=c("no transformation","sqrt transformation","log transformation"), 
       col=c("black", "blue", "red"), lwd=c(1,1,1), cex=0.8)

plot.a4.xmin <- min(c(utc.sigma2e_marg[,1], utc.sigma2e_marg.sqrt[,1], utc.sigma2e_marg.log[,1]))
plot.a4.xmax <- max(c(utc.sigma2e_marg[,1], utc.sigma2e_marg.sqrt[,1], utc.sigma2e_marg.log[,1]))
plot.a4.ymin <- min(c(utc.sigma2e_marg[,2], utc.sigma2e_marg.sqrt[,2], utc.sigma2e_marg.log[,2]))
plot.a4.ymax <- max(c(utc.sigma2e_marg[,2], utc.sigma2e_marg.sqrt[,2], utc.sigma2e_marg.log[,2]))
plot(utc.sigma2e_marg,xlim=c(plot.a4.xmin,plot.a4.xmax),ylim=c(plot.a4.ymin,plot.a4.ymax),
     col="black",t="l",xlab=expression(sigma[e]^2),ylab="")
lines(utc.sigma2e_marg.sqrt,col="blue",t="l",xlab="",ylab="")
lines(utc.sigma2e_marg.log,col="red",t="l",xlab="",ylab="")
legend(x = "topright", legend=c("no transformation","sqrt transformation","log transformation"), 
       col=c("black", "blue", "red"), lwd=c(1,1,1), cex=0.8)

plot.a5.xmin <- min(c(utc.var.nom.marg[,1], utc.var.nom.marg.sqrt[,1], utc.var.nom.marg.log[,1]))
plot.a5.xmax <- max(c(utc.var.nom.marg[,1], utc.var.nom.marg.sqrt[,1], utc.var.nom.marg.log[,1]))
plot.a5.ymin <- min(c(utc.var.nom.marg[,2], utc.var.nom.marg.sqrt[,2], utc.var.nom.marg.log[,2]))
plot.a5.ymax <- max(c(utc.var.nom.marg[,2], utc.var.nom.marg.sqrt[,2], utc.var.nom.marg.log[,2]))
plot(utc.var.nom.marg,xlim=c(plot.a5.xmin,plot.a5.xmax),ylim=c(plot.a5.ymin,plot.a5.ymax),
     col="black",t="l",xlab=expression(sigma^2),ylab="")
lines(utc.var.nom.marg.sqrt,col="blue",t="l",xlab="",ylab="")
lines(utc.var.nom.marg.log,col="red",t="l",xlab="",ylab="")
legend(x = "topright", legend=c("no transformation","sqrt transformation","log transformation"), 
       col=c("black", "blue", "red"), lwd=c(1,1,1), cex=0.8)

plot.a6.xmin <- min(c(utc.range.nom.marg[,1], utc.range.nom.marg.sqrt[,1], utc.range.nom.marg.log[,1]))
plot.a6.xmax <- max(c(utc.range.nom.marg[,1], utc.range.nom.marg.sqrt[,1], utc.range.nom.marg.log[,1]))
plot.a6.ymin <- min(c(utc.range.nom.marg[,2], utc.range.nom.marg.sqrt[,2], utc.range.nom.marg.log[,2]))
plot.a6.ymax <- max(c(utc.range.nom.marg[,2], utc.range.nom.marg.sqrt[,2], utc.range.nom.marg.log[,2]))
plot(utc.range.nom.marg,xlim=c(plot.a6.xmin,plot.a6.xmax),ylim=c(plot.a6.ymin,plot.a6.ymax),
     col="black",t="l",xlab="r",ylab="")
lines(utc.range.nom.marg.sqrt,col="blue",t="l",xlab="",ylab="")
lines(utc.range.nom.marg.log,col="red",t="l",xlab="",ylab="")
legend(x = "topright", legend=c("no transformation","sqrt transformation","log transformation"), 
       col=c("black", "blue", "red"), lwd=c(1,1,1), cex=0.8)


####################################
# decide whether covariate play an important role
####################################
plot.b1.xmin <- min(c(utc.fixed1_marg.log[,1], utc.fixed1_marg.noboth[,1]))
plot.b1.xmax <- max(c(utc.fixed1_marg.log[,1], utc.fixed1_marg.noboth[,1]))
plot.b1.ymin <- min(c(utc.fixed1_marg.log[,2], utc.fixed1_marg.noboth[,2]))
plot.b1.ymax <- max(c(utc.fixed1_marg.log[,2], utc.fixed1_marg.noboth[,2]))
plot(utc.fixed1_marg.noboth,xlim=c(plot.b1.xmin,plot.b1.xmax),ylim=c(plot.b1.ymin,plot.b1.ymax),
     col="blue",t="l",xlab=expression(beta[0]),ylab="")
lines(utc.fixed1_marg.log,col="red",t="l",xlab="",ylab="")
legend(x = "topright", legend=c("w/o covariates","w/ covariate"), 
       col=c("blue", "red"), lwd=c(1,1), cex=0.8)

plot.b2.xmin <- min(c(utc.fixed2_marg.log[,1]))
plot.b2.xmax <- max(c(utc.fixed2_marg.log[,1]))
plot.b2.ymin <- min(c(utc.fixed2_marg.log[,2]))
plot.b2.ymax <- max(c(utc.fixed2_marg.log[,2]))
plot(utc.fixed2_marg.log,xlim=c(plot.b2.xmin,plot.b2.xmax),ylim=c(plot.b2.ymin,plot.b2.ymax),
     col="red",t="l",xlab=expression(beta[1]),ylab="")
legend(x = "topright", legend=c("w/ covariate"), 
       col=c("red"), lwd=c(1), cex=0.8)

plot.b3.xmin <- min(c(utc.fixed3_marg.log[,1]))
plot.b3.xmax <- max(c(utc.fixed3_marg.log[,1]))
plot.b3.ymin <- min(c(utc.fixed3_marg.log[,2]))
plot.b3.ymax <- max(c(utc.fixed3_marg.log[,2]))
plot(utc.fixed3_marg.log,xlim=c(plot.b3.xmin,plot.b3.xmax),ylim=c(plot.b3.ymin,plot.b3.ymax),
     col="red",t="l",xlab=expression(beta[2]),ylab="")
legend(x = "topright", legend=c("w/ covariate"), 
       col=c("red"), lwd=c(1), cex=0.8)

plot.b4.xmin <- min(c(utc.sigma2e_marg.noboth[,1], utc.sigma2e_marg.log[,1]))
plot.b4.xmax <- max(c(utc.sigma2e_marg.noboth[,1], utc.sigma2e_marg.log[,1]))
plot.b4.ymin <- min(c(utc.sigma2e_marg.noboth[,2], utc.sigma2e_marg.log[,2]))
plot.b4.ymax <- max(c(utc.sigma2e_marg.noboth[,2], utc.sigma2e_marg.log[,2]))
plot(utc.sigma2e_marg.noboth,xlim=c(plot.b4.xmin,plot.b4.xmax),ylim=c(plot.b4.ymin,plot.b4.ymax),
     col="blue",t="l",xlab=expression(sigma[e]^2),ylab="")
lines(utc.sigma2e_marg.log,col="red",t="l",xlab="",ylab="")
legend(x = "topright", legend=c("w/o covariates","w/ covariate"), 
       col=c("blue", "red"), lwd=c(1,1,1), cex=0.8)

plot.b5.xmin <- min(c(utc.var.nom.marg.noboth[,1], utc.var.nom.marg.log[,1]))
plot.b5.xmax <- max(c(utc.var.nom.marg.noboth[,1], utc.var.nom.marg.log[,1]))
plot.b5.ymin <- min(c(utc.var.nom.marg.noboth[,2], utc.var.nom.marg.log[,2]))
plot.b5.ymax <- max(c(utc.var.nom.marg.noboth[,2], utc.var.nom.marg.log[,2]))
plot(utc.var.nom.marg.noboth,xlim=c(plot.b5.xmin,plot.b5.xmax),ylim=c(plot.b5.ymin,plot.b5.ymax),
     col="blue",t="l",xlab=expression(sigma^2),ylab="")
lines(utc.var.nom.marg.log,col="red",t="l",xlab="",ylab="")
legend(x = "topright", legend=c("w/o covariates","w/ covariate"), 
       col=c("blue", "red"), lwd=c(1,1), cex=0.8)

plot.b6.xmin <- min(c(utc.range.nom.marg.noboth[,1], utc.range.nom.marg.log[,1]))
plot.b6.xmax <- max(c(utc.range.nom.marg.noboth[,1], utc.range.nom.marg.log[,1]))
plot.b6.ymin <- min(c(utc.range.nom.marg.noboth[,2], utc.range.nom.marg.log[,2]))
plot.b6.ymax <- max(c(utc.range.nom.marg.noboth[,2], utc.range.nom.marg.log[,2]))
plot(utc.range.nom.marg.noboth,xlim=c(plot.b6.xmin,plot.b6.xmax),ylim=c(plot.b6.ymin,plot.b6.ymax),
     col="blue",t="l",xlab="r",ylab="")
lines(utc.range.nom.marg.log,col="red",t="l",xlab="",ylab="")
legend(x = "topright", legend=c("w/o covariates","w/ covariate"), 
       col=c("blue", "red"), lwd=c(1,1), cex=0.8)


####################################
# explain why worst case scenario with best prior never work
####################################
plot.c1.xmin <- min(c(utc.fixed1_marg.bsvp[,1], utc.fixed1_marg.bswp.n[,1], utc.fixed1_marg.bswp.p[,1], 
                      utc.fixed1_marg.wsbp.n[,1], utc.fixed1_marg.wsbp.p[,1], utc.fixed1_marg.log[,1]))
plot.c1.xmax <- max(c(utc.fixed1_marg.bsvp[,1], utc.fixed1_marg.bswp.n[,1], utc.fixed1_marg.bswp.p[,1], 
                      utc.fixed1_marg.wsbp.n[,1], utc.fixed1_marg.wsbp.p[,1], utc.fixed1_marg.log[,1]))
plot.c1.ymin <- min(c(utc.fixed1_marg.bsvp[,2], utc.fixed1_marg.bswp.n[,2], utc.fixed1_marg.bswp.p[,2], 
                      utc.fixed1_marg.wsbp.n[,2], utc.fixed1_marg.wsbp.p[,2], utc.fixed1_marg.log[,2]))
plot.c1.ymax <- max(c(utc.fixed1_marg.bsvp[,2], utc.fixed1_marg.bswp.n[,2], utc.fixed1_marg.bswp.p[,2], 
                      utc.fixed1_marg.wsbp.n[,2], utc.fixed1_marg.wsbp.p[,2], utc.fixed1_marg.log[,2]))
plot(utc.fixed1_marg.bsvp,xlim=c(plot.c1.xmin,plot.c1.xmax),ylim=c(plot.c1.ymin,plot.c1.ymax),
     col="blue",t="l",lty=3,xlab=expression(beta[0]),ylab="")
lines(utc.fixed1_marg.bswp.n,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.bswp.p,col="blue",t="l",lty=2,xlab="",ylab="")
lines(utc.fixed1_marg.wsbp.n,col="red",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.wsbp.p,col="red",t="l",lty=2,xlab="",ylab="")
lines(utc.fixed1_marg.log,col="red",t="l",lty=3,xlab="",ylab="")
legend(x = "topright", legend=c("best data w/ vague prior","best data w/ worst n. prior",
                                "best data w/ worst p. prior","worst data w/ vague prior",
                                "worst data w/ best n. prior","worst data w/ best p. prior"), 
       col=c("blue", "blue", "blue", "red", "red", "red"), lwd=c(1,1,1,1,1,1), lty=c(3,1,2,3,1,2), cex=0.8)

plot.c2.xmin <- min(c(utc.fixed2_marg.bsvp[,1], utc.fixed2_marg.bswp.n[,1], utc.fixed2_marg.bswp.p[,1], 
                      utc.fixed2_marg.wsbp.n[,1], utc.fixed2_marg.wsbp.p[,1], utc.fixed2_marg.log[,1]))
plot.c2.xmax <- max(c(utc.fixed2_marg.bsvp[,1], utc.fixed2_marg.bswp.n[,1], utc.fixed2_marg.bswp.p[,1], 
                      utc.fixed2_marg.wsbp.n[,1], utc.fixed2_marg.wsbp.p[,1], utc.fixed2_marg.log[,1]))
plot.c2.ymin <- min(c(utc.fixed2_marg.bsvp[,2], utc.fixed2_marg.bswp.n[,2], utc.fixed2_marg.bswp.p[,2], 
                      utc.fixed2_marg.wsbp.n[,2], utc.fixed2_marg.wsbp.p[,2], utc.fixed2_marg.log[,2]))
plot.c2.ymax <- max(c(utc.fixed2_marg.bsvp[,2], utc.fixed2_marg.bswp.n[,2], utc.fixed2_marg.bswp.p[,2], 
                      utc.fixed2_marg.wsbp.n[,2], utc.fixed2_marg.wsbp.p[,2], utc.fixed2_marg.log[,2]))
plot(utc.fixed2_marg.bsvp,xlim=c(plot.c2.xmin,plot.c2.xmax),ylim=c(plot.c2.ymin,plot.c2.ymax),
     col="blue",t="l",lty=3,xlab=expression(beta[1]),ylab="")
lines(utc.fixed2_marg.bswp.n,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.bswp.p,col="blue",t="l",lty=2,xlab="",ylab="")
lines(utc.fixed2_marg.wsbp.n,col="red",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.wsbp.p,col="red",t="l",lty=2,xlab="",ylab="")
lines(utc.fixed2_marg.log,col="red",t="l",lty=3,xlab="",ylab="")
legend(x = "topright", legend=c("best data w/ vague prior","best data w/ worst n. prior",
                                "best data w/ worst p. prior","worst data w/ vague prior",
                                "worst data w/ best n. prior","worst data w/ best p. prior"), 
       col=c("blue", "blue", "blue", "red", "red", "red"), lwd=c(1,1,1,1,1,1), lty=c(3,1,2,3,1,2), cex=0.8)

plot.c3.xmin <- min(c(utc.fixed3_marg.bsvp[,1], utc.fixed3_marg.bswp.n[,1], utc.fixed3_marg.bswp.p[,1], 
                      utc.fixed3_marg.wsbp.n[,1], utc.fixed3_marg.wsbp.p[,1], utc.fixed3_marg.log[,1]))
plot.c3.xmax <- max(c(utc.fixed3_marg.bsvp[,1], utc.fixed3_marg.bswp.n[,1], utc.fixed3_marg.bswp.p[,1], 
                      utc.fixed3_marg.wsbp.n[,1], utc.fixed3_marg.wsbp.p[,1], utc.fixed3_marg.log[,1]))
plot.c3.ymin <- min(c(utc.fixed3_marg.bsvp[,2], utc.fixed3_marg.bswp.n[,2], utc.fixed3_marg.bswp.p[,2], 
                      utc.fixed3_marg.wsbp.n[,2], utc.fixed3_marg.wsbp.p[,2], utc.fixed3_marg.log[,2]))
plot.c3.ymax <- max(c(utc.fixed3_marg.bsvp[,2], utc.fixed3_marg.bswp.n[,2], utc.fixed3_marg.bswp.p[,2], 
                      utc.fixed3_marg.wsbp.n[,2], utc.fixed3_marg.wsbp.p[,2], utc.fixed3_marg.log[,2]))
plot(utc.fixed3_marg.bsvp,xlim=c(plot.c3.xmin,plot.c3.xmax),ylim=c(plot.c3.ymin,plot.c3.ymax),
     col="blue",t="l",lty=3,xlab=expression(beta[2]),ylab="")
lines(utc.fixed3_marg.bswp.n,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.bswp.p,col="blue",t="l",lty=2,xlab="",ylab="")
lines(utc.fixed3_marg.wsbp.n,col="red",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.wsbp.p,col="red",t="l",lty=2,xlab="",ylab="")
lines(utc.fixed3_marg.log,col="red",t="l",lty=3,xlab="",ylab="")
legend(x = "topright", legend=c("best data w/ vague prior","best data w/ worst n. prior",
                                "best data w/ worst p. prior","worst data w/ vague prior",
                                "worst data w/ best n. prior","worst data w/ best p. prior"), 
       col=c("blue", "blue", "blue", "red", "red", "red"), lwd=c(1,1,1,1,1,1), lty=c(3,1,2,3,1,2), cex=0.8)

plot.c4.xmin <- min(c(utc.sigma2e_marg.bsvp[,1], utc.sigma2e_marg.bswp.n[,1], utc.sigma2e_marg.bswp.p[,1], 
                      utc.sigma2e_marg.wsbp.n[,1], utc.sigma2e_marg.wsbp.p[,1], utc.sigma2e_marg.log[,1]))
plot.c4.xmax <- max(c(utc.sigma2e_marg.bsvp[,1], utc.sigma2e_marg.bswp.n[,1], utc.sigma2e_marg.bswp.p[,1], 
                      utc.sigma2e_marg.wsbp.n[,1], utc.sigma2e_marg.wsbp.p[,1], utc.sigma2e_marg.log[,1]))
plot.c4.ymin <- min(c(utc.sigma2e_marg.bsvp[,2], utc.sigma2e_marg.bswp.n[,2], utc.sigma2e_marg.bswp.p[,2], 
                      utc.sigma2e_marg.wsbp.n[,2], utc.sigma2e_marg.wsbp.p[,2], utc.sigma2e_marg.log[,2]))
plot.c4.ymax <- max(c(utc.sigma2e_marg.bsvp[,2], utc.sigma2e_marg.bswp.n[,2], utc.sigma2e_marg.bswp.p[,2], 
                      utc.sigma2e_marg.wsbp.n[,2], utc.sigma2e_marg.wsbp.p[,2], utc.sigma2e_marg.log[,2]))
plot(utc.sigma2e_marg.bsvp,xlim=c(plot.c4.xmin,plot.c4.xmax),ylim=c(plot.c4.ymin,plot.c4.ymax),
     col="blue",t="l",lty=3,xlab=expression(sigma[e]^2),ylab="")
lines(utc.sigma2e_marg.bswp.n,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.bswp.p,col="blue",t="l",lty=2,xlab="",ylab="")
lines(utc.sigma2e_marg.wsbp.n,col="red",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.wsbp.p,col="red",t="l",lty=2,xlab="",ylab="")
lines(utc.sigma2e_marg.log,col="red",t="l",lty=3,xlab="",ylab="")
legend(x = "topright", legend=c("best data w/ vague prior","best data w/ worst n. prior",
                                "best data w/ worst p. prior","worst data w/ vague prior",
                                "worst data w/ best n. prior","worst data w/ best p. prior"), 
       col=c("blue", "blue", "blue", "red", "red", "red"), lwd=c(1,1,1,1,1,1), lty=c(3,1,2,3,1,2), cex=0.8)

plot.c5.xmin <- min(c(utc.var.nom.marg.bsvp[,1], utc.var.nom.marg.bswp.n[,1], utc.var.nom.marg.bswp.p[,1], 
                      utc.var.nom.marg.wsbp.n[,1], utc.var.nom.marg.wsbp.p[,1], utc.var.nom.marg.log[,1]))
plot.c5.xmax <- max(c(utc.var.nom.marg.bsvp[,1], utc.var.nom.marg.bswp.n[,1], utc.var.nom.marg.bswp.p[,1], 
                      utc.var.nom.marg.wsbp.n[,1], utc.var.nom.marg.wsbp.p[,1], utc.var.nom.marg.log[,1]))
plot.c5.ymin <- min(c(utc.var.nom.marg.bsvp[,2], utc.var.nom.marg.bswp.n[,2], utc.var.nom.marg.bswp.p[,2], 
                      utc.var.nom.marg.wsbp.n[,2], utc.var.nom.marg.wsbp.p[,2], utc.var.nom.marg.log[,2]))
plot.c5.ymax <- max(c(utc.var.nom.marg.bsvp[,2], utc.var.nom.marg.bswp.n[,2], utc.var.nom.marg.bswp.p[,2], 
                      utc.var.nom.marg.wsbp.n[,2], utc.var.nom.marg.wsbp.p[,2], utc.var.nom.marg.log[,2]))
plot(utc.var.nom.marg.bsvp,xlim=c(plot.c5.xmin,plot.c5.xmax),ylim=c(plot.c5.ymin,plot.c5.ymax),
     col="blue",t="l",lty=3,xlab=expression(sigma^2),ylab="")
lines(utc.var.nom.marg.bswp.n,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.bswp.p,col="blue",t="l",lty=2,xlab="",ylab="")
lines(utc.var.nom.marg.wsbp.n,col="red",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.wsbp.p,col="red",t="l",lty=2,xlab="",ylab="")
lines(utc.var.nom.marg.log,col="red",t="l",lty=3,xlab="",ylab="")
legend(x = "topright", legend=c("best data w/ vague prior","best data w/ worst n. prior",
                                "best data w/ worst p. prior","worst data w/ vague prior",
                                "worst data w/ best n. prior","worst data w/ best p. prior"), 
       col=c("blue", "blue", "blue", "red", "red", "red"), lwd=c(1,1,1,1,1,1), lty=c(3,1,2,3,1,2), cex=0.8)

plot.c6.xmin <- min(c(utc.range.nom.marg.bsvp[,1], utc.range.nom.marg.bswp.n[,1], utc.range.nom.marg.bswp.p[,1], 
                      utc.range.nom.marg.wsbp.n[,1], utc.range.nom.marg.wsbp.p[,1], utc.range.nom.marg.log[,1]))
plot.c6.xmax <- max(c(utc.range.nom.marg.bsvp[,1], utc.range.nom.marg.bswp.n[,1], utc.range.nom.marg.bswp.p[,1], 
                      utc.range.nom.marg.wsbp.n[,1], utc.range.nom.marg.wsbp.p[,1], utc.range.nom.marg.log[,1]))
plot.c6.ymin <- min(c(utc.range.nom.marg.bsvp[,2], utc.range.nom.marg.bswp.n[,2], utc.range.nom.marg.bswp.p[,2], 
                      utc.range.nom.marg.wsbp.n[,2], utc.range.nom.marg.wsbp.p[,2], utc.range.nom.marg.log[,2]))
plot.c6.ymax <- max(c(utc.range.nom.marg.bsvp[,2], utc.range.nom.marg.bswp.n[,2], utc.range.nom.marg.bswp.p[,2], 
                      utc.range.nom.marg.wsbp.n[,2], utc.range.nom.marg.wsbp.p[,2], utc.range.nom.marg.log[,2]))
plot(utc.range.nom.marg.bsvp,xlim=c(plot.c6.xmin,plot.c6.xmax),ylim=c(plot.c6.ymin,plot.c6.ymax),
     col="blue",t="l",lty=3,xlab="r",ylab="")
lines(utc.range.nom.marg.bswp.n,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.bswp.p,col="blue",t="l",lty=2,xlab="",ylab="")
lines(utc.range.nom.marg.wsbp.n,col="red",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.wsbp.p,col="red",t="l",lty=2,xlab="",ylab="")
lines(utc.range.nom.marg.log,col="red",t="l",lty=3,xlab="",ylab="")
legend(x = "topright", legend=c("best data w/ vague prior","best data w/ worst n. prior",
                                "best data w/ worst p. prior","worst data w/ vague prior",
                                "worst data w/ best n. prior","worst data w/ best p. prior"), 
       col=c("blue", "blue", "blue", "red", "red", "red"), lwd=c(1,1,1,1,1,1), lty=c(3,1,2,3,1,2), cex=0.8)

####################################
# show results with unweighted sampling and normal prior
####################################
plot.d1.xmin <- min(c(utc.fixed1_marg.un1[,1], utc.fixed1_marg.un2[,1], utc.fixed1_marg.un3[,1], 
                      utc.fixed1_marg.un4[,1], utc.fixed1_marg.un5[,1], utc.fixed1_marg.un6[,1], 
                      utc.fixed1_marg.bswp.n[,1]))
plot.d1.xmax <- max(c(utc.fixed1_marg.un1[,1], utc.fixed1_marg.un2[,1], utc.fixed1_marg.un3[,1], 
                      utc.fixed1_marg.un4[,1], utc.fixed1_marg.un5[,1], utc.fixed1_marg.un6[,1], 
                      utc.fixed1_marg.bswp.n[,1]))
plot.d1.ymin <- min(c(utc.fixed1_marg.un1[,2], utc.fixed1_marg.un2[,2], utc.fixed1_marg.un3[,2], 
                      utc.fixed1_marg.un4[,2], utc.fixed1_marg.un5[,2], utc.fixed1_marg.un6[,2], 
                      utc.fixed1_marg.bswp.n[,2]))
plot.d1.ymax <- max(c(utc.fixed1_marg.un1[,2], utc.fixed1_marg.un2[,2], utc.fixed1_marg.un3[,2], 
                      utc.fixed1_marg.un4[,2], utc.fixed1_marg.un5[,2], utc.fixed1_marg.un6[,2], 
                      utc.fixed1_marg.bswp.n[,2]))
plot(utc.fixed1_marg.un1,xlim=c(plot.d1.xmin,plot.d1.xmax),ylim=c(plot.d1.ymin,plot.d1.ymax),
     col="red",t="l",lty=1,xlab=expression(beta[0]),ylab="")
lines(utc.fixed1_marg.un2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.un3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.un4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.un5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.un6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.bswp.n,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ worst n. prior","2% data w/ worst n. prior",
                                "5% data w/ worst n. prior","10% data w/ worst n. prior",
                                "20% data w/ worst n. prior", "50% data w/ worst n. prior",
                                "100% data w/ worst n. prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

plot.d2.xmin <- min(c(utc.fixed2_marg.un1[,1], utc.fixed2_marg.un2[,1], utc.fixed2_marg.un3[,1], 
                      utc.fixed2_marg.un4[,1], utc.fixed2_marg.un5[,1], utc.fixed2_marg.un6[,1], 
                      utc.fixed2_marg.bswp.n[,1]))
plot.d2.xmax <- max(c(utc.fixed2_marg.un1[,1], utc.fixed2_marg.un2[,1], utc.fixed2_marg.un3[,1], 
                      utc.fixed2_marg.un4[,1], utc.fixed2_marg.un5[,1], utc.fixed2_marg.un6[,1], 
                      utc.fixed2_marg.bswp.n[,1]))
plot.d2.ymin <- min(c(utc.fixed2_marg.un1[,2], utc.fixed2_marg.un2[,2], utc.fixed2_marg.un3[,2], 
                      utc.fixed2_marg.un4[,2], utc.fixed2_marg.un5[,2], utc.fixed2_marg.un6[,2], 
                      utc.fixed2_marg.bswp.n[,2]))
plot.d2.ymax <- max(c(utc.fixed2_marg.un1[,2], utc.fixed2_marg.un2[,2], utc.fixed2_marg.un3[,2], 
                      utc.fixed2_marg.un4[,2], utc.fixed2_marg.un5[,2], utc.fixed2_marg.un6[,2], 
                      utc.fixed2_marg.bswp.n[,2]))
plot(utc.fixed2_marg.un1,xlim=c(plot.d2.xmin,plot.d2.xmax),ylim=c(plot.d2.ymin,plot.d2.ymax),
     col="red",t="l",lty=1,xlab=expression(beta[1]),ylab="")
lines(utc.fixed2_marg.un2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.un3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.un4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.un5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.un6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.bswp.n,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ worst n. prior","2% data w/ worst n. prior",
                                "5% data w/ worst n. prior","10% data w/ worst n. prior",
                                "20% data w/ worst n. prior", "50% data w/ worst n. prior",
                                "100% data w/ worst n. prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

plot.d3.xmin <- min(c(utc.fixed3_marg.un1[,1], utc.fixed3_marg.un2[,1], utc.fixed3_marg.un3[,1], 
                      utc.fixed3_marg.un4[,1], utc.fixed3_marg.un5[,1], utc.fixed3_marg.un6[,1], 
                      utc.fixed3_marg.bswp.n[,1]))
plot.d3.xmax <- max(c(utc.fixed3_marg.un1[,1], utc.fixed3_marg.un2[,1], utc.fixed3_marg.un3[,1], 
                      utc.fixed3_marg.un4[,1], utc.fixed3_marg.un5[,1], utc.fixed3_marg.un6[,1], 
                      utc.fixed3_marg.bswp.n[,1]))
plot.d3.ymin <- min(c(utc.fixed3_marg.un1[,2], utc.fixed3_marg.un2[,2], utc.fixed3_marg.un3[,2], 
                      utc.fixed3_marg.un4[,2], utc.fixed3_marg.un5[,2], utc.fixed3_marg.un6[,2], 
                      utc.fixed3_marg.bswp.n[,2]))
plot.d3.ymax <- max(c(utc.fixed3_marg.un1[,2], utc.fixed3_marg.un2[,2], utc.fixed3_marg.un3[,2], 
                      utc.fixed3_marg.un4[,2], utc.fixed3_marg.un5[,2], utc.fixed3_marg.un6[,2], 
                      utc.fixed3_marg.bswp.n[,2]))
plot(utc.fixed3_marg.un1,xlim=c(plot.d3.xmin,plot.d3.xmax),ylim=c(plot.d3.ymin,plot.d3.ymax),
     col="red",t="l",lty=1,xlab=expression(beta[2]),ylab="")
lines(utc.fixed3_marg.un2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.un3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.un4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.un5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.un6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.bswp.n,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ worst n. prior","2% data w/ worst n. prior",
                                "5% data w/ worst n. prior","10% data w/ worst n. prior",
                                "20% data w/ worst n. prior", "50% data w/ worst n. prior",
                                "100% data w/ worst n. prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

plot.d4.xmin <- min(c(utc.sigma2e_marg.un1[,1], utc.sigma2e_marg.un2[,1], utc.sigma2e_marg.un3[,1], 
                      utc.sigma2e_marg.un4[,1], utc.sigma2e_marg.un5[,1], utc.sigma2e_marg.un6[,1], 
                      utc.sigma2e_marg.bswp.n[,1]))
plot.d4.xmax <- max(c(utc.sigma2e_marg.un1[,1], utc.sigma2e_marg.un2[,1], utc.sigma2e_marg.un3[,1], 
                      utc.sigma2e_marg.un4[,1], utc.sigma2e_marg.un5[,1], utc.sigma2e_marg.un6[,1], 
                      utc.sigma2e_marg.bswp.n[,1]))
plot.d4.ymin <- min(c(utc.sigma2e_marg.un1[,2], utc.sigma2e_marg.un2[,2], utc.sigma2e_marg.un3[,2], 
                      utc.sigma2e_marg.un4[,2], utc.sigma2e_marg.un5[,2], utc.sigma2e_marg.un6[,2], 
                      utc.sigma2e_marg.bswp.n[,2]))
plot.d4.ymax <- max(c(utc.sigma2e_marg.un1[,2], utc.sigma2e_marg.un2[,2], utc.sigma2e_marg.un3[,2], 
                      utc.sigma2e_marg.un4[,2], utc.sigma2e_marg.un5[,2], utc.sigma2e_marg.un6[,2], 
                      utc.sigma2e_marg.bswp.n[,2]))
plot(utc.sigma2e_marg.un1,xlim=c(plot.d4.xmin,plot.d4.xmax),ylim=c(plot.d4.ymin,plot.d4.ymax),
     col="red",t="l",lty=1,xlab=expression(sigma[e]^2),ylab="")
lines(utc.sigma2e_marg.un2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.un3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.un4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.un5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.un6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.bswp.n,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ worst n. prior","2% data w/ worst n. prior",
                                "5% data w/ worst n. prior","10% data w/ worst n. prior",
                                "20% data w/ worst n. prior", "50% data w/ worst n. prior",
                                "100% data w/ worst n. prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

plot.d5.xmin <- min(c(utc.var.nom.marg.un1[,1], utc.var.nom.marg.un2[,1], utc.var.nom.marg.un3[,1], 
                      utc.var.nom.marg.un4[,1], utc.var.nom.marg.un5[,1], utc.var.nom.marg.un6[,1], 
                      utc.var.nom.marg.bswp.n[,1]))
plot.d5.xmax <- max(c(utc.var.nom.marg.un1[,1], utc.var.nom.marg.un2[,1], utc.var.nom.marg.un3[,1], 
                      utc.var.nom.marg.un4[,1], utc.var.nom.marg.un5[,1], utc.var.nom.marg.un6[,1], 
                      utc.var.nom.marg.bswp.n[,1]))
plot.d5.ymin <- min(c(utc.var.nom.marg.un1[,2], utc.var.nom.marg.un2[,2], utc.var.nom.marg.un3[,2], 
                      utc.var.nom.marg.un4[,2], utc.var.nom.marg.un5[,2], utc.var.nom.marg.un6[,2], 
                      utc.var.nom.marg.bswp.n[,2]))
plot.d5.ymax <- max(c(utc.var.nom.marg.un1[,2], utc.var.nom.marg.un2[,2], utc.var.nom.marg.un3[,2], 
                      utc.var.nom.marg.un4[,2], utc.var.nom.marg.un5[,2], utc.var.nom.marg.un6[,2], 
                      utc.var.nom.marg.bswp.n[,2]))
plot(utc.var.nom.marg.un1,xlim=c(plot.d5.xmin,plot.d5.xmax),ylim=c(plot.d5.ymin,plot.d5.ymax),
     col="red",t="l",lty=1,xlab=expression(sigma^2),ylab="")
lines(utc.var.nom.marg.un2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.un3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.un4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.un5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.un6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.bswp.n,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ worst n. prior","2% data w/ worst n. prior",
                                "5% data w/ worst n. prior","10% data w/ worst n. prior",
                                "20% data w/ worst n. prior", "50% data w/ worst n. prior",
                                "100% data w/ worst n. prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

plot.d6.xmin <- min(c(utc.range.nom.marg.un1[,1], utc.range.nom.marg.un2[,1], utc.range.nom.marg.un3[,1], 
                      utc.range.nom.marg.un4[,1], utc.range.nom.marg.un5[,1], utc.range.nom.marg.un6[,1], 
                      utc.range.nom.marg.bswp.n[,1]))
plot.d6.xmax <- max(c(utc.range.nom.marg.un1[,1], utc.range.nom.marg.un2[,1], utc.range.nom.marg.un3[,1], 
                      utc.range.nom.marg.un4[,1], utc.range.nom.marg.un5[,1], utc.range.nom.marg.un6[,1], 
                      utc.range.nom.marg.bswp.n[,1]))
plot.d6.ymin <- min(c(utc.range.nom.marg.un1[,2], utc.range.nom.marg.un2[,2], utc.range.nom.marg.un3[,2], 
                      utc.range.nom.marg.un4[,2], utc.range.nom.marg.un5[,2], utc.range.nom.marg.un6[,2], 
                      utc.range.nom.marg.bswp.n[,2]))
plot.d6.ymax <- max(c(utc.range.nom.marg.un1[,2], utc.range.nom.marg.un2[,2], utc.range.nom.marg.un3[,2], 
                      utc.range.nom.marg.un4[,2], utc.range.nom.marg.un5[,2], utc.range.nom.marg.un6[,2], 
                      utc.range.nom.marg.bswp.n[,2]))
plot(utc.range.nom.marg.un1,xlim=c(plot.d6.xmin,plot.d6.xmax),ylim=c(plot.d6.ymin,plot.d6.ymax),
     col="red",t="l",lty=1,xlab="r",ylab="")
lines(utc.range.nom.marg.un2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.un3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.un4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.un5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.un6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.bswp.n,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ worst n. prior","2% data w/ worst n. prior",
                                "5% data w/ worst n. prior","10% data w/ worst n. prior",
                                "20% data w/ worst n. prior", "50% data w/ worst n. prior",
                                "100% data w/ worst n. prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

####################################
# show results with unweighted sampling and pc prior
####################################
plot.e1.xmin <- min(c(utc.fixed1_marg.up1[,1], utc.fixed1_marg.up2[,1], utc.fixed1_marg.up3[,1], 
                      utc.fixed1_marg.up4[,1], utc.fixed1_marg.up5[,1], utc.fixed1_marg.up6[,1], 
                      utc.fixed1_marg.bswp.p[,1]))
plot.e1.xmax <- max(c(utc.fixed1_marg.up1[,1], utc.fixed1_marg.up2[,1], utc.fixed1_marg.up3[,1], 
                      utc.fixed1_marg.up4[,1], utc.fixed1_marg.up5[,1], utc.fixed1_marg.up6[,1], 
                      utc.fixed1_marg.bswp.p[,1]))
plot.e1.ymin <- min(c(utc.fixed1_marg.up1[,2], utc.fixed1_marg.up2[,2], utc.fixed1_marg.up3[,2], 
                      utc.fixed1_marg.up4[,2], utc.fixed1_marg.up5[,2], utc.fixed1_marg.up6[,2], 
                      utc.fixed1_marg.bswp.p[,2]))
plot.e1.ymax <- max(c(utc.fixed1_marg.up1[,2], utc.fixed1_marg.up2[,2], utc.fixed1_marg.up3[,2], 
                      utc.fixed1_marg.up4[,2], utc.fixed1_marg.up5[,2], utc.fixed1_marg.up6[,2], 
                      utc.fixed1_marg.bswp.p[,2]))
plot(utc.fixed1_marg.up1,xlim=c(plot.e1.xmin,plot.e1.xmax),ylim=c(plot.e1.ymin,plot.e1.ymax),
     col="red",t="l",lty=1,xlab=expression(beta[0]),ylab="")
lines(utc.fixed1_marg.up2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.up3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.up4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.up5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.up6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.bswp.p,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ worst p. prior","2% data w/ worst p. prior",
                                "5% data w/ worst p. prior","10% data w/ worst p. prior",
                                "20% data w/ worst p. prior", "50% data w/ worst p. prior",
                                "100% data w/ worst p. prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

plot.e2.xmin <- min(c(utc.fixed2_marg.up1[,1], utc.fixed2_marg.up2[,1], utc.fixed2_marg.up3[,1], 
                      utc.fixed2_marg.up4[,1], utc.fixed2_marg.up5[,1], utc.fixed2_marg.up6[,1], 
                      utc.fixed2_marg.bswp.p[,1]))
plot.e2.xmax <- max(c(utc.fixed2_marg.up1[,1], utc.fixed2_marg.up2[,1], utc.fixed2_marg.up3[,1], 
                      utc.fixed2_marg.up4[,1], utc.fixed2_marg.up5[,1], utc.fixed2_marg.up6[,1], 
                      utc.fixed2_marg.bswp.p[,1]))
plot.e2.ymin <- min(c(utc.fixed2_marg.up1[,2], utc.fixed2_marg.up2[,2], utc.fixed2_marg.up3[,2], 
                      utc.fixed2_marg.up4[,2], utc.fixed2_marg.up5[,2], utc.fixed2_marg.up6[,2], 
                      utc.fixed2_marg.bswp.p[,2]))
plot.e2.ymax <- max(c(utc.fixed2_marg.up1[,2], utc.fixed2_marg.up2[,2], utc.fixed2_marg.up3[,2], 
                      utc.fixed2_marg.up4[,2], utc.fixed2_marg.up5[,2], utc.fixed2_marg.up6[,2], 
                      utc.fixed2_marg.bswp.p[,2]))
plot(utc.fixed2_marg.up1,xlim=c(plot.e2.xmin,plot.e2.xmax),ylim=c(plot.e2.ymin,plot.e2.ymax),
     col="red",t="l",lty=1,xlab=expression(beta[1]),ylab="")
lines(utc.fixed2_marg.up2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.up3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.up4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.up5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.up6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.bswp.p,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ worst p. prior","2% data w/ worst p. prior",
                                "5% data w/ worst p. prior","10% data w/ worst p. prior",
                                "20% data w/ worst p. prior", "50% data w/ worst p. prior",
                                "100% data w/ worst p. prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

plot.e3.xmin <- min(c(utc.fixed3_marg.up1[,1], utc.fixed3_marg.up2[,1], utc.fixed3_marg.up3[,1], 
                      utc.fixed3_marg.up4[,1], utc.fixed3_marg.up5[,1], utc.fixed3_marg.up6[,1], 
                      utc.fixed3_marg.bswp.p[,1]))
plot.e3.xmax <- max(c(utc.fixed3_marg.up1[,1], utc.fixed3_marg.up2[,1], utc.fixed3_marg.up3[,1], 
                      utc.fixed3_marg.up4[,1], utc.fixed3_marg.up5[,1], utc.fixed3_marg.up6[,1], 
                      utc.fixed3_marg.bswp.p[,1]))
plot.e3.ymin <- min(c(utc.fixed3_marg.up1[,2], utc.fixed3_marg.up2[,2], utc.fixed3_marg.up3[,2], 
                      utc.fixed3_marg.up4[,2], utc.fixed3_marg.up5[,2], utc.fixed3_marg.up6[,2], 
                      utc.fixed3_marg.bswp.p[,2]))
plot.e3.ymax <- max(c(utc.fixed3_marg.up1[,2], utc.fixed3_marg.up2[,2], utc.fixed3_marg.up3[,2], 
                      utc.fixed3_marg.up4[,2], utc.fixed3_marg.up5[,2], utc.fixed3_marg.up6[,2], 
                      utc.fixed3_marg.bswp.p[,2]))
plot(utc.fixed3_marg.up1,xlim=c(plot.e3.xmin,plot.e3.xmax),ylim=c(plot.e3.ymin,plot.e3.ymax),
     col="red",t="l",lty=1,xlab=expression(beta[2]),ylab="")
lines(utc.fixed3_marg.up2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.up3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.up4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.up5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.up6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.bswp.p,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ worst p. prior","2% data w/ worst p. prior",
                                "5% data w/ worst p. prior","10% data w/ worst p. prior",
                                "20% data w/ worst p. prior", "50% data w/ worst p. prior",
                                "100% data w/ worst p. prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

plot.e4.xmin <- min(c(utc.sigma2e_marg.up1[,1], utc.sigma2e_marg.up2[,1], utc.sigma2e_marg.up3[,1], 
                      utc.sigma2e_marg.up4[,1], utc.sigma2e_marg.up5[,1], utc.sigma2e_marg.up6[,1], 
                      utc.sigma2e_marg.bswp.p[,1]))
plot.e4.xmax <- max(c(utc.sigma2e_marg.up1[,1], utc.sigma2e_marg.up2[,1], utc.sigma2e_marg.up3[,1], 
                      utc.sigma2e_marg.up4[,1], utc.sigma2e_marg.up5[,1], utc.sigma2e_marg.up6[,1], 
                      utc.sigma2e_marg.bswp.p[,1]))
plot.e4.ymin <- min(c(utc.sigma2e_marg.up1[,2], utc.sigma2e_marg.up2[,2], utc.sigma2e_marg.up3[,2], 
                      utc.sigma2e_marg.up4[,2], utc.sigma2e_marg.up5[,2], utc.sigma2e_marg.up6[,2], 
                      utc.sigma2e_marg.bswp.p[,2]))
plot.e4.ymax <- max(c(utc.sigma2e_marg.up1[,2], utc.sigma2e_marg.up2[,2], utc.sigma2e_marg.up3[,2], 
                      utc.sigma2e_marg.up4[,2], utc.sigma2e_marg.up5[,2], utc.sigma2e_marg.up6[,2], 
                      utc.sigma2e_marg.bswp.p[,2]))
plot(utc.sigma2e_marg.up1,xlim=c(plot.e4.xmin,plot.e4.xmax),ylim=c(plot.e4.ymin,plot.e4.ymax),
     col="red",t="l",lty=1,xlab=expression(sigma[e]^2),ylab="")
lines(utc.sigma2e_marg.up2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.up3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.up4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.up5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.up6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.bswp.p,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ worst p. prior","2% data w/ worst p. prior",
                                "5% data w/ worst p. prior","10% data w/ worst p. prior",
                                "20% data w/ worst p. prior", "50% data w/ worst p. prior",
                                "100% data w/ worst p. prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

plot.e5.xmin <- min(c(utc.var.nom.marg.up1[,1], utc.var.nom.marg.up2[,1], utc.var.nom.marg.up3[,1], 
                      utc.var.nom.marg.up4[,1], utc.var.nom.marg.up5[,1], utc.var.nom.marg.up6[,1], 
                      utc.var.nom.marg.bswp.p[,1]))
plot.e5.xmax <- max(c(utc.var.nom.marg.up1[,1], utc.var.nom.marg.up2[,1], utc.var.nom.marg.up3[,1], 
                      utc.var.nom.marg.up4[,1], utc.var.nom.marg.up5[,1], utc.var.nom.marg.up6[,1], 
                      utc.var.nom.marg.bswp.p[,1]))
plot.e5.ymin <- min(c(utc.var.nom.marg.up1[,2], utc.var.nom.marg.up2[,2], utc.var.nom.marg.up3[,2], 
                      utc.var.nom.marg.up4[,2], utc.var.nom.marg.up5[,2], utc.var.nom.marg.up6[,2], 
                      utc.var.nom.marg.bswp.p[,2]))
plot.e5.ymax <- max(c(utc.var.nom.marg.up1[,2], utc.var.nom.marg.up2[,2], utc.var.nom.marg.up3[,2], 
                      utc.var.nom.marg.up4[,2], utc.var.nom.marg.up5[,2], utc.var.nom.marg.up6[,2], 
                      utc.var.nom.marg.bswp.p[,2]))
plot(utc.var.nom.marg.up1,xlim=c(plot.e5.xmin,plot.e5.xmax),ylim=c(plot.e5.ymin,plot.e5.ymax),
     col="red",t="l",lty=1,xlab=expression(sigma^2),ylab="")
lines(utc.var.nom.marg.up2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.up3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.up4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.up5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.up6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.bswp.p,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ worst p. prior","2% data w/ worst p. prior",
                                "5% data w/ worst p. prior","10% data w/ worst p. prior",
                                "20% data w/ worst p. prior", "50% data w/ worst p. prior",
                                "100% data w/ worst p. prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

plot.e6.xmin <- min(c(utc.range.nom.marg.up1[,1], utc.range.nom.marg.up2[,1], utc.range.nom.marg.up3[,1], 
                      utc.range.nom.marg.up4[,1], utc.range.nom.marg.up5[,1], utc.range.nom.marg.up6[,1], 
                      utc.range.nom.marg.bswp.p[,1]))
plot.e6.xmax <- max(c(utc.range.nom.marg.up1[,1], utc.range.nom.marg.up2[,1], utc.range.nom.marg.up3[,1], 
                      utc.range.nom.marg.up4[,1], utc.range.nom.marg.up5[,1], utc.range.nom.marg.up6[,1], 
                      utc.range.nom.marg.bswp.p[,1]))
plot.e6.ymin <- min(c(utc.range.nom.marg.up1[,2], utc.range.nom.marg.up2[,2], utc.range.nom.marg.up3[,2], 
                      utc.range.nom.marg.up4[,2], utc.range.nom.marg.up5[,2], utc.range.nom.marg.up6[,2], 
                      utc.range.nom.marg.bswp.p[,2]))
plot.e6.ymax <- max(c(utc.range.nom.marg.up1[,2], utc.range.nom.marg.up2[,2], utc.range.nom.marg.up3[,2], 
                      utc.range.nom.marg.up4[,2], utc.range.nom.marg.up5[,2], utc.range.nom.marg.up6[,2], 
                      utc.range.nom.marg.bswp.p[,2]))
plot(utc.range.nom.marg.up1,xlim=c(plot.e6.xmin,plot.e6.xmax),ylim=c(plot.e6.ymin,plot.e6.ymax),
     col="red",t="l",lty=1,xlab="r",ylab="")
lines(utc.range.nom.marg.up2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.up3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.up4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.up5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.up6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.bswp.p,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ worst p. prior","2% data w/ worst p. prior",
                                "5% data w/ worst p. prior","10% data w/ worst p. prior",
                                "20% data w/ worst p. prior", "50% data w/ worst p. prior",
                                "100% data w/ worst p. prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

####################################
# show results with unweighted sampling and vague prior
####################################
plot.f1.xmin <- min(c(utc.fixed1_marg.uv1[,1], utc.fixed1_marg.uv2[,1], utc.fixed1_marg.uv3[,1], 
                      utc.fixed1_marg.uv4[,1], utc.fixed1_marg.uv5[,1], utc.fixed1_marg.uv6[,1], 
                      utc.fixed1_marg.bsvp[,1]))
plot.f1.xmax <- max(c(utc.fixed1_marg.uv1[,1], utc.fixed1_marg.uv2[,1], utc.fixed1_marg.uv3[,1], 
                      utc.fixed1_marg.uv4[,1], utc.fixed1_marg.uv5[,1], utc.fixed1_marg.uv6[,1], 
                      utc.fixed1_marg.bsvp[,1]))
plot.f1.ymin <- min(c(utc.fixed1_marg.uv1[,2], utc.fixed1_marg.uv2[,2], utc.fixed1_marg.uv3[,2], 
                      utc.fixed1_marg.uv4[,2], utc.fixed1_marg.uv5[,2], utc.fixed1_marg.uv6[,2], 
                      utc.fixed1_marg.bsvp[,2]))
plot.f1.ymax <- max(c(utc.fixed1_marg.uv1[,2], utc.fixed1_marg.uv2[,2], utc.fixed1_marg.uv3[,2], 
                      utc.fixed1_marg.uv4[,2], utc.fixed1_marg.uv5[,2], utc.fixed1_marg.uv6[,2], 
                      utc.fixed1_marg.bsvp[,2]))
plot(utc.fixed1_marg.uv1,xlim=c(plot.f1.xmin,plot.f1.xmax),ylim=c(plot.f1.ymin,plot.f1.ymax),
     col="red",t="l",lty=1,xlab=expression(beta[0]),ylab="")
lines(utc.fixed1_marg.uv2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.uv3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.uv4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.uv5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.uv6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.bsvp,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ vague prior","2% data w/ vague prior",
                                "5% data w/ vague prior","10% data w/ vague prior",
                                "20% data w/ vague prior", "50% data w/ vague prior",
                                "100% data w/ vague prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

plot.f2.xmin <- min(c(utc.fixed2_marg.uv1[,1], utc.fixed2_marg.uv2[,1], utc.fixed2_marg.uv3[,1], 
                      utc.fixed2_marg.uv4[,1], utc.fixed2_marg.uv5[,1], utc.fixed2_marg.uv6[,1], 
                      utc.fixed2_marg.bsvp[,1]))
plot.f2.xmax <- max(c(utc.fixed2_marg.uv1[,1], utc.fixed2_marg.uv2[,1], utc.fixed2_marg.uv3[,1], 
                      utc.fixed2_marg.uv4[,1], utc.fixed2_marg.uv5[,1], utc.fixed2_marg.uv6[,1], 
                      utc.fixed2_marg.bsvp[,1]))
plot.f2.ymin <- min(c(utc.fixed2_marg.uv1[,2], utc.fixed2_marg.uv2[,2], utc.fixed2_marg.uv3[,2], 
                      utc.fixed2_marg.uv4[,2], utc.fixed2_marg.uv5[,2], utc.fixed2_marg.uv6[,2], 
                      utc.fixed2_marg.bsvp[,2]))
plot.f2.ymax <- max(c(utc.fixed2_marg.uv1[,2], utc.fixed2_marg.uv2[,2], utc.fixed2_marg.uv3[,2], 
                      utc.fixed2_marg.uv4[,2], utc.fixed2_marg.uv5[,2], utc.fixed2_marg.uv6[,2], 
                      utc.fixed2_marg.bsvp[,2]))
plot(utc.fixed2_marg.uv1,xlim=c(plot.f2.xmin,plot.f2.xmax),ylim=c(plot.f2.ymin,plot.f2.ymax),
     col="red",t="l",lty=1,xlab=expression(beta[1]),ylab="")
lines(utc.fixed2_marg.uv2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.uv3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.uv4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.uv5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.uv6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.bsvp,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ vague prior","2% data w/ vague prior",
                                "5% data w/ vague prior","10% data w/ vague prior",
                                "20% data w/ vague prior", "50% data w/ vague prior",
                                "100% data w/ vague prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

plot.f3.xmin <- min(c(utc.fixed3_marg.uv1[,1], utc.fixed3_marg.uv2[,1], utc.fixed3_marg.uv3[,1], 
                      utc.fixed3_marg.uv4[,1], utc.fixed3_marg.uv5[,1], utc.fixed3_marg.uv6[,1], 
                      utc.fixed3_marg.bsvp[,1]))
plot.f3.xmax <- max(c(utc.fixed3_marg.uv1[,1], utc.fixed3_marg.uv2[,1], utc.fixed3_marg.uv3[,1], 
                      utc.fixed3_marg.uv4[,1], utc.fixed3_marg.uv5[,1], utc.fixed3_marg.uv6[,1], 
                      utc.fixed3_marg.bsvp[,1]))
plot.f3.ymin <- min(c(utc.fixed3_marg.uv1[,2], utc.fixed3_marg.uv2[,2], utc.fixed3_marg.uv3[,2], 
                      utc.fixed3_marg.uv4[,2], utc.fixed3_marg.uv5[,2], utc.fixed3_marg.uv6[,2], 
                      utc.fixed3_marg.bsvp[,2]))
plot.f3.ymax <- max(c(utc.fixed3_marg.uv1[,2], utc.fixed3_marg.uv2[,2], utc.fixed3_marg.uv3[,2], 
                      utc.fixed3_marg.uv4[,2], utc.fixed3_marg.uv5[,2], utc.fixed3_marg.uv6[,2], 
                      utc.fixed3_marg.bsvp[,2]))
plot(utc.fixed3_marg.uv1,xlim=c(plot.f3.xmin,plot.f3.xmax),ylim=c(plot.f3.ymin,plot.f3.ymax),
     col="red",t="l",lty=1,xlab=expression(beta[2]),ylab="")
lines(utc.fixed3_marg.uv2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.uv3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.uv4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.uv5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.uv6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.bsvp,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ vague prior","2% data w/ vague prior",
                                "5% data w/ vague prior","10% data w/ vague prior",
                                "20% data w/ vague prior", "50% data w/ vague prior",
                                "100% data w/ vague prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

plot.f4.xmin <- min(c(utc.sigma2e_marg.uv1[,1], utc.sigma2e_marg.uv2[,1], utc.sigma2e_marg.uv3[,1], 
                      utc.sigma2e_marg.uv4[,1], utc.sigma2e_marg.uv5[,1], utc.sigma2e_marg.uv6[,1], 
                      utc.sigma2e_marg.bsvp[,1]))
plot.f4.xmax <- max(c(utc.sigma2e_marg.uv1[,1], utc.sigma2e_marg.uv2[,1], utc.sigma2e_marg.uv3[,1], 
                      utc.sigma2e_marg.uv4[,1], utc.sigma2e_marg.uv5[,1], utc.sigma2e_marg.uv6[,1], 
                      utc.sigma2e_marg.bsvp[,1]))
plot.f4.ymin <- min(c(utc.sigma2e_marg.uv1[,2], utc.sigma2e_marg.uv2[,2], utc.sigma2e_marg.uv3[,2], 
                      utc.sigma2e_marg.uv4[,2], utc.sigma2e_marg.uv5[,2], utc.sigma2e_marg.uv6[,2], 
                      utc.sigma2e_marg.bsvp[,2]))
plot.f4.ymax <- max(c(utc.sigma2e_marg.uv1[,2], utc.sigma2e_marg.uv2[,2], utc.sigma2e_marg.uv3[,2], 
                      utc.sigma2e_marg.uv4[,2], utc.sigma2e_marg.uv5[,2], utc.sigma2e_marg.uv6[,2], 
                      utc.sigma2e_marg.bsvp[,2]))
plot(utc.sigma2e_marg.uv1,xlim=c(plot.f4.xmin,plot.f4.xmax),ylim=c(plot.f4.ymin,plot.f4.ymax),
     col="red",t="l",lty=1,xlab=expression(sigma[e]^2),ylab="")
lines(utc.sigma2e_marg.uv2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.uv3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.uv4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.uv5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.uv6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.bsvp,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ vague prior","2% data w/ vague prior",
                                "5% data w/ vague prior","10% data w/ vague prior",
                                "20% data w/ vague prior", "50% data w/ vague prior",
                                "100% data w/ vague prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

plot.f5.xmin <- min(c(utc.var.nom.marg.uv1[,1], utc.var.nom.marg.uv2[,1], utc.var.nom.marg.uv3[,1], 
                      utc.var.nom.marg.uv4[,1], utc.var.nom.marg.uv5[,1], utc.var.nom.marg.uv6[,1], 
                      utc.var.nom.marg.bsvp[,1]))
plot.f5.xmax <- max(c(utc.var.nom.marg.uv1[,1], utc.var.nom.marg.uv2[,1], utc.var.nom.marg.uv3[,1], 
                      utc.var.nom.marg.uv4[,1], utc.var.nom.marg.uv5[,1], utc.var.nom.marg.uv6[,1], 
                      utc.var.nom.marg.bsvp[,1]))
plot.f5.ymin <- min(c(utc.var.nom.marg.uv1[,2], utc.var.nom.marg.uv2[,2], utc.var.nom.marg.uv3[,2], 
                      utc.var.nom.marg.uv4[,2], utc.var.nom.marg.uv5[,2], utc.var.nom.marg.uv6[,2], 
                      utc.var.nom.marg.bsvp[,2]))
plot.f5.ymax <- max(c(utc.var.nom.marg.uv1[,2], utc.var.nom.marg.uv2[,2], utc.var.nom.marg.uv3[,2], 
                      utc.var.nom.marg.uv4[,2], utc.var.nom.marg.uv5[,2], utc.var.nom.marg.uv6[,2], 
                      utc.var.nom.marg.bsvp[,2]))
plot(utc.var.nom.marg.uv1,xlim=c(plot.f5.xmin,plot.f5.xmax),ylim=c(plot.f5.ymin,plot.f5.ymax),
     col="red",t="l",lty=1,xlab=expression(sigma^2),ylab="")
lines(utc.var.nom.marg.uv2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.uv3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.uv4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.uv5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.uv6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.bsvp,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ vague prior","2% data w/ vague prior",
                                "5% data w/ vague prior","10% data w/ vague prior",
                                "20% data w/ vague prior", "50% data w/ vague prior",
                                "100% data w/ vague prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

plot.f6.xmin <- min(c(utc.range.nom.marg.uv1[,1], utc.range.nom.marg.uv2[,1], utc.range.nom.marg.uv3[,1], 
                      utc.range.nom.marg.uv4[,1], utc.range.nom.marg.uv5[,1], utc.range.nom.marg.uv6[,1], 
                      utc.range.nom.marg.bsvp[,1]))
plot.f6.xmax <- max(c(utc.range.nom.marg.uv1[,1], utc.range.nom.marg.uv2[,1], utc.range.nom.marg.uv3[,1], 
                      utc.range.nom.marg.uv4[,1], utc.range.nom.marg.uv5[,1], utc.range.nom.marg.uv6[,1], 
                      utc.range.nom.marg.bsvp[,1]))
plot.f6.ymin <- min(c(utc.range.nom.marg.uv1[,2], utc.range.nom.marg.uv2[,2], utc.range.nom.marg.uv3[,2], 
                      utc.range.nom.marg.uv4[,2], utc.range.nom.marg.uv5[,2], utc.range.nom.marg.uv6[,2], 
                      utc.range.nom.marg.bsvp[,2]))
plot.f6.ymax <- max(c(utc.range.nom.marg.uv1[,2], utc.range.nom.marg.uv2[,2], utc.range.nom.marg.uv3[,2], 
                      utc.range.nom.marg.uv4[,2], utc.range.nom.marg.uv5[,2], utc.range.nom.marg.uv6[,2], 
                      utc.range.nom.marg.bsvp[,2]))
plot(utc.range.nom.marg.uv1,xlim=c(plot.f6.xmin,plot.f6.xmax),ylim=c(plot.f6.ymin,plot.f6.ymax),
     col="red",t="l",lty=1,xlab="r",ylab="")
lines(utc.range.nom.marg.uv2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.uv3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.uv4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.uv5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.uv6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.bsvp,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ vague prior","2% data w/ vague prior",
                                "5% data w/ vague prior","10% data w/ vague prior",
                                "20% data w/ vague prior", "50% data w/ vague prior",
                                "100% data w/ vague prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

####################################
# show results with weighted sampling and normal prior
####################################
plot.g1.xmin <- min(c(utc.fixed1_marg.wn1[,1], utc.fixed1_marg.wn2[,1], utc.fixed1_marg.wn3[,1], 
                      utc.fixed1_marg.wn4[,1], utc.fixed1_marg.wn5[,1], utc.fixed1_marg.wn6[,1], 
                      utc.fixed1_marg.bswp.n[,1]))
plot.g1.xmax <- max(c(utc.fixed1_marg.wn1[,1], utc.fixed1_marg.wn2[,1], utc.fixed1_marg.wn3[,1], 
                      utc.fixed1_marg.wn4[,1], utc.fixed1_marg.wn5[,1], utc.fixed1_marg.wn6[,1], 
                      utc.fixed1_marg.bswp.n[,1]))
plot.g1.ymin <- min(c(utc.fixed1_marg.wn1[,2], utc.fixed1_marg.wn2[,2], utc.fixed1_marg.wn3[,2], 
                      utc.fixed1_marg.wn4[,2], utc.fixed1_marg.wn5[,2], utc.fixed1_marg.wn6[,2], 
                      utc.fixed1_marg.bswp.n[,2]))
plot.g1.ymax <- max(c(utc.fixed1_marg.wn1[,2], utc.fixed1_marg.wn2[,2], utc.fixed1_marg.wn3[,2], 
                      utc.fixed1_marg.wn4[,2], utc.fixed1_marg.wn5[,2], utc.fixed1_marg.wn6[,2], 
                      utc.fixed1_marg.bswp.n[,2]))
plot(utc.fixed1_marg.wn1,xlim=c(plot.g1.xmin,plot.g1.xmax),ylim=c(plot.g1.ymin,plot.g1.ymax),
     col="red",t="l",lty=1,xlab=expression(beta[0]),ylab="")
lines(utc.fixed1_marg.wn2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.wn3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.wn4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.wn5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.wn6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.bswp.n,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ worst n. prior","2% data w/ worst n. prior",
                                "5% data w/ worst n. prior","10% data w/ worst n. prior",
                                "20% data w/ worst n. prior", "50% data w/ worst n. prior",
                                "100% data w/ worst n. prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

plot.g2.xmin <- min(c(utc.fixed2_marg.wn1[,1], utc.fixed2_marg.wn2[,1], utc.fixed2_marg.wn3[,1], 
                      utc.fixed2_marg.wn4[,1], utc.fixed2_marg.wn5[,1], utc.fixed2_marg.wn6[,1], 
                      utc.fixed2_marg.bswp.n[,1]))
plot.g2.xmax <- max(c(utc.fixed2_marg.wn1[,1], utc.fixed2_marg.wn2[,1], utc.fixed2_marg.wn3[,1], 
                      utc.fixed2_marg.wn4[,1], utc.fixed2_marg.wn5[,1], utc.fixed2_marg.wn6[,1], 
                      utc.fixed2_marg.bswp.n[,1]))
plot.g2.ymin <- min(c(utc.fixed2_marg.wn1[,2], utc.fixed2_marg.wn2[,2], utc.fixed2_marg.wn3[,2], 
                      utc.fixed2_marg.wn4[,2], utc.fixed2_marg.wn5[,2], utc.fixed2_marg.wn6[,2], 
                      utc.fixed2_marg.bswp.n[,2]))
plot.g2.ymax <- max(c(utc.fixed2_marg.wn1[,2], utc.fixed2_marg.wn2[,2], utc.fixed2_marg.wn3[,2], 
                      utc.fixed2_marg.wn4[,2], utc.fixed2_marg.wn5[,2], utc.fixed2_marg.wn6[,2], 
                      utc.fixed2_marg.bswp.n[,2]))
plot(utc.fixed2_marg.wn1,xlim=c(plot.g2.xmin,plot.g2.xmax),ylim=c(plot.g2.ymin,plot.g2.ymax),
     col="red",t="l",lty=1,xlab=expression(beta[1]),ylab="")
lines(utc.fixed2_marg.wn2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.wn3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.wn4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.wn5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.wn6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.bswp.n,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ worst n. prior","2% data w/ worst n. prior",
                                "5% data w/ worst n. prior","10% data w/ worst n. prior",
                                "20% data w/ worst n. prior", "50% data w/ worst n. prior",
                                "100% data w/ worst n. prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

plot.g3.xmin <- min(c(utc.fixed3_marg.wn1[,1], utc.fixed3_marg.wn2[,1], utc.fixed3_marg.wn3[,1], 
                      utc.fixed3_marg.wn4[,1], utc.fixed3_marg.wn5[,1], utc.fixed3_marg.wn6[,1], 
                      utc.fixed3_marg.bswp.n[,1]))
plot.g3.xmax <- max(c(utc.fixed3_marg.wn1[,1], utc.fixed3_marg.wn2[,1], utc.fixed3_marg.wn3[,1], 
                      utc.fixed3_marg.wn4[,1], utc.fixed3_marg.wn5[,1], utc.fixed3_marg.wn6[,1], 
                      utc.fixed3_marg.bswp.n[,1]))
plot.g3.ymin <- min(c(utc.fixed3_marg.wn1[,2], utc.fixed3_marg.wn2[,2], utc.fixed3_marg.wn3[,2], 
                      utc.fixed3_marg.wn4[,2], utc.fixed3_marg.wn5[,2], utc.fixed3_marg.wn6[,2], 
                      utc.fixed3_marg.bswp.n[,2]))
plot.g3.ymax <- max(c(utc.fixed3_marg.wn1[,2], utc.fixed3_marg.wn2[,2], utc.fixed3_marg.wn3[,2], 
                      utc.fixed3_marg.wn4[,2], utc.fixed3_marg.wn5[,2], utc.fixed3_marg.wn6[,2], 
                      utc.fixed3_marg.bswp.n[,2]))
plot(utc.fixed3_marg.wn1,xlim=c(plot.g3.xmin,plot.g3.xmax),ylim=c(plot.g3.ymin,plot.g3.ymax),
     col="red",t="l",lty=1,xlab=expression(beta[2]),ylab="")
lines(utc.fixed3_marg.wn2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.wn3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.wn4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.wn5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.wn6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.bswp.n,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ worst n. prior","2% data w/ worst n. prior",
                                "5% data w/ worst n. prior","10% data w/ worst n. prior",
                                "20% data w/ worst n. prior", "50% data w/ worst n. prior",
                                "100% data w/ worst n. prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

plot.g4.xmin <- min(c(utc.sigma2e_marg.wn1[,1], utc.sigma2e_marg.wn2[,1], utc.sigma2e_marg.wn3[,1], 
                      utc.sigma2e_marg.wn4[,1], utc.sigma2e_marg.wn5[,1], utc.sigma2e_marg.wn6[,1], 
                      utc.sigma2e_marg.bswp.n[,1]))
plot.g4.xmax <- max(c(utc.sigma2e_marg.wn1[,1], utc.sigma2e_marg.wn2[,1], utc.sigma2e_marg.wn3[,1], 
                      utc.sigma2e_marg.wn4[,1], utc.sigma2e_marg.wn5[,1], utc.sigma2e_marg.wn6[,1], 
                      utc.sigma2e_marg.bswp.n[,1]))
plot.g4.ymin <- min(c(utc.sigma2e_marg.wn1[,2], utc.sigma2e_marg.wn2[,2], utc.sigma2e_marg.wn3[,2], 
                      utc.sigma2e_marg.wn4[,2], utc.sigma2e_marg.wn5[,2], utc.sigma2e_marg.wn6[,2], 
                      utc.sigma2e_marg.bswp.n[,2]))
plot.g4.ymax <- max(c(utc.sigma2e_marg.wn1[,2], utc.sigma2e_marg.wn2[,2], utc.sigma2e_marg.wn3[,2], 
                      utc.sigma2e_marg.wn4[,2], utc.sigma2e_marg.wn5[,2], utc.sigma2e_marg.wn6[,2], 
                      utc.sigma2e_marg.bswp.n[,2]))
plot(utc.sigma2e_marg.wn1,xlim=c(plot.g4.xmin,plot.g4.xmax),ylim=c(plot.g4.ymin,plot.g4.ymax),
     col="red",t="l",lty=1,xlab=expression(sigma[e]^2),ylab="")
lines(utc.sigma2e_marg.wn2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.wn3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.wn4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.wn5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.wn6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.bswp.n,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ worst n. prior","2% data w/ worst n. prior",
                                "5% data w/ worst n. prior","10% data w/ worst n. prior",
                                "20% data w/ worst n. prior", "50% data w/ worst n. prior",
                                "100% data w/ worst n. prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

plot.g5.xmin <- min(c(utc.var.nom.marg.wn1[,1], utc.var.nom.marg.wn2[,1], utc.var.nom.marg.wn3[,1], 
                      utc.var.nom.marg.wn4[,1], utc.var.nom.marg.wn5[,1], utc.var.nom.marg.wn6[,1], 
                      utc.var.nom.marg.bswp.n[,1]))
plot.g5.xmax <- max(c(utc.var.nom.marg.wn1[,1], utc.var.nom.marg.wn2[,1], utc.var.nom.marg.wn3[,1], 
                      utc.var.nom.marg.wn4[,1], utc.var.nom.marg.wn5[,1], utc.var.nom.marg.wn6[,1], 
                      utc.var.nom.marg.bswp.n[,1]))
plot.g5.ymin <- min(c(utc.var.nom.marg.wn1[,2], utc.var.nom.marg.wn2[,2], utc.var.nom.marg.wn3[,2], 
                      utc.var.nom.marg.wn4[,2], utc.var.nom.marg.wn5[,2], utc.var.nom.marg.wn6[,2], 
                      utc.var.nom.marg.bswp.n[,2]))
plot.g5.ymax <- max(c(utc.var.nom.marg.wn1[,2], utc.var.nom.marg.wn2[,2], utc.var.nom.marg.wn3[,2], 
                      utc.var.nom.marg.wn4[,2], utc.var.nom.marg.wn5[,2], utc.var.nom.marg.wn6[,2], 
                      utc.var.nom.marg.bswp.n[,2]))
plot(utc.var.nom.marg.wn1,xlim=c(plot.g5.xmin,plot.g5.xmax),ylim=c(plot.g5.ymin,plot.g5.ymax),
     col="red",t="l",lty=1,xlab=expression(sigma^2),ylab="")
lines(utc.var.nom.marg.wn2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.wn3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.wn4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.wn5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.wn6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.bswp.n,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ worst n. prior","2% data w/ worst n. prior",
                                "5% data w/ worst n. prior","10% data w/ worst n. prior",
                                "20% data w/ worst n. prior", "50% data w/ worst n. prior",
                                "100% data w/ worst n. prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

plot.g6.xmin <- min(c(utc.range.nom.marg.wn1[,1], utc.range.nom.marg.wn2[,1], utc.range.nom.marg.wn3[,1], 
                      utc.range.nom.marg.wn4[,1], utc.range.nom.marg.wn5[,1], utc.range.nom.marg.wn6[,1], 
                      utc.range.nom.marg.bswp.n[,1]))
plot.g6.xmax <- max(c(utc.range.nom.marg.wn1[,1], utc.range.nom.marg.wn2[,1], utc.range.nom.marg.wn3[,1], 
                      utc.range.nom.marg.wn4[,1], utc.range.nom.marg.wn5[,1], utc.range.nom.marg.wn6[,1], 
                      utc.range.nom.marg.bswp.n[,1]))
plot.g6.ymin <- min(c(utc.range.nom.marg.wn1[,2], utc.range.nom.marg.wn2[,2], utc.range.nom.marg.wn3[,2], 
                      utc.range.nom.marg.wn4[,2], utc.range.nom.marg.wn5[,2], utc.range.nom.marg.wn6[,2], 
                      utc.range.nom.marg.bswp.n[,2]))
plot.g6.ymax <- max(c(utc.range.nom.marg.wn1[,2], utc.range.nom.marg.wn2[,2], utc.range.nom.marg.wn3[,2], 
                      utc.range.nom.marg.wn4[,2], utc.range.nom.marg.wn5[,2], utc.range.nom.marg.wn6[,2], 
                      utc.range.nom.marg.bswp.n[,2]))
plot(utc.range.nom.marg.wn1,xlim=c(plot.g6.xmin,plot.g6.xmax),ylim=c(plot.g6.ymin,plot.g6.ymax),
     col="red",t="l",lty=1,xlab="r",ylab="")
lines(utc.range.nom.marg.wn2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.wn3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.wn4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.wn5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.wn6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.bswp.n,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ worst n. prior","2% data w/ worst n. prior",
                                "5% data w/ worst n. prior","10% data w/ worst n. prior",
                                "20% data w/ worst n. prior", "50% data w/ worst n. prior",
                                "100% data w/ worst n. prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

####################################
# show results with weighted sampling and pc prior
####################################
plot.h1.xmin <- min(c(utc.fixed1_marg.wp1[,1], utc.fixed1_marg.wp2[,1], utc.fixed1_marg.wp3[,1], 
                      utc.fixed1_marg.wp4[,1], utc.fixed1_marg.wp5[,1], utc.fixed1_marg.wp6[,1], 
                      utc.fixed1_marg.bswp.p[,1]))
plot.h1.xmax <- max(c(utc.fixed1_marg.wp1[,1], utc.fixed1_marg.wp2[,1], utc.fixed1_marg.wp3[,1], 
                      utc.fixed1_marg.wp4[,1], utc.fixed1_marg.wp5[,1], utc.fixed1_marg.wp6[,1], 
                      utc.fixed1_marg.bswp.p[,1]))
plot.h1.ymin <- min(c(utc.fixed1_marg.wp1[,2], utc.fixed1_marg.wp2[,2], utc.fixed1_marg.wp3[,2], 
                      utc.fixed1_marg.wp4[,2], utc.fixed1_marg.wp5[,2], utc.fixed1_marg.wp6[,2], 
                      utc.fixed1_marg.bswp.p[,2]))
plot.h1.ymax <- max(c(utc.fixed1_marg.wp1[,2], utc.fixed1_marg.wp2[,2], utc.fixed1_marg.wp3[,2], 
                      utc.fixed1_marg.wp4[,2], utc.fixed1_marg.wp5[,2], utc.fixed1_marg.wp6[,2], 
                      utc.fixed1_marg.bswp.p[,2]))
plot(utc.fixed1_marg.wp1,xlim=c(plot.h1.xmin,plot.h1.xmax),ylim=c(plot.h1.ymin,plot.h1.ymax),
     col="red",t="l",lty=1,xlab=expression(beta[0]),ylab="")
lines(utc.fixed1_marg.wp2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.wp3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.wp4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.wp5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.wp6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.bswp.p,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ worst p. prior","2% data w/ worst p. prior",
                                "5% data w/ worst p. prior","10% data w/ worst p. prior",
                                "20% data w/ worst p. prior", "50% data w/ worst p. prior",
                                "100% data w/ worst p. prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

plot.h2.xmin <- min(c(utc.fixed2_marg.wp1[,1], utc.fixed2_marg.wp2[,1], utc.fixed2_marg.wp3[,1], 
                      utc.fixed2_marg.wp4[,1], utc.fixed2_marg.wp5[,1], utc.fixed2_marg.wp6[,1], 
                      utc.fixed2_marg.bswp.p[,1]))
plot.h2.xmax <- max(c(utc.fixed2_marg.wp1[,1], utc.fixed2_marg.wp2[,1], utc.fixed2_marg.wp3[,1], 
                      utc.fixed2_marg.wp4[,1], utc.fixed2_marg.wp5[,1], utc.fixed2_marg.wp6[,1], 
                      utc.fixed2_marg.bswp.p[,1]))
plot.h2.ymin <- min(c(utc.fixed2_marg.wp1[,2], utc.fixed2_marg.wp2[,2], utc.fixed2_marg.wp3[,2], 
                      utc.fixed2_marg.wp4[,2], utc.fixed2_marg.wp5[,2], utc.fixed2_marg.wp6[,2], 
                      utc.fixed2_marg.bswp.p[,2]))
plot.h2.ymax <- max(c(utc.fixed2_marg.wp1[,2], utc.fixed2_marg.wp2[,2], utc.fixed2_marg.wp3[,2], 
                      utc.fixed2_marg.wp4[,2], utc.fixed2_marg.wp5[,2], utc.fixed2_marg.wp6[,2], 
                      utc.fixed2_marg.bswp.p[,2]))
plot(utc.fixed2_marg.wp1,xlim=c(plot.h2.xmin,plot.h2.xmax),ylim=c(plot.h2.ymin,plot.h2.ymax),
     col="red",t="l",lty=1,xlab=expression(beta[1]),ylab="")
lines(utc.fixed2_marg.wp2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.wp3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.wp4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.wp5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.wp6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.bswp.p,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ worst p. prior","2% data w/ worst p. prior",
                                "5% data w/ worst p. prior","10% data w/ worst p. prior",
                                "20% data w/ worst p. prior", "50% data w/ worst p. prior",
                                "100% data w/ worst p. prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

plot.h3.xmin <- min(c(utc.fixed3_marg.wp1[,1], utc.fixed3_marg.wp2[,1], utc.fixed3_marg.wp3[,1], 
                      utc.fixed3_marg.wp4[,1], utc.fixed3_marg.wp5[,1], utc.fixed3_marg.wp6[,1], 
                      utc.fixed3_marg.bswp.p[,1]))
plot.h3.xmax <- max(c(utc.fixed3_marg.wp1[,1], utc.fixed3_marg.wp2[,1], utc.fixed3_marg.wp3[,1], 
                      utc.fixed3_marg.wp4[,1], utc.fixed3_marg.wp5[,1], utc.fixed3_marg.wp6[,1], 
                      utc.fixed3_marg.bswp.p[,1]))
plot.h3.ymin <- min(c(utc.fixed3_marg.wp1[,2], utc.fixed3_marg.wp2[,2], utc.fixed3_marg.wp3[,2], 
                      utc.fixed3_marg.wp4[,2], utc.fixed3_marg.wp5[,2], utc.fixed3_marg.wp6[,2], 
                      utc.fixed3_marg.bswp.p[,2]))
plot.h3.ymax <- max(c(utc.fixed3_marg.wp1[,2], utc.fixed3_marg.wp2[,2], utc.fixed3_marg.wp3[,2], 
                      utc.fixed3_marg.wp4[,2], utc.fixed3_marg.wp5[,2], utc.fixed3_marg.wp6[,2], 
                      utc.fixed3_marg.bswp.p[,2]))
plot(utc.fixed3_marg.wp1,xlim=c(plot.h3.xmin,plot.h3.xmax),ylim=c(plot.h3.ymin,plot.h3.ymax),
     col="red",t="l",lty=1,xlab=expression(beta[2]),ylab="")
lines(utc.fixed3_marg.wp2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.wp3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.wp4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.wp5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.wp6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.bswp.p,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ worst p. prior","2% data w/ worst p. prior",
                                "5% data w/ worst p. prior","10% data w/ worst p. prior",
                                "20% data w/ worst p. prior", "50% data w/ worst p. prior",
                                "100% data w/ worst p. prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

plot.h4.xmin <- min(c(utc.sigma2e_marg.wp1[,1], utc.sigma2e_marg.wp2[,1], utc.sigma2e_marg.wp3[,1], 
                      utc.sigma2e_marg.wp4[,1], utc.sigma2e_marg.wp5[,1], utc.sigma2e_marg.wp6[,1], 
                      utc.sigma2e_marg.bswp.p[,1]))
plot.h4.xmax <- max(c(utc.sigma2e_marg.wp1[,1], utc.sigma2e_marg.wp2[,1], utc.sigma2e_marg.wp3[,1], 
                      utc.sigma2e_marg.wp4[,1], utc.sigma2e_marg.wp5[,1], utc.sigma2e_marg.wp6[,1], 
                      utc.sigma2e_marg.bswp.p[,1]))
plot.h4.ymin <- min(c(utc.sigma2e_marg.wp1[,2], utc.sigma2e_marg.wp2[,2], utc.sigma2e_marg.wp3[,2], 
                      utc.sigma2e_marg.wp4[,2], utc.sigma2e_marg.wp5[,2], utc.sigma2e_marg.wp6[,2], 
                      utc.sigma2e_marg.bswp.p[,2]))
plot.h4.ymax <- max(c(utc.sigma2e_marg.wp1[,2], utc.sigma2e_marg.wp2[,2], utc.sigma2e_marg.wp3[,2], 
                      utc.sigma2e_marg.wp4[,2], utc.sigma2e_marg.wp5[,2], utc.sigma2e_marg.wp6[,2], 
                      utc.sigma2e_marg.bswp.p[,2]))
plot(utc.sigma2e_marg.wp1,xlim=c(plot.h4.xmin,plot.h4.xmax),ylim=c(plot.h4.ymin,plot.h4.ymax),
     col="red",t="l",lty=1,xlab=expression(sigma[e]^2),ylab="")
lines(utc.sigma2e_marg.wp2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.wp3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.wp4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.wp5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.wp6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.bswp.p,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ worst p. prior","2% data w/ worst p. prior",
                                "5% data w/ worst p. prior","10% data w/ worst p. prior",
                                "20% data w/ worst p. prior", "50% data w/ worst p. prior",
                                "100% data w/ worst p. prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

plot.h5.xmin <- min(c(utc.var.nom.marg.wp1[,1], utc.var.nom.marg.wp2[,1], utc.var.nom.marg.wp3[,1], 
                      utc.var.nom.marg.wp4[,1], utc.var.nom.marg.wp5[,1], utc.var.nom.marg.wp6[,1], 
                      utc.var.nom.marg.bswp.p[,1]))
plot.h5.xmax <- max(c(utc.var.nom.marg.wp1[,1], utc.var.nom.marg.wp2[,1], utc.var.nom.marg.wp3[,1], 
                      utc.var.nom.marg.wp4[,1], utc.var.nom.marg.wp5[,1], utc.var.nom.marg.wp6[,1], 
                      utc.var.nom.marg.bswp.p[,1]))
plot.h5.ymin <- min(c(utc.var.nom.marg.wp1[,2], utc.var.nom.marg.wp2[,2], utc.var.nom.marg.wp3[,2], 
                      utc.var.nom.marg.wp4[,2], utc.var.nom.marg.wp5[,2], utc.var.nom.marg.wp6[,2], 
                      utc.var.nom.marg.bswp.p[,2]))
plot.h5.ymax <- max(c(utc.var.nom.marg.wp1[,2], utc.var.nom.marg.wp2[,2], utc.var.nom.marg.wp3[,2], 
                      utc.var.nom.marg.wp4[,2], utc.var.nom.marg.wp5[,2], utc.var.nom.marg.wp6[,2], 
                      utc.var.nom.marg.bswp.p[,2]))
plot(utc.var.nom.marg.wp1,xlim=c(plot.h5.xmin,plot.h5.xmax),ylim=c(plot.h5.ymin,plot.h5.ymax),
     col="red",t="l",lty=1,xlab=expression(sigma^2),ylab="")
lines(utc.var.nom.marg.wp2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.wp3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.wp4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.wp5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.wp6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.bswp.p,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ worst p. prior","2% data w/ worst p. prior",
                                "5% data w/ worst p. prior","10% data w/ worst p. prior",
                                "20% data w/ worst p. prior", "50% data w/ worst p. prior",
                                "100% data w/ worst p. prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

plot.h6.xmin <- min(c(utc.range.nom.marg.wp1[,1], utc.range.nom.marg.wp2[,1], utc.range.nom.marg.wp3[,1], 
                      utc.range.nom.marg.wp4[,1], utc.range.nom.marg.wp5[,1], utc.range.nom.marg.wp6[,1], 
                      utc.range.nom.marg.bswp.p[,1]))
plot.h6.xmax <- max(c(utc.range.nom.marg.wp1[,1], utc.range.nom.marg.wp2[,1], utc.range.nom.marg.wp3[,1], 
                      utc.range.nom.marg.wp4[,1], utc.range.nom.marg.wp5[,1], utc.range.nom.marg.wp6[,1], 
                      utc.range.nom.marg.bswp.p[,1]))
plot.h6.ymin <- min(c(utc.range.nom.marg.wp1[,2], utc.range.nom.marg.wp2[,2], utc.range.nom.marg.wp3[,2], 
                      utc.range.nom.marg.wp4[,2], utc.range.nom.marg.wp5[,2], utc.range.nom.marg.wp6[,2], 
                      utc.range.nom.marg.bswp.p[,2]))
plot.h6.ymax <- max(c(utc.range.nom.marg.wp1[,2], utc.range.nom.marg.wp2[,2], utc.range.nom.marg.wp3[,2], 
                      utc.range.nom.marg.wp4[,2], utc.range.nom.marg.wp5[,2], utc.range.nom.marg.wp6[,2], 
                      utc.range.nom.marg.bswp.p[,2]))
plot(utc.range.nom.marg.wp1,xlim=c(plot.h6.xmin,plot.h6.xmax),ylim=c(plot.h6.ymin,plot.h6.ymax),
     col="red",t="l",lty=1,xlab="r",ylab="")
lines(utc.range.nom.marg.wp2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.wp3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.wp4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.wp5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.wp6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.bswp.p,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ worst p. prior","2% data w/ worst p. prior",
                                "5% data w/ worst p. prior","10% data w/ worst p. prior",
                                "20% data w/ worst p. prior", "50% data w/ worst p. prior",
                                "100% data w/ worst p. prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

####################################
# show results with weighted sampling and vague prior
####################################
plot.i1.xmin <- min(c(utc.fixed1_marg.wv1[,1], utc.fixed1_marg.wv2[,1], utc.fixed1_marg.wv3[,1], 
                      utc.fixed1_marg.wv4[,1], utc.fixed1_marg.wv5[,1], utc.fixed1_marg.wv6[,1], 
                      utc.fixed1_marg.bsvp[,1]))
plot.i1.xmax <- max(c(utc.fixed1_marg.wv1[,1], utc.fixed1_marg.wv2[,1], utc.fixed1_marg.wv3[,1], 
                      utc.fixed1_marg.wv4[,1], utc.fixed1_marg.wv5[,1], utc.fixed1_marg.wv6[,1], 
                      utc.fixed1_marg.bsvp[,1]))
plot.i1.ymin <- min(c(utc.fixed1_marg.wv1[,2], utc.fixed1_marg.wv2[,2], utc.fixed1_marg.wv3[,2], 
                      utc.fixed1_marg.wv4[,2], utc.fixed1_marg.wv5[,2], utc.fixed1_marg.wv6[,2], 
                      utc.fixed1_marg.bsvp[,2]))
plot.i1.ymax <- max(c(utc.fixed1_marg.wv1[,2], utc.fixed1_marg.wv2[,2], utc.fixed1_marg.wv3[,2], 
                      utc.fixed1_marg.wv4[,2], utc.fixed1_marg.wv5[,2], utc.fixed1_marg.wv6[,2], 
                      utc.fixed1_marg.bsvp[,2]))
plot(utc.fixed1_marg.wv1,xlim=c(plot.i1.xmin,plot.i1.xmax),ylim=c(plot.i1.ymin,plot.i1.ymax),
     col="red",t="l",lty=1,xlab=expression(beta[0]),ylab="")
lines(utc.fixed1_marg.wv2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.wv3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.wv4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.wv5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.wv6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed1_marg.bsvp,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ vague prior","2% data w/ vague prior",
                                "5% data w/ vague prior","10% data w/ vague prior",
                                "20% data w/ vague prior", "50% data w/ vague prior",
                                "100% data w/ vague prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

plot.i2.xmin <- min(c(utc.fixed2_marg.wv1[,1], utc.fixed2_marg.wv2[,1], utc.fixed2_marg.wv3[,1], 
                      utc.fixed2_marg.wv4[,1], utc.fixed2_marg.wv5[,1], utc.fixed2_marg.wv6[,1], 
                      utc.fixed2_marg.bsvp[,1]))
plot.i2.xmax <- max(c(utc.fixed2_marg.wv1[,1], utc.fixed2_marg.wv2[,1], utc.fixed2_marg.wv3[,1], 
                      utc.fixed2_marg.wv4[,1], utc.fixed2_marg.wv5[,1], utc.fixed2_marg.wv6[,1], 
                      utc.fixed2_marg.bsvp[,1]))
plot.i2.ymin <- min(c(utc.fixed2_marg.wv1[,2], utc.fixed2_marg.wv2[,2], utc.fixed2_marg.wv3[,2], 
                      utc.fixed2_marg.wv4[,2], utc.fixed2_marg.wv5[,2], utc.fixed2_marg.wv6[,2], 
                      utc.fixed2_marg.bsvp[,2]))
plot.i2.ymax <- max(c(utc.fixed2_marg.wv1[,2], utc.fixed2_marg.wv2[,2], utc.fixed2_marg.wv3[,2], 
                      utc.fixed2_marg.wv4[,2], utc.fixed2_marg.wv5[,2], utc.fixed2_marg.wv6[,2], 
                      utc.fixed2_marg.bsvp[,2]))
plot(utc.fixed2_marg.wv1,xlim=c(plot.i2.xmin,plot.i2.xmax),ylim=c(plot.i2.ymin,plot.i2.ymax),
     col="red",t="l",lty=1,xlab=expression(beta[1]),ylab="")
lines(utc.fixed2_marg.wv2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.wv3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.wv4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.wv5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.wv6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed2_marg.bsvp,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ vague prior","2% data w/ vague prior",
                                "5% data w/ vague prior","10% data w/ vague prior",
                                "20% data w/ vague prior", "50% data w/ vague prior",
                                "100% data w/ vague prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

plot.i3.xmin <- min(c(utc.fixed3_marg.wv1[,1], utc.fixed3_marg.wv2[,1], utc.fixed3_marg.wv3[,1], 
                      utc.fixed3_marg.wv4[,1], utc.fixed3_marg.wv5[,1], utc.fixed3_marg.wv6[,1], 
                      utc.fixed3_marg.bsvp[,1]))
plot.i3.xmax <- max(c(utc.fixed3_marg.wv1[,1], utc.fixed3_marg.wv2[,1], utc.fixed3_marg.wv3[,1], 
                      utc.fixed3_marg.wv4[,1], utc.fixed3_marg.wv5[,1], utc.fixed3_marg.wv6[,1], 
                      utc.fixed3_marg.bsvp[,1]))
plot.i3.ymin <- min(c(utc.fixed3_marg.wv1[,2], utc.fixed3_marg.wv2[,2], utc.fixed3_marg.wv3[,2], 
                      utc.fixed3_marg.wv4[,2], utc.fixed3_marg.wv5[,2], utc.fixed3_marg.wv6[,2], 
                      utc.fixed3_marg.bsvp[,2]))
plot.i3.ymax <- max(c(utc.fixed3_marg.wv1[,2], utc.fixed3_marg.wv2[,2], utc.fixed3_marg.wv3[,2], 
                      utc.fixed3_marg.wv4[,2], utc.fixed3_marg.wv5[,2], utc.fixed3_marg.wv6[,2], 
                      utc.fixed3_marg.bsvp[,2]))
plot(utc.fixed3_marg.wv1,xlim=c(plot.i3.xmin,plot.i3.xmax),ylim=c(plot.i3.ymin,plot.i3.ymax),
     col="red",t="l",lty=1,xlab=expression(beta[2]),ylab="")
lines(utc.fixed3_marg.wv2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.wv3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.wv4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.wv5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.wv6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.fixed3_marg.bsvp,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ vague prior","2% data w/ vague prior",
                                "5% data w/ vague prior","10% data w/ vague prior",
                                "20% data w/ vague prior", "50% data w/ vague prior",
                                "100% data w/ vague prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

plot.i4.xmin <- min(c(utc.sigma2e_marg.wv1[,1], utc.sigma2e_marg.wv2[,1], utc.sigma2e_marg.wv3[,1], 
                      utc.sigma2e_marg.wv4[,1], utc.sigma2e_marg.wv5[,1], utc.sigma2e_marg.wv6[,1], 
                      utc.sigma2e_marg.bsvp[,1]))
plot.i4.xmax <- max(c(utc.sigma2e_marg.wv1[,1], utc.sigma2e_marg.wv2[,1], utc.sigma2e_marg.wv3[,1], 
                      utc.sigma2e_marg.wv4[,1], utc.sigma2e_marg.wv5[,1], utc.sigma2e_marg.wv6[,1], 
                      utc.sigma2e_marg.bsvp[,1]))
plot.i4.ymin <- min(c(utc.sigma2e_marg.wv1[,2], utc.sigma2e_marg.wv2[,2], utc.sigma2e_marg.wv3[,2], 
                      utc.sigma2e_marg.wv4[,2], utc.sigma2e_marg.wv5[,2], utc.sigma2e_marg.wv6[,2], 
                      utc.sigma2e_marg.bsvp[,2]))
plot.i4.ymax <- max(c(utc.sigma2e_marg.wv1[,2], utc.sigma2e_marg.wv2[,2], utc.sigma2e_marg.wv3[,2], 
                      utc.sigma2e_marg.wv4[,2], utc.sigma2e_marg.wv5[,2], utc.sigma2e_marg.wv6[,2], 
                      utc.sigma2e_marg.bsvp[,2]))
plot(utc.sigma2e_marg.wv1,xlim=c(plot.i4.xmin,plot.i4.xmax),ylim=c(plot.i4.ymin,plot.i4.ymax),
     col="red",t="l",lty=1,xlab=expression(sigma[e]^2),ylab="")
lines(utc.sigma2e_marg.wv2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.wv3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.wv4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.wv5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.wv6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.sigma2e_marg.bsvp,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ vague prior","2% data w/ vague prior",
                                "5% data w/ vague prior","10% data w/ vague prior",
                                "20% data w/ vague prior", "50% data w/ vague prior",
                                "100% data w/ vague prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

plot.i5.xmin <- min(c(utc.var.nom.marg.wv1[,1], utc.var.nom.marg.wv2[,1], utc.var.nom.marg.wv3[,1], 
                      utc.var.nom.marg.wv4[,1], utc.var.nom.marg.wv5[,1], utc.var.nom.marg.wv6[,1], 
                      utc.var.nom.marg.bsvp[,1]))
plot.i5.xmax <- max(c(utc.var.nom.marg.wv1[,1], utc.var.nom.marg.wv2[,1], utc.var.nom.marg.wv3[,1], 
                      utc.var.nom.marg.wv4[,1], utc.var.nom.marg.wv5[,1], utc.var.nom.marg.wv6[,1], 
                      utc.var.nom.marg.bsvp[,1]))
plot.i5.ymin <- min(c(utc.var.nom.marg.wv1[,2], utc.var.nom.marg.wv2[,2], utc.var.nom.marg.wv3[,2], 
                      utc.var.nom.marg.wv4[,2], utc.var.nom.marg.wv5[,2], utc.var.nom.marg.wv6[,2], 
                      utc.var.nom.marg.bsvp[,2]))
plot.i5.ymax <- max(c(utc.var.nom.marg.wv1[,2], utc.var.nom.marg.wv2[,2], utc.var.nom.marg.wv3[,2], 
                      utc.var.nom.marg.wv4[,2], utc.var.nom.marg.wv5[,2], utc.var.nom.marg.wv6[,2], 
                      utc.var.nom.marg.bsvp[,2]))
plot(utc.var.nom.marg.wv1,xlim=c(plot.i5.xmin,plot.i5.xmax),ylim=c(plot.i5.ymin,plot.i5.ymax),
     col="red",t="l",lty=1,xlab=expression(sigma^2),ylab="")
lines(utc.var.nom.marg.wv2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.wv3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.wv4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.wv5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.wv6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.var.nom.marg.bsvp,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ vague prior","2% data w/ vague prior",
                                "5% data w/ vague prior","10% data w/ vague prior",
                                "20% data w/ vague prior", "50% data w/ vague prior",
                                "100% data w/ vague prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)

plot.i6.xmin <- min(c(utc.range.nom.marg.wv1[,1], utc.range.nom.marg.wv2[,1], utc.range.nom.marg.wv3[,1], 
                      utc.range.nom.marg.wv4[,1], utc.range.nom.marg.wv5[,1], utc.range.nom.marg.wv6[,1], 
                      utc.range.nom.marg.bsvp[,1]))
plot.i6.xmax <- max(c(utc.range.nom.marg.wv1[,1], utc.range.nom.marg.wv2[,1], utc.range.nom.marg.wv3[,1], 
                      utc.range.nom.marg.wv4[,1], utc.range.nom.marg.wv5[,1], utc.range.nom.marg.wv6[,1], 
                      utc.range.nom.marg.bsvp[,1]))
plot.i6.ymin <- min(c(utc.range.nom.marg.wv1[,2], utc.range.nom.marg.wv2[,2], utc.range.nom.marg.wv3[,2], 
                      utc.range.nom.marg.wv4[,2], utc.range.nom.marg.wv5[,2], utc.range.nom.marg.wv6[,2], 
                      utc.range.nom.marg.bsvp[,2]))
plot.i6.ymax <- max(c(utc.range.nom.marg.wv1[,2], utc.range.nom.marg.wv2[,2], utc.range.nom.marg.wv3[,2], 
                      utc.range.nom.marg.wv4[,2], utc.range.nom.marg.wv5[,2], utc.range.nom.marg.wv6[,2], 
                      utc.range.nom.marg.bsvp[,2]))
plot(utc.range.nom.marg.wv1,xlim=c(plot.i6.xmin,plot.i6.xmax),ylim=c(plot.i6.ymin,plot.i6.ymax),
     col="red",t="l",lty=1,xlab="r",ylab="")
lines(utc.range.nom.marg.wv2,col="orange",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.wv3,col="green",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.wv4,col="cyan",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.wv5,col="blue",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.wv6,col="purple",t="l",lty=1,xlab="",ylab="")
lines(utc.range.nom.marg.bsvp,col="black",t="l",lty=1,xlab="",ylab="")
legend(x = "topright", legend=c("1% data w/ vague prior","2% data w/ vague prior",
                                "5% data w/ vague prior","10% data w/ vague prior",
                                "20% data w/ vague prior", "50% data w/ vague prior",
                                "100% data w/ vague prior"), 
       col=c("red", "orange", "green", "cyan", "blue", "purple", "black"), lwd=c(1,1,1,1,1,1,1), lty=c(1,1,1,1,1,1,1), cex=0.8)


####################################
# show predictive performances with pcc and rmse (fake)
####################################
percentage_study_data <- read_excel("~/percentage_study_data_fake.xlsx")
percentage_study_data$percentage <- as.factor(percentage_study_data$percentage)

pcc.percentage.un <- subset(percentage_study_data, design=="un", select=c(percentage,pcc))
pcc.percentage.up <- subset(percentage_study_data, design=="up", select=c(percentage,pcc))
pcc.percentage.uv <- subset(percentage_study_data, design=="uv", select=c(percentage,pcc))
pcc.percentage.wn <- subset(percentage_study_data, design=="wn", select=c(percentage,pcc))
pcc.percentage.wp <- subset(percentage_study_data, design=="wp", select=c(percentage,pcc))
pcc.percentage.wv <- subset(percentage_study_data, design=="wv", select=c(percentage,pcc))

rmse.percentage.un <- subset(percentage_study_data, design=="un", select=c(percentage,rmse))
rmse.percentage.up <- subset(percentage_study_data, design=="up", select=c(percentage,rmse))
rmse.percentage.uv <- subset(percentage_study_data, design=="uv", select=c(percentage,rmse))
rmse.percentage.wn <- subset(percentage_study_data, design=="wn", select=c(percentage,rmse))
rmse.percentage.wp <- subset(percentage_study_data, design=="wp", select=c(percentage,rmse))
rmse.percentage.wv <- subset(percentage_study_data, design=="wv", select=c(percentage,rmse))

par(mfrow=c(1,1), mar = c(5, 5, 2, 2))
plot(seq_along(pcc.percentage.un$percentage), pcc.percentage.un$pcc, 
     ylim=c(0.485,0.905),xaxt='n',col="red",t="o",ylab="PCC", xlab="Percentage data acquired by survey", lty=1, pch=1)
axis(1, at=pcc.percentage.un$percentage, labels=c("1%", "2%", "5%","10%","20%","50%","100%"))
lines(pcc.percentage.up,t="o",xlab="",ylab="",col="blue", lty=1, pch=1)
lines(pcc.percentage.uv,t="o",xlab="",ylab="",col="green", lty=1, pch=1)
lines(pcc.percentage.wn,t="o",xlab="",ylab="",col="red", lty=2, pch=1)
lines(pcc.percentage.wp,t="o",xlab="",ylab="",col="blue", lty=2, pch=1)
lines(pcc.percentage.wv,t="o",xlab="",ylab="",col="green", lty=2, pch=1)
abline(h=correl.15as.log, lty=3,col="black")
legend(x = "bottomright", legend=c("worst n. prior + unweighted sample", "worst p. prior + unweighted sample", 
                                   "vague prior + unweighted sample", "worst n. prior + weighted sample",
                                   "worst p. prior + weighted sample", "vague prior + weighted sample", 
                                   "worst data w/ vague prior"), 
       col=c("red","blue", "green", "red","blue", "green","black"), 
       lty=c(1,1,1,2,2,2,3), pch=c(1,1,1,1,1,1,NA), lwd=c(1,1,1,1,1,1,1), cex=0.8)

plot(seq_along(rmse.percentage.un$percentage), rmse.percentage.un$rmse, 
     ylim=c(121.090,257.520),xaxt='n',col="red",t="o",ylab="RMSE", xlab="Percentage data acquired by survey", lty=1, pch=1)
axis(1, at=rmse.percentage.un$percentage, labels=c("1%", "2%", "5%","10%","20%","50%","100%"))
lines(rmse.percentage.up,t="o",xlab="",ylab="",col="blue", lty=1, pch=1)
lines(rmse.percentage.uv,t="o",xlab="",ylab="",col="green", lty=1, pch=1)
lines(rmse.percentage.wn,t="o",xlab="",ylab="",col="red", lty=2, pch=1)
lines(rmse.percentage.wp,t="o",xlab="",ylab="",col="blue", lty=2, pch=1)
lines(rmse.percentage.wv,t="o",xlab="",ylab="",col="green", lty=2, pch=1)
abline(h=RMSE.15as.log, lty=3,col="black")
legend(x = "topright", legend=c("worst n. prior + unweighted sample", "worst p. prior + unweighted sample", 
                                "vague prior + unweighted sample", "worst n. prior + weighted sample",
                                "worst p. prior + weighted sample", "vague prior + weighted sample", 
                                "worst data w/ vague prior"), 
       col=c("red","blue", "green", "red","blue", "green","black"), 
       lty=c(1,1,1,2,2,2,3), pch=c(1,1,1,1,1,1,NA), lwd=c(1,1,1,1,1,1,1), cex=0.8)


####################################
# show predictive performances with pcc and rmse (real)
####################################
percentage_study_data <- read_excel("~/percentage_study_data_real.xlsx")
percentage_study_data$percentage <- as.factor(percentage_study_data$percentage)

pcc.percentage.un <- subset(percentage_study_data, design=="un", select=c(percentage,pcc))
pcc.percentage.up <- subset(percentage_study_data, design=="up", select=c(percentage,pcc))
pcc.percentage.uv <- subset(percentage_study_data, design=="uv", select=c(percentage,pcc))
pcc.percentage.ub <- subset(percentage_study_data, design=="ub", select=c(percentage,pcc))
pcc.percentage.wn <- subset(percentage_study_data, design=="wn", select=c(percentage,pcc))
pcc.percentage.wp <- subset(percentage_study_data, design=="wp", select=c(percentage,pcc))
pcc.percentage.wv <- subset(percentage_study_data, design=="wv", select=c(percentage,pcc))
pcc.percentage.wb <- subset(percentage_study_data, design=="wb", select=c(percentage,pcc))

rmse.percentage.un <- subset(percentage_study_data, design=="un", select=c(percentage,rmse))
rmse.percentage.up <- subset(percentage_study_data, design=="up", select=c(percentage,rmse))
rmse.percentage.uv <- subset(percentage_study_data, design=="uv", select=c(percentage,rmse))
rmse.percentage.ub <- subset(percentage_study_data, design=="ub", select=c(percentage,rmse))
rmse.percentage.wn <- subset(percentage_study_data, design=="wn", select=c(percentage,rmse))
rmse.percentage.wp <- subset(percentage_study_data, design=="wp", select=c(percentage,rmse))
rmse.percentage.wv <- subset(percentage_study_data, design=="wv", select=c(percentage,rmse))
rmse.percentage.wb <- subset(percentage_study_data, design=="wb", select=c(percentage,rmse))

par(mfrow=c(1,1), mar = c(5, 5, 2, 2))
plot(seq_along(pcc.percentage.un$percentage), pcc.percentage.un$pcc, 
     ylim=c(0.055,0.727),xaxt='n',col="red",t="o",ylab="PCC", xlab="Percentage data acquired by survey", lty=1, pch=1)
axis(1, at=pcc.percentage.un$percentage, labels=c("1%", "2%", "5%","10%","20%","50%","100%"))
lines(pcc.percentage.up,t="o",xlab="",ylab="",col="blue", lty=1, pch=1)
lines(pcc.percentage.uv,t="o",xlab="",ylab="",col="green", lty=1, pch=1)
lines(pcc.percentage.ub,t="o",xlab="",ylab="",col="black", lty=1, pch=1)
lines(pcc.percentage.wn,t="o",xlab="",ylab="",col="red", lty=2, pch=1)
lines(pcc.percentage.wp,t="o",xlab="",ylab="",col="blue", lty=2, pch=1)
lines(pcc.percentage.wv,t="o",xlab="",ylab="",col="green", lty=2, pch=1)
lines(pcc.percentage.wb,t="o",xlab="",ylab="",col="black", lty=2, pch=1)
abline(h=correl.15as.log, lty=3,col="black")
legend(x = "bottomright", legend=c("worst n. prior + unweighted sample", "worst p. prior + unweighted sample", 
                                   "vague prior + unweighted sample", "bottom-up + unweighted sample", 
                                   "worst n. prior + weighted sample", "worst p. prior + weighted sample", 
                                   "vague prior + weighted sample", "bottom-up + weighted sample", 
                                   "worst data w/ vague prior"), 
       col=c("red","blue", "green", "black", "red","blue", "green","black", "black"), 
       lty=c(1,1,1,1,2,2,2,2,3), pch=c(1,1,1,1,1,1,1,1,NA), lwd=c(1,1,1,1,1,1,1,1,1), cex=0.8)

plot(seq_along(rmse.percentage.un$percentage), rmse.percentage.un$rmse, 
     ylim=c(366.631,600),xaxt='n',col="red",t="o",ylab="RMSE", xlab="Percentage data acquired by survey", lty=1, pch=1)
axis(1, at=rmse.percentage.un$percentage, labels=c("1%", "2%", "5%","10%","20%","50%","100%"))
lines(rmse.percentage.up,t="o",xlab="",ylab="",col="blue", lty=1, pch=1)
lines(rmse.percentage.uv,t="o",xlab="",ylab="",col="green", lty=1, pch=1)
lines(rmse.percentage.ub,t="o",xlab="",ylab="",col="black", lty=1, pch=1)
lines(rmse.percentage.wn,t="o",xlab="",ylab="",col="red", lty=2, pch=1)
lines(rmse.percentage.wp,t="o",xlab="",ylab="",col="blue", lty=2, pch=1)
lines(rmse.percentage.wv,t="o",xlab="",ylab="",col="green", lty=2, pch=1)
lines(rmse.percentage.wb,t="o",xlab="",ylab="",col="black", lty=2, pch=1)
abline(h=RMSE.15as.log, lty=3,col="black")
legend(x = "topright", legend=c("worst n. prior + unweighted sample", "worst p. prior + unweighted sample", 
                                "vague prior + unweighted sample", "bottom-up + unweighted sample", 
                                "worst n. prior + weighted sample", "worst p. prior + weighted sample", 
                                "vague prior + weighted sample", "bottom-up + weighted sample", 
                                "worst data w/ vague prior"), 
       col=c("red","blue", "green", "black", "red","blue", "green", "black","black"), 
       lty=c(1,1,1,1,2,2,2,2,3), pch=c(1,1,1,1,1,1,1,1,NA), lwd=c(1,1,1,1,1,1,1,1,1), cex=0.8)

####################################
# min and max for plotting posterior mean and sd
####################################
max(log(as.numeric(stacked_utc_pred_15as_fit$pop_utc_15_15as)+1),
    log(as.numeric(stacked_utc_30as_fit$pop_utc_15_30as)+1),
    log(utc.pred.response.mean.log.mod/4+1),
    log(utc.pred.response.mean.bsvp.mod/4+1), 
    log(utc.pred.response.mean.bswp.n.mod/4+1), 
    log(utc.pred.response.mean.bswp.p.mod/4+1),
    log(utc.pred.response.mean.un1.mod/4+1),
    log(utc.pred.response.mean.un2.mod/4+1),
    log(utc.pred.response.mean.un3.mod/4+1),
    log(utc.pred.response.mean.un4.mod/4+1),
    log(utc.pred.response.mean.un5.mod/4+1),
    log(utc.pred.response.mean.un6.mod/4+1),
    log(utc.pred.response.mean.up1.mod/4+1),
    log(utc.pred.response.mean.up2.mod/4+1),
    log(utc.pred.response.mean.up3.mod/4+1),
    log(utc.pred.response.mean.up4.mod/4+1),
    log(utc.pred.response.mean.up5.mod/4+1),
    log(utc.pred.response.mean.up6.mod/4+1),
    log(utc.pred.response.mean.uv1.mod/4+1),
    log(utc.pred.response.mean.uv2.mod/4+1),
    log(utc.pred.response.mean.uv3.mod/4+1),
    log(utc.pred.response.mean.uv4.mod/4+1),
    log(utc.pred.response.mean.uv5.mod/4+1),
    log(utc.pred.response.mean.uv6.mod/4+1),
    log(utc.pred.response.mean.wn1.mod/4+1),
    log(utc.pred.response.mean.wn2.mod/4+1),
    log(utc.pred.response.mean.wn3.mod/4+1),
    log(utc.pred.response.mean.wn4.mod/4+1),
    log(utc.pred.response.mean.wn5.mod/4+1),
    log(utc.pred.response.mean.wn6.mod/4+1),
    log(utc.pred.response.mean.wp1.mod/4+1),
    log(utc.pred.response.mean.wp2.mod/4+1),
    log(utc.pred.response.mean.wp3.mod/4+1),
    log(utc.pred.response.mean.wp4.mod/4+1),
    log(utc.pred.response.mean.wp5.mod/4+1),
    log(utc.pred.response.mean.wp6.mod/4+1),
    log(utc.pred.response.mean.wv1.mod/4+1),
    log(utc.pred.response.mean.wv2.mod/4+1),
    log(utc.pred.response.mean.wv3.mod/4+1),
    log(utc.pred.response.mean.wv4.mod/4+1),
    log(utc.pred.response.mean.wv5.mod/4+1),
    log(utc.pred.response.mean.wv6.mod/4+1))

max(utc.pred.response.sd.log,
    utc.pred.response.sd.bsvp, 
    utc.pred.response.sd.bswp.n, 
    utc.pred.response.sd.bswp.p,
    utc.pred.response.sd.un1,
    utc.pred.response.sd.un2,
    utc.pred.response.sd.un3,
    utc.pred.response.sd.un4,
    utc.pred.response.sd.un5,
    utc.pred.response.sd.un6,
    utc.pred.response.sd.up1,
    utc.pred.response.sd.up2,
    utc.pred.response.sd.up3,
    utc.pred.response.sd.up4,
    utc.pred.response.sd.up5,
    utc.pred.response.sd.up6,
    utc.pred.response.sd.uv1,
    utc.pred.response.sd.uv2,
    utc.pred.response.sd.uv3,
    utc.pred.response.sd.uv4,
    utc.pred.response.sd.uv5,
    utc.pred.response.sd.uv6,
    utc.pred.response.sd.wn1,
    utc.pred.response.sd.wn2,
    utc.pred.response.sd.wn3,
    utc.pred.response.sd.wn4,
    utc.pred.response.sd.wn5,
    utc.pred.response.sd.wn6,
    utc.pred.response.sd.wp1,
    utc.pred.response.sd.wp2,
    utc.pred.response.sd.wp3,
    utc.pred.response.sd.wp4,
    utc.pred.response.sd.wp5,
    utc.pred.response.sd.wp6,
    utc.pred.response.sd.wv1,
    utc.pred.response.sd.wv2,
    utc.pred.response.sd.wv3,
    utc.pred.response.sd.wv4,
    utc.pred.response.sd.wv5,
    utc.pred.response.sd.wv6)

max(log(exp(utc.pred.response.mean.bt.best)+1),
    log(exp(utc.pred.response.mean.bt.u1)+1),
    log(exp(utc.pred.response.mean.bt.u2)+1),
    log(exp(utc.pred.response.mean.bt.u3)+1),
    log(exp(utc.pred.response.mean.bt.u4)+1),
    log(exp(utc.pred.response.mean.bt.u5)+1),
    log(exp(utc.pred.response.mean.bt.u6)+1),
    log(exp(utc.pred.response.mean.bt.w1)+1),
    log(exp(utc.pred.response.mean.bt.w2)+1),
    log(exp(utc.pred.response.mean.bt.w3)+1),
    log(exp(utc.pred.response.mean.bt.w4)+1),
    log(exp(utc.pred.response.mean.bt.w5)+1),
    log(exp(utc.pred.response.mean.bt.w6)+1))

####################################
### exploratory data analysis for nairobi: assess assumptions
### a temporary extension
####################################
# data for nab
pop_nab_01_30as <- raster("~/data/wppop_nab_2001_30arcsec.tif")
pop_nab_04_30as <- raster("~/data/wppop_nab_2004_30arcsec.tif")
pop_nab_07_30as <- raster("~/data/wppop_nab_2007_30arcsec.tif")
pop_nab_10_30as <- raster("~/data/wppop_nab_2010_30arcsec.tif")
pop_nab_13_30as <- raster("~/data/wppop_nab_2013_30arcsec.tif")
pop_nab_16_30as <- raster("~/data/wppop_nab_2016_30arcsec.tif")

pop_nab_02_30as <- raster("~/data/wppop_nab_2002_30arcsec.tif")
pop_nab_05_30as <- raster("~/data/wppop_nab_2005_30arcsec.tif")
pop_nab_08_30as <- raster("~/data/wppop_nab_2008_30arcsec.tif")
pop_nab_11_30as <- raster("~/data/wppop_nab_2011_30arcsec.tif")
pop_nab_14_30as <- raster("~/data/wppop_nab_2014_30arcsec.tif")
pop_nab_17_30as <- raster("~/data/wppop_nab_2017_30arcsec.tif")

pop_nab_03_30as <- raster("~/data/wppop_nab_2003_30arcsec.tif")
pop_nab_06_30as <- raster("~/data/wppop_nab_2006_30arcsec.tif")
pop_nab_09_30as <- raster("~/data/wppop_nab_2009_30arcsec.tif")
pop_nab_12_30as <- raster("~/data/wppop_nab_2012_30arcsec.tif")
pop_nab_15_30as <- raster("~/data/wppop_nab_2015_30arcsec.tif")
pop_nab_18_30as <- raster("~/data/wppop_nab_2018_30arcsec.tif")

stacked_nab_eda <- stack(pop_nab_01_30as, pop_nab_04_30as, pop_nab_07_30as, pop_nab_10_30as, pop_nab_13_30as, pop_nab_16_30as,
                         pop_nab_02_30as, pop_nab_05_30as, pop_nab_08_30as, pop_nab_11_30as, pop_nab_14_30as, pop_nab_17_30as, 
                         pop_nab_03_30as, pop_nab_06_30as, pop_nab_09_30as, pop_nab_12_30as, pop_nab_15_30as, pop_nab_18_30as)

stacked_nab_eda_df <- as.data.frame(stacked_nab_eda, xy = TRUE)

colnames(stacked_nab_eda_df) <- c("x", "y", 
                                  "pop_nab_01_30as", "pop_nab_04_30as", "pop_nab_07_30as", "pop_nab_10_30as", "pop_nab_13_30as", "pop_nab_16_30as",
                                  "pop_nab_02_30as", "pop_nab_05_30as", "pop_nab_08_30as", "pop_nab_11_30as", "pop_nab_14_30as", "pop_nab_17_30as",  
                                  "pop_nab_03_30as", "pop_nab_06_30as", "pop_nab_09_30as", "pop_nab_12_30as", "pop_nab_15_30as", "pop_nab_18_30as")
stacked_nab_eda_df[is.na(stacked_nab_eda_df)] <- 0

nab_outline <- readOGR(dsn="~/data",layer="nairobi_outline_2007")
crs.geo <- CRS("+proj=longlat +datum=WGS84 +no_defs")
nab_outline  <- spTransform(nab_outline, crs.geo)
nab.ov.30as <- over(SpatialPoints(stacked_nab_eda_df[,1:2], proj4string=crs.geo), nab_outline)
stacked_nab_eda_df <- cbind(stacked_nab_eda_df, nab.ov.30as)
stacked_nab_eda_df <- stacked_nab_eda_df[!is.na(stacked_nab_eda_df$FID),]

nab.loc_eda <- cbind(stacked_nab_eda_df$x, stacked_nab_eda_df$y)

nab.01.sp <- SpatialPointsDataFrame(SpatialPoints(nab.loc_eda), data = data.frame(counts = log(stacked_nab_eda_df$pop_nab_01_30as+1)))
nab.04.sp <- SpatialPointsDataFrame(SpatialPoints(nab.loc_eda), data = data.frame(counts = log(stacked_nab_eda_df$pop_nab_04_30as+1)))
nab.07.sp <- SpatialPointsDataFrame(SpatialPoints(nab.loc_eda), data = data.frame(counts = log(stacked_nab_eda_df$pop_nab_07_30as+1)))
nab.10.sp <- SpatialPointsDataFrame(SpatialPoints(nab.loc_eda), data = data.frame(counts = log(stacked_nab_eda_df$pop_nab_10_30as+1)))
nab.13.sp <- SpatialPointsDataFrame(SpatialPoints(nab.loc_eda), data = data.frame(counts = log(stacked_nab_eda_df$pop_nab_13_30as+1)))
nab.16.sp <- SpatialPointsDataFrame(SpatialPoints(nab.loc_eda), data = data.frame(counts = log(stacked_nab_eda_df$pop_nab_16_30as+1)))

nab.01.g <- gstat(id = "counts", formula = counts ~ 1, data = nab.01.sp)
nab.01.evgm <- variogram(nab.01.g, cutoff = 0.12, width = 0.003)
nab.01.revgm <- variogram(nab.01.g, cutoff = 0.12, width = 0.003, cressie = TRUE)
nab.01.aevgm <- variogram(nab.01.g, cutoff = 0.12, width = 0.003, alpha = c(0, 45, 90, 135))

nab.04.g <- gstat(id = "counts", formula = counts ~ 1, data = nab.04.sp)
nab.04.evgm <- variogram(nab.04.g, cutoff = 0.12, width = 0.003)
nab.04.revgm <- variogram(nab.04.g, cutoff = 0.12, width = 0.003, cressie = TRUE)
nab.04.aevgm <- variogram(nab.04.g, cutoff = 0.12, width = 0.003, alpha = c(0, 45, 90, 135))

nab.07.g <- gstat(id = "counts", formula = counts ~ 1, data = nab.07.sp)
nab.07.evgm <- variogram(nab.07.g, cutoff = 0.12, width = 0.003)
nab.07.revgm <- variogram(nab.07.g, cutoff = 0.12, width = 0.003, cressie = TRUE)
nab.07.aevgm <- variogram(nab.07.g, cutoff = 0.12, width = 0.003, alpha = c(0, 45, 90, 135))

nab.10.g <- gstat(id = "counts", formula = counts ~ 1, data = nab.10.sp)
nab.10.evgm <- variogram(nab.10.g, cutoff = 0.12, width = 0.003)
nab.10.revgm <- variogram(nab.10.g, cutoff = 0.12, width = 0.003, cressie = TRUE)
nab.10.aevgm <- variogram(nab.10.g, cutoff = 0.12, width = 0.003, alpha = c(0, 45, 90, 135))

nab.13.g <- gstat(id = "counts", formula = counts ~ 1, data = nab.13.sp)
nab.13.evgm <- variogram(nab.13.g, cutoff = 0.12, width = 0.003)
nab.13.revgm <- variogram(nab.13.g, cutoff = 0.12, width = 0.003, cressie = TRUE)
nab.13.aevgm <- variogram(nab.13.g, cutoff = 0.12, width = 0.003, alpha = c(0, 45, 90, 135))

nab.16.g <- gstat(id = "counts", formula = counts ~ 1, data = nab.16.sp)
nab.16.evgm <- variogram(nab.16.g, cutoff = 0.12, width = 0.003)
nab.16.revgm <- variogram(nab.16.g, cutoff = 0.12, width = 0.003, cressie = TRUE)
nab.16.aevgm <- variogram(nab.16.g, cutoff = 0.12, width = 0.003, alpha = c(0, 45, 90, 135))

# plot semivariogram for nab
# check assumption of second-order stationary
par(mfrow=c(1,1), mar = c(5, 5, 2, 2))
plot(cbind(nab.01.evgm$dist, nab.01.evgm$gamma), 
     xlim=c(0,0.125),ylim=c(0,6),col="red",t="o",ylab=expression(hat(gamma)), xlab="Distance (degree)", lty=1, pch=20)
lines(cbind(nab.01.revgm$dist, nab.01.revgm$gamma),t="o",xlab="",ylab="",col="red", lty=3, pch=20)
lines(cbind(nab.04.evgm$dist, nab.04.evgm$gamma),t="o",xlab="",ylab="",col="orange", lty=1, pch=20)
lines(cbind(nab.04.revgm$dist, nab.04.revgm$gamma),t="o",xlab="",ylab="",col="orange", lty=3, pch=20)
lines(cbind(nab.07.evgm$dist, nab.07.evgm$gamma),t="o",xlab="",ylab="",col="green", lty=1, pch=20)
lines(cbind(nab.07.revgm$dist, nab.07.revgm$gamma),t="o",xlab="",ylab="",col="green", lty=3, pch=20)
lines(cbind(nab.10.evgm$dist, nab.10.evgm$gamma),t="o",xlab="",ylab="",col="blue", lty=1, pch=20)
lines(cbind(nab.10.revgm$dist, nab.10.revgm$gamma),t="o",xlab="",ylab="",col="blue", lty=3, pch=20)
lines(cbind(nab.13.evgm$dist, nab.13.evgm$gamma),t="o",xlab="",ylab="",col="cyan", lty=1, pch=20)
lines(cbind(nab.13.revgm$dist, nab.13.revgm$gamma),t="o",xlab="",ylab="",col="cyan", lty=3, pch=20)
lines(cbind(nab.16.evgm$dist, nab.16.evgm$gamma),t="o",xlab="",ylab="",col="purple", lty=1, pch=20)
lines(cbind(nab.16.revgm$dist, nab.16.revgm$gamma),t="o",xlab="",ylab="",col="purple", lty=3, pch=20)
legend(x = "bottomright", legend=c("2001, classic", "2001, robust", "2004, classic","2004, robust",
                                   "2007, classic", "2007, robust", "2010, classic","2010, robust",
                                   "2013, classic", "2013, robust", "2016, classic","2016, robust"), 
       col=c("red", "red", "orange","orange","green","green","blue","blue","cyan","cyan","purple","purple"), 
       lty=c(1,3,1,3,1,3,1,3,1,3,1,3), pch=c(20,20,20,20,20,20,20,20,20,20,20,20), lwd=c(1,1), cex=0.8)

# check assumption of isotropic
nab.01.aevgm.0 <- subset(nab.01.aevgm, dir.hor==0, select=c(dist,gamma))
nab.01.aevgm.45 <- subset(nab.01.aevgm, dir.hor==45, select=c(dist,gamma))
nab.01.aevgm.90 <- subset(nab.01.aevgm, dir.hor==90, select=c(dist,gamma))
nab.01.aevgm.135 <- subset(nab.01.aevgm, dir.hor==135, select=c(dist,gamma))

nab.04.aevgm.0 <- subset(nab.04.aevgm, dir.hor==0, select=c(dist,gamma))
nab.04.aevgm.45 <- subset(nab.04.aevgm, dir.hor==45, select=c(dist,gamma))
nab.04.aevgm.90 <- subset(nab.04.aevgm, dir.hor==90, select=c(dist,gamma))
nab.04.aevgm.135 <- subset(nab.04.aevgm, dir.hor==135, select=c(dist,gamma))

nab.07.aevgm.0 <- subset(nab.07.aevgm, dir.hor==0, select=c(dist,gamma))
nab.07.aevgm.45 <- subset(nab.07.aevgm, dir.hor==45, select=c(dist,gamma))
nab.07.aevgm.90 <- subset(nab.07.aevgm, dir.hor==90, select=c(dist,gamma))
nab.07.aevgm.135 <- subset(nab.07.aevgm, dir.hor==135, select=c(dist,gamma))

nab.10.aevgm.0 <- subset(nab.10.aevgm, dir.hor==0, select=c(dist,gamma))
nab.10.aevgm.45 <- subset(nab.10.aevgm, dir.hor==45, select=c(dist,gamma))
nab.10.aevgm.90 <- subset(nab.10.aevgm, dir.hor==90, select=c(dist,gamma))
nab.10.aevgm.135 <- subset(nab.10.aevgm, dir.hor==135, select=c(dist,gamma))

nab.13.aevgm.0 <- subset(nab.13.aevgm, dir.hor==0, select=c(dist,gamma))
nab.13.aevgm.45 <- subset(nab.13.aevgm, dir.hor==45, select=c(dist,gamma))
nab.13.aevgm.90 <- subset(nab.13.aevgm, dir.hor==90, select=c(dist,gamma))
nab.13.aevgm.135 <- subset(nab.13.aevgm, dir.hor==135, select=c(dist,gamma))

nab.16.aevgm.0 <- subset(nab.16.aevgm, dir.hor==0, select=c(dist,gamma))
nab.16.aevgm.45 <- subset(nab.16.aevgm, dir.hor==45, select=c(dist,gamma))
nab.16.aevgm.90 <- subset(nab.16.aevgm, dir.hor==90, select=c(dist,gamma))
nab.16.aevgm.135 <- subset(nab.16.aevgm, dir.hor==135, select=c(dist,gamma))

plot(cbind(nab.01.aevgm.0$dist, nab.01.aevgm.0$gamma), 
     xlim=c(0,0.125),ylim=c(0,9),col="red",t="o",ylab=expression(hat(gamma)), xlab="Distance (degree)", lty=1, pch=15)
lines(cbind(nab.01.aevgm.45$dist, nab.01.aevgm.45$gamma),t="o",xlab="",ylab="",col="red", lty=1, pch=16)
lines(cbind(nab.01.aevgm.90$dist, nab.01.aevgm.90$gamma),t="o",xlab="",ylab="",col="red", lty=1, pch=17)
lines(cbind(nab.01.aevgm.135$dist, nab.01.aevgm.135$gamma),t="o",xlab="",ylab="",col="red", lty=1, pch=18)
lines(cbind(nab.04.aevgm.0$dist, nab.04.aevgm.0$gamma),t="o",xlab="",ylab="",col="orange", lty=1, pch=15)
lines(cbind(nab.04.aevgm.45$dist, nab.04.aevgm.45$gamma),t="o",xlab="",ylab="",col="orange", lty=1, pch=16)
lines(cbind(nab.04.aevgm.90$dist, nab.04.aevgm.90$gamma),t="o",xlab="",ylab="",col="orange", lty=1, pch=17)
lines(cbind(nab.04.aevgm.135$dist, nab.04.aevgm.135$gamma),t="o",xlab="",ylab="",col="orange", lty=1, pch=18)
lines(cbind(nab.07.aevgm.0$dist, nab.07.aevgm.0$gamma),t="o",xlab="",ylab="",col="green", lty=1, pch=15)
lines(cbind(nab.07.aevgm.45$dist, nab.07.aevgm.45$gamma),t="o",xlab="",ylab="",col="green", lty=1, pch=16)
lines(cbind(nab.07.aevgm.90$dist, nab.07.aevgm.90$gamma),t="o",xlab="",ylab="",col="green", lty=1, pch=17)
lines(cbind(nab.07.aevgm.135$dist, nab.07.aevgm.135$gamma),t="o",xlab="",ylab="",col="green", lty=1, pch=18)
lines(cbind(nab.10.aevgm.0$dist, nab.10.aevgm.0$gamma),t="o",xlab="",ylab="",col="blue", lty=1, pch=15)
lines(cbind(nab.10.aevgm.45$dist, nab.10.aevgm.45$gamma),t="o",xlab="",ylab="",col="blue", lty=1, pch=16)
lines(cbind(nab.10.aevgm.90$dist, nab.10.aevgm.90$gamma),t="o",xlab="",ylab="",col="blue", lty=1, pch=17)
lines(cbind(nab.10.aevgm.135$dist, nab.10.aevgm.135$gamma),t="o",xlab="",ylab="",col="blue", lty=1, pch=18)
lines(cbind(nab.13.aevgm.0$dist, nab.13.aevgm.0$gamma),t="o",xlab="",ylab="",col="cyan", lty=1, pch=15)
lines(cbind(nab.13.aevgm.45$dist, nab.13.aevgm.45$gamma),t="o",xlab="",ylab="",col="cyan", lty=1, pch=16)
lines(cbind(nab.13.aevgm.90$dist, nab.13.aevgm.90$gamma),t="o",xlab="",ylab="",col="cyan", lty=1, pch=17)
lines(cbind(nab.13.aevgm.135$dist, nab.13.aevgm.135$gamma),t="o",xlab="",ylab="",col="cyan", lty=1, pch=18)
lines(cbind(nab.16.aevgm.0$dist, nab.16.aevgm.0$gamma),t="o",xlab="",ylab="",col="purple", lty=1, pch=15)
lines(cbind(nab.16.aevgm.45$dist, nab.16.aevgm.45$gamma),t="o",xlab="",ylab="",col="purple", lty=1, pch=16)
lines(cbind(nab.16.aevgm.90$dist, nab.16.aevgm.90$gamma),t="o",xlab="",ylab="",col="purple", lty=1, pch=17)
lines(cbind(nab.16.aevgm.135$dist, nab.16.aevgm.135$gamma),t="o",xlab="",ylab="",col="purple", lty=1, pch=18)
legend(x = "topleft", legend=c("2001, 0°", "2001, 45°", "2001, 90°","2001, 135°", 
                               "2004, 0°", "2004, 45°", "2004, 90°","2004, 135°",
                               "2007, 0°", "2007, 45°", "2007, 90°","2007, 135°",
                               "2010, 0°", "2010, 45°", "2010, 90°","2010, 135°",
                               "2013, 0°", "2013, 45°", "2013, 90°","2013, 135°",
                               "2016, 0°", "2016, 45°", "2016, 90°","2016, 135°"), 
       col=c("red","red", "red","red","orange","orange","orange","orange",
             "green","green","green","green","blue", "blue", "blue", "blue",
             "cyan","cyan","cyan","cyan","purple","purple","purple","purple"), 
       lty=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), pch=c(15,16,17,18, 15,16,17,18, 15,16,17,18, 
                                                                     15,16,17,18, 15,16,17,18, 15,16,17,18), lwd=c(1,1),ncol=2, cex=0.8)


# boxplot for nab
nab.pop.full.eda <- as.vector(log(as.matrix(stacked_nab_eda_df[,3:20])+1))
nab.year.span <- c("2001-2003", "2004-2006", "2007-2009", "2010-2012", "2013-2015", "2016-2018")
nab.year.eda <- c(rep(nab.year.span, each=dim(stacked_nab_eda_df)[1]), 
                  rep(nab.year.span, each=dim(stacked_nab_eda_df)[1]), 
                  rep(nab.year.span, each=dim(stacked_nab_eda_df)[1]))
k <- 6
nab.group <- c(rep("est group", dim(stacked_nab_eda_df)[1]*k), 
               rep("pred group 1", dim(stacked_nab_eda_df)[1]*k), 
               rep("pred group 2", dim(stacked_nab_eda_df)[1]*k))
nab.eda.data <- data.frame(nab.year.eda, nab.pop.full.eda, nab.group)

# boxplot for nab
boxplot(nab.pop.full.eda~nab.group+
          factor(nab.year.eda, levels=c("2001-2003", "2004-2006", "2007-2009", "2010-2012", "2013-2015", "2016-2018")), 
        data=nab.eda.data, col=c("red","green","blue"), at=c(1:3,5:7,9:11,13:15,17:19,21:23), xaxt = 'n',
        ylab="Population counts on the transformed scale", xlab="Years")
axis(side=1,at=c(2,6,10,14,18,22),
     labels=c("2001-2003", "2004-2006", "2007-2009", "2010-2012", "2013-2015", "2016-2018"))
legend("topleft", c("Estimation 0","Prediction 1","Prediction 2"), fill=c("red", "green","blue"), horiz=TRUE, cex=0.8)


####################################
###  import data of narirobi for estimation group
####################################
# import nightlight map
ntl_nab_01_30as <- raster("~/data/harmntl_nab_2001_30arcsec.tif")
ntl_nab_04_30as <- raster("~/data/harmntl_nab_2004_30arcsec.tif")
ntl_nab_07_30as <- raster("~/data/harmntl_nab_2007_30arcsec.tif")
ntl_nab_10_30as <- raster("~/data/harmntl_nab_2010_30arcsec.tif")
ntl_nab_13_30as <- raster("~/data/harmntl_nab_2013_30arcsec.tif")
ntl_nab_16_30as <- raster("~/data/harmntl_nab_2016_30arcsec.tif")

# import dem map
slope_nab_3as <- raster("~/data/slope_nab_3arcsec.tif")

# resample nightlight
ntl_nab_01_resampled <- resample(ntl_nab_01_30as, pop_nab_01_30as, method = 'ngb')
ntl_nab_04_resampled <- resample(ntl_nab_04_30as, pop_nab_04_30as, method = 'ngb') 
ntl_nab_07_resampled <- resample(ntl_nab_07_30as, pop_nab_07_30as, method = 'ngb') 
ntl_nab_10_resampled <- resample(ntl_nab_10_30as, pop_nab_10_30as, method = 'ngb') 
ntl_nab_13_resampled <- resample(ntl_nab_13_30as, pop_nab_13_30as, method = 'ngb') 
ntl_nab_16_resampled <- resample(ntl_nab_16_30as, pop_nab_16_30as, method = 'ngb') 

# resample dem
slope_nab_resampled <- resample(slope_nab_3as, pop_nab_16_30as, method = 'ngb') 

stacked_nab <- stack(pop_nab_01_30as, pop_nab_04_30as, pop_nab_07_30as, pop_nab_10_30as, pop_nab_13_30as, pop_nab_16_30as, 
                     ntl_nab_01_resampled, ntl_nab_04_resampled, ntl_nab_07_resampled, ntl_nab_10_resampled, ntl_nab_13_resampled, ntl_nab_16_resampled, 
                     slope_nab_resampled)

stacked_nab_df <- as.data.frame(stacked_nab, xy = TRUE)

colnames(stacked_nab_df) <- c("x", "y", "pop_nab_01_30as", "pop_nab_04_30as", "pop_nab_07_30as", "pop_nab_10_30as", "pop_nab_13_30as", "pop_nab_16_30as",  
                              "ntl_nab_01_30as", "ntl_nab_04_30as", "ntl_nab_07_30as", "ntl_nab_10_30as", "ntl_nab_13_30as", "ntl_nab_16_30as",  
                              "slope_nab_30as")
stacked_nab_df[is.na(stacked_nab_df)] <- 0

stacked_nab_df <- cbind(stacked_nab_df, nab.ov.30as)
stacked_nab_df <- stacked_nab_df[!is.na(stacked_nab_df$FID),]

####################################
### import data of narirobi for prediction group 1 and 2
####################################
# import ntl grid for nairobi
ntl_nab_02_30as_grid <- raster("~/data/harmntl_nab_2002_30arcsec_grid.tif")
ntl_nab_03_30as_grid <- raster("~/data/harmntl_nab_2003_30arcsec_grid.tif")
ntl_nab_05_30as_grid <- raster("~/data/harmntl_nab_2005_30arcsec_grid.tif")
ntl_nab_06_30as_grid <- raster("~/data/harmntl_nab_2006_30arcsec_grid.tif")
ntl_nab_08_30as_grid <- raster("~/data/harmntl_nab_2008_30arcsec_grid.tif")
ntl_nab_09_30as_grid <- raster("~/data/harmntl_nab_2009_30arcsec_grid.tif")
ntl_nab_11_30as_grid <- raster("~/data/harmntl_nab_2011_30arcsec_grid.tif")
ntl_nab_12_30as_grid <- raster("~/data/harmntl_nab_2012_30arcsec_grid.tif")
ntl_nab_14_30as_grid <- raster("~/data/harmntl_nab_2014_30arcsec_grid.tif")
ntl_nab_15_30as_grid <- raster("~/data/harmntl_nab_2015_30arcsec_grid.tif")
ntl_nab_17_30as_grid <- raster("~/data/harmntl_nab_2017_30arcsec_grid.tif")
ntl_nab_18_30as_grid <- raster("~/data/harmntl_nab_2018_30arcsec_grid.tif")

# import slope grid for nairobi
slope_nab_3as_grid <- raster("~/data/slope_nab_3arcsec_grid.tif")

# import pop grid for nairobi
pop_nab_02_30as_grid <- raster("~/data/wppop_nab_2002_30arcsec_grid.tif")
pop_nab_03_30as_grid <- raster("~/data/wppop_nab_2003_30arcsec_grid.tif")
pop_nab_05_30as_grid <- raster("~/data/wppop_nab_2005_30arcsec_grid.tif")
pop_nab_06_30as_grid <- raster("~/data/wppop_nab_2006_30arcsec_grid.tif")
pop_nab_08_30as_grid <- raster("~/data/wppop_nab_2008_30arcsec_grid.tif")
pop_nab_09_30as_grid <- raster("~/data/wppop_nab_2009_30arcsec_grid.tif")
pop_nab_11_30as_grid <- raster("~/data/wppop_nab_2011_30arcsec_grid.tif")
pop_nab_12_30as_grid <- raster("~/data/wppop_nab_2012_30arcsec_grid.tif")
pop_nab_14_30as_grid <- raster("~/data/wppop_nab_2014_30arcsec_grid.tif")
pop_nab_15_30as_grid <- raster("~/data/wppop_nab_2015_30arcsec_grid.tif")
pop_nab_17_30as_grid <- raster("~/data/wppop_nab_2017_30arcsec_grid.tif")
pop_nab_18_30as_grid <- raster("~/data/wppop_nab_2018_30arcsec_grid.tif")

# resample pop grid
ntl_nab_02_resampled_pred <- resample(ntl_nab_02_30as_grid, pop_nab_02_30as_grid, method = 'ngb')
ntl_nab_03_resampled_pred <- resample(ntl_nab_03_30as_grid, pop_nab_03_30as_grid, method = 'ngb')
ntl_nab_05_resampled_pred <- resample(ntl_nab_05_30as_grid, pop_nab_05_30as_grid, method = 'ngb')
ntl_nab_06_resampled_pred <- resample(ntl_nab_06_30as_grid, pop_nab_06_30as_grid, method = 'ngb')
ntl_nab_08_resampled_pred <- resample(ntl_nab_08_30as_grid, pop_nab_08_30as_grid, method = 'ngb')
ntl_nab_09_resampled_pred <- resample(ntl_nab_09_30as_grid, pop_nab_09_30as_grid, method = 'ngb')
ntl_nab_11_resampled_pred <- resample(ntl_nab_11_30as_grid, pop_nab_11_30as_grid, method = 'ngb')
ntl_nab_12_resampled_pred <- resample(ntl_nab_12_30as_grid, pop_nab_12_30as_grid, method = 'ngb')
ntl_nab_14_resampled_pred <- resample(ntl_nab_14_30as_grid, pop_nab_14_30as_grid, method = 'ngb')
ntl_nab_15_resampled_pred <- resample(ntl_nab_15_30as_grid, pop_nab_15_30as_grid, method = 'ngb')
ntl_nab_17_resampled_pred <- resample(ntl_nab_17_30as_grid, pop_nab_17_30as_grid, method = 'ngb')
ntl_nab_18_resampled_pred <- resample(ntl_nab_18_30as_grid, pop_nab_18_30as_grid, method = 'ngb')

# resample dem grid
slope_nab_resampled_pred <- resample(slope_nab_3as_grid, pop_nab_18_30as_grid, method = 'ngb')



stacked_nab_pred_set1 <- stack(pop_nab_02_30as_grid, pop_nab_05_30as_grid, pop_nab_08_30as_grid, pop_nab_11_30as_grid, pop_nab_14_30as_grid, pop_nab_17_30as_grid, 
                               ntl_nab_02_resampled_pred, ntl_nab_05_resampled_pred, ntl_nab_08_resampled_pred, ntl_nab_11_resampled_pred, 
                               ntl_nab_14_resampled_pred, ntl_nab_17_resampled_pred, slope_nab_resampled_pred)

stacked_nab_pred_set2 <- stack(pop_nab_03_30as_grid, pop_nab_06_30as_grid, pop_nab_09_30as_grid, pop_nab_12_30as_grid, pop_nab_15_30as_grid, pop_nab_18_30as_grid, 
                               ntl_nab_03_resampled_pred, ntl_nab_06_resampled_pred, ntl_nab_09_resampled_pred, ntl_nab_12_resampled_pred, 
                               ntl_nab_15_resampled_pred, ntl_nab_18_resampled_pred, slope_nab_resampled_pred)

stacked_nab_df_pred_set1 <- as.data.frame(stacked_nab_pred_set1, xy = TRUE)
stacked_nab_df_pred_set2 <- as.data.frame(stacked_nab_pred_set2, xy = TRUE)


colnames(stacked_nab_df_pred_set1) <- c("x", "y", "pop_nab_02_30as_grid", "pop_nab_05_30as_grid", "pop_nab_08_30as_grid", 
                                        "pop_nab_11_30as_grid", "pop_nab_14_30as_grid", "pop_nab_17_30as_grid", 
                                        "ntl_nab_02_30as_grid", "ntl_nab_05_30as_grid", "ntl_nab_08_30as_grid", "ntl_nab_11_30as_grid", 
                                        "ntl_nab_14_30as_grid", "ntl_nab_17_30as_grid", "slope_nab_30as_grid")

colnames(stacked_nab_df_pred_set2) <- c("x", "y", "pop_nab_03_30as_grid", "pop_nab_06_30as_grid", "pop_nab_09_30as_grid", 
                                        "pop_nab_12_30as_grid", "pop_nab_15_30as_grid", "pop_nab_18_30as_grid", 
                                        "ntl_nab_03_30as_grid", "ntl_nab_06_30as_grid", "ntl_nab_09_30as_grid", "ntl_nab_12_30as_grid_grid", 
                                        "ntl_nab_15_30as_grid", "ntl_nab_18_30as_grid", "slope_nab_30as_grid")

stacked_nab_df_pred_set1[is.na(stacked_nab_df_pred_set1)] <- 0
stacked_nab_df_pred_set2[is.na(stacked_nab_df_pred_set2)] <- 0


# import border map
nab_border <- readOGR(dsn="~/data",layer="nab_border")
nab_borderline <- readOGR(dsn="~/data",layer="nab_borderline")
nab_border_df <- as.data.frame(nab_border, xy = TRUE)
nab_border_df <- subset(nab_border_df, select = c(coords.x1, coords.x2))

### fitting spatiotemporal model
# make mesh
nab.loc <- cbind(stacked_nab_df$x, stacked_nab_df$y)
nab.bound <- inla.nonconvex.hull(nab.loc, concave = 0.1, resolution = c(20, 20))
nab.mesh <- inla.mesh.2d(boundary=nab.bound, max.edge=c(0.02, 0.04), cutoff = 0.01) 
nab.spde <- inla.spde2.matern(mesh=nab.mesh, alpha=2, 
                              B.tau=matrix(c(0, 1, 0),nrow=1,ncol=3),
                              B.kappa=matrix(c(0, 0, 1),nrow=1,ncol=3), 
                              prior.variance.nominal = 1,
                              theta.prior.prec=0.1)

nab.fixed.vague <- list(mean.intercept = 0, 
                        prec.intercept = 0, 
                        mean = 0,
                        prec = 0.001)

nab.prec.vague <- list(initial=4, prior="loggamma", 
                       param=c(1, 0.00005), fixed=FALSE)

nab.rho.vague <- list(rho = list(prior = 'normal', param = c(0, 0.15)))

par(mfrow=c(1,1), mar = c(1, 1, 1, 1))
plot(nab.mesh,main="",asp=1)
points(cbind(stacked_nab_df$x, stacked_nab_df$y), pch=19, col=2, cex=0.5) 
lines(nab_borderline,lwd=4,col=1)

k <- 6 # no. time stages
n <- as.numeric(dim(stacked_nab_df)[1])

nab.s.index.0 <- inla.spde.make.index(name="nab.spatial.field.0", n.spde=nab.spde$n.spde, n.group = k) 
nab.s.index.1 <- inla.spde.make.index(name="nab.spatial.field.1", n.spde=nab.spde$n.spde, n.group = k) 

nab.A_0 <- inla.spde.make.A(nab.mesh, 
                            cbind(rep(stacked_nab_df[, 1], k), rep(stacked_nab_df[, 2], k)),
                            group = rep(1:k, each = n))
nab.A_1 <- inla.spde.make.A(nab.mesh,
                            cbind(rep(stacked_nab_df[, 1], k), rep(stacked_nab_df[, 2], k)),
                            group = rep(1:k, each = n), weights = as.vector(as.matrix(stacked_nab_df[,9:14])))

nab.x.res.pred <- dim(pop_nab_18_30as_grid)[2]
nab.y.res.pred <- dim(pop_nab_18_30as_grid)[1] 

nab.seq.x.grid.pred <- seq(from=pop_nab_18_30as_grid@extent[1],to=pop_nab_18_30as_grid@extent[2],length=nab.x.res.pred)
nab.seq.y.grid.pred <- seq(from=pop_nab_18_30as_grid@extent[3],to=pop_nab_18_30as_grid@extent[4],length=nab.y.res.pred)
nab.pred.grid.pred <- as.matrix(expand.grid(x=nab.seq.x.grid.pred,y=nab.seq.y.grid.pred))

n_pred <- as.numeric(dim(stacked_nab_df_pred_set1)[1]) # as set 1 and 2 have the same number of data points

nab.A_0_pred <- inla.spde.make.A(nab.mesh, 
                                 cbind(rep(nab.pred.grid.pred[, 1], k), rep(nab.pred.grid.pred[, 2], k)),
                                 group = rep(1:k, each = n_pred))

nab.A_1_pred_1 <- inla.spde.make.A(nab.mesh, 
                                   cbind(rep(nab.pred.grid.pred[, 1], k), rep(nab.pred.grid.pred[, 2], k)),
                                   group = rep(1:k, each = n_pred), weights = as.vector(as.matrix(stacked_nab_df_pred_set1[,9:14])))

nab.A_1_pred_2 <- inla.spde.make.A(nab.mesh, 
                                   cbind(rep(nab.pred.grid.pred[, 1], k), rep(nab.pred.grid.pred[, 2], k)),
                                   group = rep(1:k, each = n_pred), weights = as.vector(as.matrix(stacked_nab_df_pred_set2[,9:14])))

####################################
# log transformation + w/ ntl + w/ slope + w/ field.0 + w/ field.1
####################################
nab.stack.est <- inla.stack(
  data = list(logpop = log(as.vector(as.matrix(stacked_nab_df[,3:8]))+1)), 
  A = list(nab.A_0, nab.A_1, 1), 
  effects = list(nab.s.index.0, nab.s.index.1, 
                 data.frame(Intercept = 1, 
                            ntl = as.vector(as.matrix(stacked_nab_df[,9:14])), 
                            slope = rep(stacked_nab_df$slope_nab_30as,k))),
  tag = "est") 

nab.stack.pred.set1 <- inla.stack(
  data = list(logpop = NA), 
  A = list(nab.A_0_pred, nab.A_1_pred_1, 1), 
  effects = list(nab.s.index.0, nab.s.index.1,
                 data.frame(Intercept = 1, 
                            ntl = as.vector(as.matrix(stacked_nab_df_pred_set1[,9:14])),
                            slope = rep(stacked_nab_df_pred_set1$slope_nab_30as_grid,k))),
  tag = "pred.response.1") 

nab.stack.pred.set2 <- inla.stack(
  data = list(logpop = NA), 
  A = list(nab.A_0_pred, nab.A_1_pred_2, 1), 
  effects = list(nab.s.index.0, nab.s.index.1,
                 data.frame(Intercept = 1, 
                            ntl = as.vector(as.matrix(stacked_nab_df_pred_set2[,9:14])),
                            slope = rep(stacked_nab_df_pred_set2$slope_nab_30as_grid,k))),
  tag = "pred.response.2")

nab.join.stack.pred <- inla.stack(nab.stack.est, nab.stack.pred.set1, nab.stack.pred.set2)

nab.formula <- logpop ~ -1 + Intercept + ntl + slope +
  f(nab.spatial.field.0, model = nab.spde, group = nab.spatial.field.0.group, 
    control.group = list(model = 'ar1', hyper = nab.rho.vague)) + 
  f(nab.spatial.field.1, model = nab.spde, group = nab.spatial.field.1.group, 
    control.group = list(model = 'ar1', hyper = nab.rho.vague))

nab.output.pred <- inla(nab.formula, family="gaussian",
                        data = inla.stack.data(nab.join.stack.pred, spde=nab.spde), 
                        control.fixed=nab.fixed.vague, 
                        control.family=list(hyper=list(prec=nab.prec.vague)), 
                        control.predictor = list(A = inla.stack.A(nab.join.stack.pred), compute=TRUE), 
                        control.inla = list(int.strategy="eb"),
                        control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE), verbose=TRUE)

summary(nab.output.pred)
nab.output.pred$dic$dic # 4364.769
nab.output.pred$cpu.used
slcpo(nab.output.pred)

# extract parameters for final model
nab.fixed.out <- round(nab.output.pred$summary.fixed[,1:5],3)
nab.fixed1_marg <- nab.output.pred$marginals.fixed[[1]]
nab.fixed2_marg <- nab.output.pred$marginals.fixed[[2]]
nab.fixed3_marg <- nab.output.pred$marginals.fixed[[3]]

# variance of unstructured residual
nab.sigma2e_marg <- inla.tmarginal(function(x) 1/x, nab.output.pred$marginals.hyperpar[[1]])
nab.sigma2e_m1 <- inla.emarginal(function(x) x, nab.sigma2e_marg)
nab.sigma2e_m2 <- inla.emarginal(function(x) x^2, nab.sigma2e_marg)
nab.sigma2e_stdev <- sqrt(nab.sigma2e_m2 - nab.sigma2e_m1^2)
nab.sigma2e_quantiles <- inla.qmarginal(c(0.025, 0.5, 0.975), nab.sigma2e_marg)

# temporal effect for field 0
nab.0.rho_marg <- nab.output.pred$marginals.hyperpar[[4]]
nab.0.rho_m1 <- inla.emarginal(function(x) x, nab.0.rho_marg)
nab.0.rho_m2 <- inla.emarginal(function(x) x^2, nab.0.rho_marg)
nab.0.rho_stdev <- sqrt(nab.0.rho_m2 - nab.0.rho_m1^2)
nab.0.rho_quantiles <- inla.qmarginal(c(0.025, 0.5, 0.975), nab.0.rho_marg)

# temporal effect for field 1
nab.1.rho_marg <- nab.output.pred$marginals.hyperpar[[7]]
nab.1.rho_m1 <- inla.emarginal(function(x) x, nab.1.rho_marg)
nab.1.rho_m2 <- inla.emarginal(function(x) x^2, nab.1.rho_marg)
nab.1.rho_stdev <- sqrt(nab.1.rho_m2 - nab.1.rho_m1^2)
nab.1.rho_quantiles <- inla.qmarginal(c(0.025, 0.5, 0.975), nab.1.rho_marg)

# field 0
nab.spde.result.0 <- inla.spde2.result(nab.output.pred, name="nab.spatial.field.0", nab.spde)

# variance of spatial effects
nab.0.var.nom.marg <- nab.spde.result.0$marginals.variance.nominal[[1]]
nab.0.var.nom.m1 <- inla.emarginal(function(x) x, nab.0.var.nom.marg)
nab.0.var.nom.m2 <- inla.emarginal(function(x) x^2, nab.0.var.nom.marg)
nab.0.var.nom.stdev <- sqrt(nab.0.var.nom.m2 - nab.0.var.nom.m1^2)
nab.0.var.nom.quantiles <- inla.qmarginal(c(0.025, 0.5, 0.975), nab.0.var.nom.marg)

nab.0.range.nom.marg <- nab.spde.result.0$marginals.range.nominal[[1]]
nab.0.range.nom.m1 <- inla.emarginal(function(x) x, nab.0.range.nom.marg)
nab.0.range.nom.m2 <- inla.emarginal(function(x) x^2, nab.0.range.nom.marg)
nab.0.range.nom.stdev <- sqrt(nab.0.range.nom.m2 - nab.0.range.nom.m1^2)
nab.0.range.nom.quantiles <- inla.qmarginal(c(0.025, 0.5, 0.975), nab.0.range.nom.marg)


# field 1, if applied
nab.spde.result.1 <- inla.spde2.result(nab.output.pred, name="nab.spatial.field.1", nab.spde)

# variance of spatial effects
nab.1.var.nom.marg <- nab.spde.result.1$marginals.variance.nominal[[1]]
nab.1.var.nom.m1 <- inla.emarginal(function(x) x, nab.1.var.nom.marg)
nab.1.var.nom.m2 <- inla.emarginal(function(x) x^2, nab.1.var.nom.marg)
nab.1.var.nom.stdev <- sqrt(nab.1.var.nom.m2 - nab.1.var.nom.m1^2)
nab.1.var.nom.quantiles <- inla.qmarginal(c(0.025, 0.5, 0.975), nab.1.var.nom.marg)

nab.1.range.nom.marg <- nab.spde.result.1$marginals.range.nominal[[1]]
nab.1.range.nom.m1 <- inla.emarginal(function(x) x, nab.1.range.nom.marg)
nab.1.range.nom.m2 <- inla.emarginal(function(x) x^2, nab.1.range.nom.marg)
nab.1.range.nom.stdev <- sqrt(nab.1.range.nom.m2 - nab.1.range.nom.m1^2)
nab.1.range.nom.quantiles <- inla.qmarginal(c(0.025, 0.5, 0.975), nab.1.range.nom.marg)

# plot posterior distribution of parameters
par(mfrow=c(5,2), mar = c(5, 5, 2, 2))

plot.j1.xmin <- min(nab.fixed1_marg[,1])
plot.j1.xmax <- max(nab.fixed1_marg[,1])
plot.j1.ymin <- min(nab.fixed1_marg[,2])
plot.j1.ymax <- max(nab.fixed1_marg[,2])
plot(nab.fixed1_marg,xlim=c(plot.j1.xmin,plot.j1.xmax),ylim=c(plot.j1.ymin,plot.j1.ymax),
     col="red",t="l",lty=1,xlab=expression(beta[0]),ylab="")

plot.j2.xmin <- min(nab.fixed2_marg[,1])
plot.j2.xmax <- max(nab.fixed2_marg[,1])
plot.j2.ymin <- min(nab.fixed2_marg[,2])
plot.j2.ymax <- max(nab.fixed2_marg[,2])
plot(nab.fixed2_marg,xlim=c(plot.j2.xmin,plot.j2.xmax),ylim=c(plot.j2.ymin,plot.j2.ymax),
     col="red",t="l",lty=1,xlab=expression(beta[1]),ylab="")

plot.j3.xmin <- min(nab.fixed3_marg[,1])
plot.j3.xmax <- max(nab.fixed3_marg[,1])
plot.j3.ymin <- min(nab.fixed3_marg[,2])
plot.j3.ymax <- max(nab.fixed3_marg[,2])
plot(nab.fixed3_marg,xlim=c(plot.j3.xmin,plot.j3.xmax),ylim=c(plot.j3.ymin,plot.j3.ymax),
     col="red",t="l",lty=1,xlab=expression(beta[2]),ylab="")

plot.j4.xmin <- min(nab.sigma2e_marg[,1])
plot.j4.xmax <- max(nab.sigma2e_marg[,1])
plot.j4.ymin <- min(nab.sigma2e_marg[,2])
plot.j4.ymax <- max(nab.sigma2e_marg[,2])
plot(nab.sigma2e_marg,xlim=c(plot.j4.xmin,plot.j4.xmax),ylim=c(plot.j4.ymin,plot.j4.ymax),
     col="red",t="l",lty=1,xlab=expression(sigma[e]^2),ylab="")

plot.j5.xmin <- min(nab.0.rho_marg[,1])
plot.j5.xmax <- max(nab.0.rho_marg[,1])
plot.j5.ymin <- min(nab.0.rho_marg[,2])
plot.j5.ymax <- max(nab.0.rho_marg[,2])
plot(nab.0.rho_marg,xlim=c(plot.j5.xmin,plot.j5.xmax),ylim=c(plot.j5.ymin,plot.j5.ymax),
     col="red",t="l",lty=1,xlab=expression(rho[0]),ylab="")

plot.j8.xmin <- min(nab.1.rho_marg[,1])
plot.j8.xmax <- max(nab.1.rho_marg[,1])
plot.j8.ymin <- min(nab.1.rho_marg[,2])
plot.j8.ymax <- max(nab.1.rho_marg[,2])
plot(nab.1.rho_marg,xlim=c(plot.j8.xmin,plot.j8.xmax),ylim=c(plot.j8.ymin,plot.j8.ymax),
     col="red",t="l",lty=1,xlab=expression(rho[1]),ylab="")

plot.j6.xmin <- min(nab.0.var.nom.marg[,1])
plot.j6.xmax <- max(nab.0.var.nom.marg[,1])
plot.j6.ymin <- min(nab.0.var.nom.marg[,2])
plot.j6.ymax <- max(nab.0.var.nom.marg[,2])
plot(nab.0.var.nom.marg,xlim=c(plot.j6.xmin,plot.j6.xmax),ylim=c(plot.j6.ymin,plot.j6.ymax),
     col="red",t="l",lty=1,xlab=expression(sigma[0]^2),ylab="")

plot.j9.xmin <- min(nab.1.var.nom.marg[,1])
plot.j9.xmax <- max(nab.1.var.nom.marg[,1])
plot.j9.ymin <- min(nab.1.var.nom.marg[,2])
plot.j9.ymax <- max(nab.1.var.nom.marg[,2])
plot(nab.1.var.nom.marg,xlim=c(plot.j9.xmin,plot.j9.xmax),ylim=c(plot.j9.ymin,plot.j9.ymax),
     col="red",t="l",lty=1,xlab=expression(sigma[1]^2),ylab="")

plot.j7.xmin <- min(nab.0.range.nom.marg[,1])
plot.j7.xmax <- max(nab.0.range.nom.marg[,1])
plot.j7.ymin <- min(nab.0.range.nom.marg[,2])
plot.j7.ymax <- max(nab.0.range.nom.marg[,2])
plot(nab.0.range.nom.marg,xlim=c(plot.j7.xmin,plot.j7.xmax),ylim=c(plot.j7.ymin,plot.j7.ymax),
     col="red",t="l",lty=1,xlab=expression(r[0]),ylab="")

plot.j10.xmin <- min(nab.1.range.nom.marg[,1])
plot.j10.xmax <- max(nab.1.range.nom.marg[,1])
plot.j10.ymin <- min(nab.1.range.nom.marg[,2])
plot.j10.ymax <- max(nab.1.range.nom.marg[,2])
plot(nab.1.range.nom.marg,xlim=c(plot.j10.xmin,plot.j10.xmax),ylim=c(plot.j10.ymin,plot.j10.ymax),
     col="red",t="l",lty=1,xlab=expression(r[1]),ylab="")



nab.x.res <- dim(pop_nab_18_30as)[2]
nab.y.res <- dim(pop_nab_18_30as)[1] 
nab.seq.x.grid <- seq(from=pop_nab_18_30as@extent[1],to=pop_nab_18_30as@extent[2],length=nab.x.res)
nab.seq.y.grid <- seq(from=pop_nab_18_30as@extent[3],to=pop_nab_18_30as@extent[4],length=nab.y.res)
nab.pred.grid <- as.matrix(expand.grid(x=nab.seq.x.grid,y=nab.seq.y.grid))


# plot fitted latent fields at each time point
nab.prj <- inla.mesh.projector(nab.mesh, xlim = range(nab.pred.grid[, 1]),
                               ylim = range(nab.pred.grid[, 2]), dims = c(nab.x.res,nab.y.res)) 


nab.ov.proj <- over(SpatialPoints(nab.prj$lattice$loc, proj4string=crs.geo), nab_outline)

n.spde <- as.numeric(nab.spde$n.spde)

nab.spde.prj.field.0 <- lapply(1:k, function(j) {
  r <- inla.mesh.project(nab.prj,
                         nab.output.pred$summary.random$nab.spatial.field.0$mean[1:n.spde + (j - 1) * n.spde])
  r[is.na(nab.ov.proj)] <- NA
  return(r) 
})

nab.spde.prj.field.1 <- lapply(1:k, function(j) {
  r <- inla.mesh.project(nab.prj,
                         nab.output.pred$summary.random$nab.spatial.field.1$mean[1:n.spde + (j - 1) * n.spde])
  r[is.na(nab.ov.proj)] <- NA
  return(r) 
})

nab.zlm.field.0 <- range(unlist(nab.spde.prj.field.0), na.rm = TRUE)
nab.zlm.field.1 <- range(unlist(nab.spde.prj.field.1), na.rm = TRUE)

par(mfrow = c(3, 2), mar = c(0, 0, 1.5, 4))

for (j in 1:k) {
  image(x = nab.prj$x, y = nab.prj$y, z = nab.spde.prj.field.0[[j]], asp = 1, 
        xlab = '', zlim = nab.zlm.field.0, axes = FALSE, col = mapping.color.c(64),
        main = paste0("Year: ", 2000+j*3-2))
  image.plot(legend.only = TRUE, zlim = nab.zlm.field.0, legend.width = 0.8, legend.shrink = 0.5, col = mapping.color.c())
}


for (j in 1:k) {
  image(x = nab.prj$x, y = nab.prj$y, z = nab.spde.prj.field.1[[j]], asp = 1, 
        xlab = '', zlim = nab.zlm.field.1, axes = FALSE, col = mapping.color.c(64),
        main = paste0("Year: ", 2000+j*2-1))
  image.plot(legend.only = TRUE, zlim = nab.zlm.field.1, legend.width = 0.8, legend.shrink = 0.5, col = mapping.color.c())
}


####################################
# outcomes for group 1
####################################
nab.index.pred.set1 <- inla.stack.index(nab.join.stack.pred, "pred.response.1")$data
nab.post.mean.pred.set1 <- nab.output.pred$summary.linear.predictor[nab.index.pred.set1, "mean"] ## 1802 --> 1579
nab.post.sd.pred.set1 <- nab.output.pred$summary.linear.predictor[nab.index.pred.set1, "sd"]

nab.post.mean.pred.grid.set1 <- matrix(as.matrix(nab.post.mean.pred.set1), nrow = nab.x.res.pred)
nab.post.sd.pred.grid.set1 <- matrix(as.matrix(nab.post.sd.pred.set1), nrow = nab.x.res.pred)

# plot real pop of nairobi
nab.ov.30as.pred <- over(SpatialPoints(stacked_nab_df_pred_set1[,1:2], proj4string=crs.geo), nab_outline)
stacked_nab_df_pred_set1 <- cbind(stacked_nab_df_pred_set1, nab.ov.30as.pred)
stacked_nab_df_pred_set1_na <- stacked_nab_df_pred_set1[is.na(stacked_nab_df_pred_set1$FID),]
stacked_nab_df_pred_set1_na[,3:8] <- NA
stacked_nab_df_pred_set1_fit <- stacked_nab_df_pred_set1[!is.na(stacked_nab_df_pred_set1$FID),]
stacked_nab_df_pred_set1_full <- rbind(stacked_nab_df_pred_set1_fit,stacked_nab_df_pred_set1_na)

stack_nab_df_pred_sub_set1 <- stacked_nab_df_pred_set1_full[,c(1,2,3:8)]
stack_nab_df_pred_sub_order_set1 <- orderBy(~y + x, data=stack_nab_df_pred_sub_set1)

nab_pop_30as_grid.mat_set1 <- matrix(nrow = nab.x.res.pred , ncol = nab.y.res.pred*k) #empty matrix

# set 1 real pop
for (j in 1:k) {
  nab_pop_30as_grid.mat_set1[,((j-1)*nab.y.res.pred+1):(j*nab.y.res.pred)] <- 
    log(matrix(as.matrix(stack_nab_df_pred_sub_order_set1[,j+2]), nrow = nab.x.res.pred)+1)
}

par(mfrow=c(1,1), mar = c(5, 5, 2, 2))

for (j in 1:k) {
  image.plot(nab.seq.x.grid.pred,nab.seq.y.grid.pred,nab_pop_30as_grid.mat_set1[,((j-1)*nab.y.res.pred+1):(j*nab.y.res.pred)],
             xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T,legend.width = 0.8,legend.mar=5,
             zlim = range(c(0, 11.802)) #, main = paste0("Time: ", 1999+j*3)
  )
  contour(nab.seq.x.grid.pred,nab.seq.y.grid.pred,nab_pop_30as_grid.mat_set1[,((j-1)*nab.y.res.pred+1):(j*nab.y.res.pred)], add=T, lwd=2, labcex=1)
  
  lines(nab_borderline,lwd=4)
} 


is.na(nab.post.mean.pred.grid.set1) <- is.na(nab_pop_30as_grid.mat_set1)

par(mfrow=c(1,1), mar = c(5, 5, 2, 2))

for (j in 1:k) {
  image.plot(nab.seq.x.grid.pred,nab.seq.y.grid.pred,nab.post.mean.pred.grid.set1[,((j-1)*nab.y.res.pred+1):(j*nab.y.res.pred)],
             xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T,legend.width = 0.8,legend.mar=5,
             zlim = range(c(0, 11.802)) #, main = paste0("Time: ", 1999+j*3)
  )
  contour(nab.seq.x.grid.pred,nab.seq.y.grid.pred,nab.post.mean.pred.grid.set1[,((j-1)*nab.y.res.pred+1):(j*nab.y.res.pred)], add=T, lwd=2, labcex=1)
  lines(nab_borderline,lwd=4)
}

is.na(nab.post.sd.pred.grid.set1) <- is.na(nab_pop_30as_grid.mat_set1)

for (j in 1:k) {
  image.plot(nab.seq.x.grid.pred,nab.seq.y.grid.pred,nab.post.sd.pred.grid.set1[,((j-1)*nab.y.res.pred+1):(j*nab.y.res.pred)],
             xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T,legend.width = 0.8,legend.mar=5,
             zlim = range(c(0,0.5)) #, main = paste0("Time: ", 1999+j*3)
  )
  contour(nab.seq.x.grid.pred,nab.seq.y.grid.pred,nab.post.sd.pred.grid.set1[,((j-1)*nab.y.res.pred+1):(j*nab.y.res.pred)], add=T, lwd=2, labcex=1)
  lines(nab_borderline,lwd=4)
}


# calculate correlation and rmse
nab.pop_30as_grid.df_set1 <- data.frame(matrix(ncol = k, 
                                               nrow = nab.x.res.pred*nab.y.res.pred))

for (j in 1:k) {
  nab.pop_30as_grid.df_set1[,j] <- 
    as.vector(nab_pop_30as_grid.mat_set1[,((j-1)*nab.y.res.pred+1):(j*nab.y.res.pred)])
}

nab.pop_30as_grid.df_set1 <- nab.pop_30as_grid.df_set1[!is.na(nab.pop_30as_grid.df_set1$X1),]


nab.post.mean.pred.df.set1 <- data.frame(matrix(ncol = k, 
                                                nrow = nab.x.res.pred*nab.y.res.pred))

for (j in 1:k) {
  nab.post.mean.pred.df.set1[,j] <- 
    as.vector(nab.post.mean.pred.grid.set1[,((j-1)*nab.y.res.pred+1):(j*nab.y.res.pred)])
}

nab.post.mean.pred.df.set1 <- nab.post.mean.pred.df.set1[!is.na(nab.post.mean.pred.df.set1$X1),]

nab.post.mean.pred.df.set1.mod <- exp(nab.post.mean.pred.df.set1)-1
nab.post.mean.pred.df.set1.mod[nab.post.mean.pred.df.set1.mod<0]<-0

nab.correl.set1 <- data.frame(matrix(nrow=1,ncol=k))
for(j in 1:k){
  nab.correl.set1[,j] <- cor(nab.post.mean.pred.df.set1.mod[,j], (exp(nab.pop_30as_grid.df_set1[,j])-1))
}
nab.correl.set1

nab.RMSE.set1 <- data.frame(matrix(nrow=1,ncol=k))
for(j in 1:k){
  nab.RMSE.set1[,j] <- sqrt(mean(((exp(nab.pop_30as_grid.df_set1[,j])-1) - nab.post.mean.pred.df.set1.mod[,j])^2))
}
nab.RMSE.set1


####################################
# outcomes for group 2
####################################
nab.index.pred.set2 <- inla.stack.index(nab.join.stack.pred, "pred.response.2")$data
nab.post.mean.pred.set2 <- nab.output.pred$summary.linear.predictor[nab.index.pred.set2, "mean"] ## 1802 --> 1579
nab.post.sd.pred.set2 <- nab.output.pred$summary.linear.predictor[nab.index.pred.set2, "sd"]

nab.post.mean.pred.grid.set2 <- matrix(as.matrix(nab.post.mean.pred.set2), nrow = nab.x.res.pred)
nab.post.sd.pred.grid.set2 <- matrix(as.matrix(nab.post.sd.pred.set2), nrow = nab.x.res.pred)

# plot real pop of nairobi
stacked_nab_df_pred_set2 <- cbind(stacked_nab_df_pred_set2, nab.ov.30as.pred)
stacked_nab_df_pred_set2_na <- stacked_nab_df_pred_set2[is.na(stacked_nab_df_pred_set2$FID),]
stacked_nab_df_pred_set2_na[,3:8] <- NA
stacked_nab_df_pred_set2_fit <- stacked_nab_df_pred_set2[!is.na(stacked_nab_df_pred_set2$FID),]
stacked_nab_df_pred_set2_full <- rbind(stacked_nab_df_pred_set2_fit,stacked_nab_df_pred_set2_na)


stack_nab_df_pred_sub_set2 <- stacked_nab_df_pred_set2_full[,c(1,2,3:8)]
stack_nab_df_pred_sub_order_set2 <- orderBy(~y + x, data=stack_nab_df_pred_sub_set2)

nab_pop_30as_grid.mat_set2 <- matrix(nrow = nab.x.res.pred , ncol = nab.y.res.pred*k) #empty matrix


# set 2 real pop
for (j in 1:k) {
  nab_pop_30as_grid.mat_set2[,((j-1)*nab.y.res.pred+1):(j*nab.y.res.pred)] <- 
    log(matrix(as.matrix(stack_nab_df_pred_sub_order_set2[,j+2]), nrow = nab.x.res.pred)+1)
}

par(mfrow=c(1,1), mar = c(5, 5, 2, 2))

for (j in 1:k) {
  image.plot(nab.seq.x.grid.pred,nab.seq.y.grid.pred,nab_pop_30as_grid.mat_set2[,((j-1)*nab.y.res.pred+1):(j*nab.y.res.pred)],
             xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T,legend.width = 0.8,legend.mar=5,
             zlim = range(c(0, 11.802)) #, main = paste0("Time: ", 1999+j*3)
  )
  contour(nab.seq.x.grid.pred,nab.seq.y.grid.pred,nab_pop_30as_grid.mat_set2[,((j-1)*nab.y.res.pred+1):(j*nab.y.res.pred)], add=T, lwd=2, labcex=1)
  
  lines(nab_borderline,lwd=4)
} 


is.na(nab.post.mean.pred.grid.set2) <- is.na(nab_pop_30as_grid.mat_set2)

par(mfrow=c(1,1), mar = c(5, 5, 2, 2))

for (j in 1:k) {
  image.plot(nab.seq.x.grid.pred,nab.seq.y.grid.pred,nab.post.mean.pred.grid.set2[,((j-1)*nab.y.res.pred+1):(j*nab.y.res.pred)],
             xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T,legend.width = 0.8,legend.mar=5,
             zlim = range(c(0, 11.802)) #, main = paste0("Time: ", 1999+j*3)
  )
  contour(nab.seq.x.grid.pred,nab.seq.y.grid.pred,nab.post.mean.pred.grid.set2[,((j-1)*nab.y.res.pred+1):(j*nab.y.res.pred)], add=T, lwd=2, labcex=1)
  lines(nab_borderline,lwd=4)
}

is.na(nab.post.sd.pred.grid.set2) <- is.na(nab_pop_30as_grid.mat_set2)

for (j in 1:k) {
  image.plot(nab.seq.x.grid.pred,nab.seq.y.grid.pred,nab.post.sd.pred.grid.set2[,((j-1)*nab.y.res.pred+1):(j*nab.y.res.pred)],
             xlab="W-E (degree)", ylab="N-S (degree)", col=mapping.color.c(64), axes=T,legend.width = 0.8,legend.mar=5,
             zlim = range(c(0,0.5)) #, main = paste0("Time: ", 1999+j*3)
  )
  contour(nab.seq.x.grid.pred,nab.seq.y.grid.pred,nab.post.sd.pred.grid.set2[,((j-1)*nab.y.res.pred+1):(j*nab.y.res.pred)], add=T, lwd=2, labcex=1)
  lines(nab_borderline,lwd=4)
}


# calculate correlation and rmse
nab.pop_30as_grid.df_set2 <- data.frame(matrix(ncol = k, 
                                               nrow = nab.x.res.pred*nab.y.res.pred))

for (j in 1:k) {
  nab.pop_30as_grid.df_set2[,j] <- 
    as.vector(nab_pop_30as_grid.mat_set2[,((j-1)*nab.y.res.pred+1):(j*nab.y.res.pred)])
}

nab.pop_30as_grid.df_set2 <- nab.pop_30as_grid.df_set2[!is.na(nab.pop_30as_grid.df_set2$X1),]


nab.post.mean.pred.df.set2 <- data.frame(matrix(ncol = k, 
                                                nrow = nab.x.res.pred*nab.y.res.pred))

for (j in 1:k) {
  nab.post.mean.pred.df.set2[,j] <- 
    as.vector(nab.post.mean.pred.grid.set2[,((j-1)*nab.y.res.pred+1):(j*nab.y.res.pred)])
}

nab.post.mean.pred.df.set2 <- nab.post.mean.pred.df.set2[!is.na(nab.post.mean.pred.df.set2$X1),]

nab.post.mean.pred.df.set2.mod <- exp(nab.post.mean.pred.df.set2)-1
nab.post.mean.pred.df.set2.mod[nab.post.mean.pred.df.set2.mod<0]<-0

nab.correl.set2 <- data.frame(matrix(nrow=1,ncol=k))
for(j in 1:k){
  nab.correl.set2[,j] <- cor(nab.post.mean.pred.df.set2.mod[,j], (exp(nab.pop_30as_grid.df_set2[,j])-1))
}
nab.correl.set2

nab.RMSE.set2 <- data.frame(matrix(nrow=1,ncol=k))
for(j in 1:k){
  nab.RMSE.set2[,j] <- sqrt(mean(((exp(nab.pop_30as_grid.df_set2[,j])-1) - nab.post.mean.pred.df.set2.mod[,j])^2))
}
nab.RMSE.set2

####################################
# The End of This Project
####################################
