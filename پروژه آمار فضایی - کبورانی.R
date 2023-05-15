#Authors:
#Alireza Kabourani

###########################
###########################
library(geoR)
data('parana')
# load(file='C://Workspace//Datasets//parana.rda')
# parana
df = data.frame(x = parana$coords[,1], y = parana$coords[,2], z = parana$data)
head(df)

plot(parana$border, type = 'l')
points(df$x, df$y, type = 'n', xlab = "X", ylab = "Y")
text(df$x, df$y, df$z, cex = 0.6)

#check normality
stem(df$z)

library(ggplot2)
ggplot(as.data.frame(df),aes(x = z)) +
  geom_histogram(aes(y = ..density..),alpha = 0.8, fill = 'purple', bins = 20) +
  stat_function(fun = dnorm,
                args = list(mean = mean(df$z), sd = sd(df$z)),
                lwd = 1, lty = 'dashed', color = 'red') +
  theme_minimal()

qqnorm(df$z) 
qqline(df$z, col = 'red', lwd = 2)

shapiro.test(df$z) #Normal ast

#check stationary in mean
plot(df$x, df$z, ylab = "Z(x,y)", xlab = "X")
abline(lm(df$z ~ df$x), lty = 'dashed', lwd = 2, col = 'red')

plot(df$z, df$y, xlab = "Z(x,y)", ylab = "Y")
abline(lm(df$y ~ df$z), lty = 'dashed', lwd = 2, col = 'red')

model1 = lm(df$z ~ df$x + df$y)
summary(model1)
df_poly = df
df_poly$x2 = df_poly$x ^ 2
df_poly$y2 = df_poly$y ^ 2
df_poly$xy = df_poly$x * df_poly$y

model2 = lm(z~x+x2+y+y2, data = df_poly)
model3 = lm(z~., data = df_poly)

#models 
summary(model1)[9]
summary(model2)[9]
summary(model3)[9] #behtarin model

#data with Trend
df_detrend = df
par(mfrow = c(2,2))
plot(df$x, df$z, ylab = "Z(x,y)", xlab = "X", main = "data with trend")
abline(lm(df$z ~ df$x), lty = 'dashed', lwd = 2, col = 'red')
plot(df$z, df$y, xlab = "Z(x,y)", ylab = "Y", main = "data with trend")
abline(lm(df$y ~ df$z), lty = 'dashed', lwd = 2, col = 'red')
#remove Trend
df_detrend$z = df_detrend$z - predict(model3,df_poly[-3])
#data without trend
plot(df_detrend$x, df_detrend$z, ylab = "Z(x,y)", xlab = "X", main = "Detrended Data")
abline(lm(df_detrend$z ~ df_detrend$x), lty = 'dashed', lwd = 2, col = 'blue')
plot(df_detrend$z, df_detrend$y, xlab = "Z(x,y)", ylab = "Y", main = "Detrended Data")
abline(lm(df_detrend$y ~ df_detrend$z), lty = 'dashed', lwd = 2, col = 'blue')

par(mfrow = c(1,2))
plot(df_detrend$x, df_detrend$y, type = 'n', main = "Scatterplot Detrended Data")
text(df_detrend$x, df_detrend$y, round(df_detrend$z,3), cex = 0.5)
hist(df_detrend$z, col = 'purple', breaks = 15, main = "Histogram Detrended Data")

dev.off()

#spatiality Test
#################
#Moran Test
df_dist_inv = 1/as.matrix(dist(cbind(df_detrend$x, df_detrend$y)))
diag(df_dist_inv) = 0
df_dist_inv[1:4,1:4]
library(ape)
Moran.I(df_detrend$z, weight = df_dist_inv)

#Based on these results, we can reject the null hypothesis that there is
#zero spatial autocorrelation present in the variable parana at alpha = .05.
#that means parana data is Spatially Correlated

#stationary in variogram and covariogram 
library(sp)
coordinates(df) =~x+y
coordinates(df_detrend) =~x+y

library(gstat)
plot(variogram(z~1, df, cloud=TRUE))
plot(variogram(z~1, df_detrend, cloud=TRUE))

Upper = quantile(df$z, 0.75) + 1.5*(quantile(df$z, 0.75) - quantile(df$z, 0.25))
Lower = quantile(df$z, 0.25) - 1.5*(quantile(df$z, 0.75) - quantile(df$z, 0.25))

which(df$z > Upper)
which(df$z < Lower)

library(gstat)
hscat(z~1, df_detrend,seq(0,700,175), mirror = T)

#baresi hamsangardi
library(georob)
par(mfrow = c(2,2))
j = 0
angles = c(20, 45, 90, 135, 180)
for(i in angles){
  if(i == 180){break}
  j = j + 1
  v = variogram(z~1, data = df_detrend,  tol.ver = i,
                tol.hor = i)
  rv.da = sample.variogram(z~1, data = df_detrend,
                           locations = ~x + y, lag.dist.def = seq(0, 200, by = 10),
                           xy.angle.def = c(angles[j] , angles[j+1]))
  
  plot(v$dist, v$gamma, col = j, type = 'l',
       xlab = "Distance", ylab = "SemiVariogram",
       main = paste("Angle: ", i), lwd = 2, ylim = c(390, 900))
  points(rv.da$lag.dist, rv.da$gamma, type = "l", col = j, lty = "dashed", lwd = 1.5)
}

ggplot(data = as.data.frame(parana$border), aes(x = east, y = north)) +
  geom_path() +
  geom_point(data = data.frame(df),aes(x = x, y = y, color = z),
             size = 7, pch = 18) +
  scale_color_gradientn(colors = c("yellow",'blue')) +
  theme_void()

dev.off()
# variogram with WLS
parana.vario = variog(as.geodata(df_detrend), max.dist = 400)
plot(parana.vario)
parana.exp = variofit(parana.vario, cov.model = "exponential")
parana.mat  = variofit(parana.vario, cov.model = "matern")
parana.sph  = variofit(parana.vario, cov.model = "spherical")
parana.lin  = variofit(parana.vario, cov.model = "linear")
parana.wav  = variofit(parana.vario, cov.model = "wave")
parana.pex  = variofit(parana.vario, cov.model = "powered.exponential")
parana.gau = variofit(parana.vario, cov.model = "gaussian")

par(mfrow = c(2,3))
plot(parana.vario, main = "exponential")
lines(parana.exp, lwd = 2, col = 1)
text(210, 0,col = 'red', cex = 0.9, paste("minimised weighted sum of squares =",round(parana.exp$value)))

plot(parana.vario, main = "matern")
lines(parana.mat, lwd = 2, col = 2)
text(210, 0,col = 'red', cex = 0.9, paste("minimised weighted sum of squares =",round(parana.mat$value)))

plot(parana.vario, main = "spherical")
lines(parana.sph, lwd = 2, col = 3)
text(210, 0,col = 'red', cex = 0.9, paste("minimised weighted sum of squares =",round(parana.sph$value)))

plot(parana.vario, main = "linear")
lines(parana.lin, lwd = 2, col = 4)
text(210, 0,col = 'red', cex = 0.9, paste("minimised weighted sum of squares =",round(parana.lin$value)))

plot(parana.vario, main = "wave")
lines(parana.wav, lwd = 2, col = 5)
text(210, 0,col = 'red', cex = 0.9, paste("minimised weighted sum of squares =",round(parana.wav$value)))

plot(parana.vario, main = "powered.exponential")
lines(parana.pex, lwd = 2, col = 6)
text(210, 0,col = 'red', cex = 0.9, paste("minimised weighted sum of squares =",round(parana.pex$value)))
# pas Spherical behtarin model
dev.off()

vario_model_wls = vgm(psill = parana.sph$nugget + parana.sph$cov.pars[1],
                    range = parana.sph$cov.pars[2],
                    nugget = parana.sph$nugget,
                    model = "Sph")

#variogram with MLE
parana.ml.exp = likfit(as.geodata(df_detrend), ini = c(1000, 50), nug = 100,
                    cov.model = "exponential")
parana.ml.sph = likfit(as.geodata(df_detrend), ini = c(1000, 50), nug = 100,
                       cov.model = "spherical")
parana.ml.lin = likfit(as.geodata(df_detrend), ini = c(1000, 50), nug = 100,
                       cov.model = "linear")
parana.ml.wav = likfit(as.geodata(df_detrend), ini = c(1000, 50), nug = 100,
                       cov.model = "wave")
parana.ml.pex = likfit(as.geodata(df_detrend), ini = c(1000, 50), nug = 100,
                       cov.model = "powered.exponential")
parana.ml.gau = likfit(as.geodata(df_detrend), ini = c(1000, 50), nug = 100,
                       cov.model = "gaussian")

parana.ml.exp$cov.model
parana.ml.exp$BIC

parana.ml.sph$cov.model
parana.ml.sph$BIC

parana.ml.lin$cov.model
parana.ml.lin$BIC

parana.ml.wav$cov.model
parana.ml.wav$BIC

parana.ml.pex$cov.model
parana.ml.pex$BIC

parana.ml.gau$cov.model
parana.ml.gau$BIC

vario_model_mle = vgm(psill = parana.ml.sph$nugget + parana.ml.sph$cov.pars[1],
                      range = parana.ml.sph$cov.pars[2],
                      nugget = parana.ml.sph$nugget,
                      model = "Sph")

# vario_model_mle = vgm(psill = parana.ml.gau$nugget + parana.ml.gau$cov.pars[1],
#                       range = parana.ml.gau$cov.pars[2],
#                       nugget = parana.ml.gau$nugget,
#                       model = "Gau")
####robust model
library(georob)
rv = sample.variogram(z~1, data = df_detrend,
                      locations =~ x + y, lag.dist.def = seq(0, 300, by = 25))
plot(rv, type = "b")

#########################################
########### Kriging #####################
#########################################
head(df_detrend)
vario_model_wls
vario_model_mle

## Calculating CVMSE for ordinary method for data with Trend (WLS)
pred = c()
for (i in 1:length(df$z)){
  ordinary_method1 = krige(
    z~1,
    df[-i,],
    df[i,0],
    model = vario_model_wls       
  )
  pred = append(pred, ordinary_method1$var1.pred)
}

CVMSE1 = sum((pred - df$z)^2)/length(df$z)

## Calculating CVMSE for ordinary method for data without Trend (WLS)
pred = c()
for (i in 1:length(df_detrend$z)){
  ordinary_method1 = krige(
    z~1,
    df_detrend[-i,],
    df_detrend[i,0],
    model = vario_model_wls       
  )
  pred = append(pred, ordinary_method1$var1.pred)
}

CVMSE2 = sum((pred - df_detrend$z)^2)/length(df_detrend$z)

## Calculating CVMSE for university method for data without Trend (WLS)
pred = c()
for (i in 1:length(df$z)){
  uni_method1 = krige(
    z~x+y,
    df[-i,],
    df[i,0],
    model = vario_model_wls       
  )
  pred = append(pred, uni_method1$var1.pred)
}
CVMSE3 = sum((pred - df$z)^2)/length(df$z)

## Calculating CVMSE for ordinary method for data with Trend (MLE)
pred = c()
for (i in 1:length(df$z)){
  ordinary_method1 = krige(
    z~1,
    df[-i,],
    df[i,0],
    model = vario_model_mle       
  )
  pred = append(pred, ordinary_method1$var1.pred)
}

CVMSE4 = sum((pred - df$z)^2)/length(df$z)

## Calculating CVMSE for ordinary method for data without Trend (MLE)
pred = c()
for (i in 1:length(df_detrend$z)){
  ordinary_method1 = krige(
    z~1,
    df_detrend[-i,],
    df_detrend[i,0],
    model = vario_model_mle       
  )
  pred = append(pred, ordinary_method1$var1.pred)
}

CVMSE5 = sum((pred - df_detrend$z)^2)/length(df_detrend$z)

## Calculating CVMSE for universal method for data without Trend (MLE)
pred = c()
for (i in 1:length(df$z)){
  uni_method1 = krige(
    z~x+y,
    df[-i,],
    df[i,0],
    model = vario_model_mle       
  )
  pred = append(pred, uni_method1$var1.pred)
}
CVMSE6 = sum((pred - df$z)^2)/length(df$z)

CVMSE1 #ordinary with raw Data (WLS)
CVMSE2 #ordinary with detrended data (WLS)
CVMSE3 #Universal (WLS)
CVMSE4 #ordinary with raw Data (MLE)
CVMSE5 #ordinary with detrended data (MLE)
CVMSE6 #Universal (MLE)
######################################
## pishbini noghat jadid ba behtarin krig ta inja
new_data = data.frame(x = round(runif(5, min(df$x), max(df$x))),
                      y = round(runif(5, min(df$y), max(df$y))))
new_data$x2 = new_data$x ^ 2
new_data$y2 = new_data$y ^ 2
new_data$xy = new_data$x * new_data$y
coordinates(new_data) = ~ x + y 

ordinary_method = krige(
  z~1,
  df_detrend,
  new_data[0],
  model = vario_model_mle   
)
ordinary_method
#ezafekardan dobare ravand
data.frame(x = new_data$x, y = new_data$y,
           predict = ordinary_method$var1.pred + predict(model3,new_data),
           var = ordinary_method$var1.var
)


######################################### Transform Krigging
shapiro.test(df$z)
shapiro.test(df_detrend$z)
#nyiaz be tranform kriging nist chon data normal ast

##############################
##### PLOT parana Kriging ####
##############################
grid = as.data.frame(df_detrend)[1:2]
grid[144,] = c(800, 550)
grid[145,] = c(100,50)
coordinates(grid) =~ x + y 
grid = makegrid(grid, cellsize = 3) 
colnames(grid) = c('x','y')
grid = locations.inside(grid, parana$border)
grid = SpatialPixelsDataFrame(points = grid[c("x", "y")] ,data = grid)

uni_method_vis = krige(
  z~x+y,
  df,
  grid,
  model = vario_model_wls       
)

library(ggplot2)
p1 = ggplot(data = as.data.frame(parana$border), aes(x = east, y = north)) +
  geom_polygon(fill = "white") + 
      geom_point(data = as.data.frame(uni_method_vis),
                aes(x = x, y = y, color = var1.pred),size = 2) +
       scale_color_gradient(low = c('#caf0f8', '#ade8f4'), high = c("#0077b6",'#001845')) +
        labs(title = "Rainfall prediction for Parana State, Brasil",
             subtitle = "Method: Universal Kriging with WLS variogram") +
       theme_void()
p2 = ggplot(data = as.data.frame(parana$border), aes(x = east, y = north)) +
  geom_polygon(fill = "white") + 
  geom_point(data = as.data.frame(uni_method_vis),
             aes(x = x, y = y, color = var1.var),size = 2) +
  scale_color_gradient(low = c('#00296b', '#00509d'), high = c("#fdc500",'#ffd500')) +
  ggtitle("Rainfall prediction var") +
  theme_void()
library(gridExtra)
grid.arrange(p1,p2, ncol = 2)

ord_method_vis = krige(
  z~1,
  df_detrend,
  grid,
  model = vario_model_mle       
)
grid2 = as.data.frame(grid)[1:2]
grid2$x2 = grid2$x ^ 2
grid2$y2 = grid2$y ^ 2
grid2$xy = grid2$x * grid2$y
ord_method_vis = data.frame(x = grid$x, y = grid$y,
           predicts = ord_method_vis$var1.pred + predict(model3,grid2),
           vars = ord_method_vis$var1.var
)

library(ggplot2)
p1 = ggplot(data = as.data.frame(parana$border), aes(x = east, y = north)) +
  geom_polygon(fill = "white") + 
  geom_point(data = as.data.frame(ord_method_vis),
             aes(x = x, y = y, color = predicts),size = 2) +
  scale_color_gradient(low = c('#caf0f8', '#ade8f4'), high = c("#0077b6",'#001845')) +
  labs(title = "Rainfall prediction for Parana State, Brasil",
       subtitle = "Method: Ordinary Kriging with MLE variogram") +
  theme_void()
p2 = ggplot(data = as.data.frame(parana$border), aes(x = east, y = north)) +
  geom_polygon(fill = "white") + 
  geom_point(data = as.data.frame(ord_method_vis),
             aes(x = x, y = y, color = vars),size = 2) +
  scale_color_gradient(low = c('#00296b', '#00509d'), high = c("#fdc500",'#ffd500')) +
  ggtitle("Rainfall prediction var") +
  theme_void()
library(gridExtra)
grid.arrange(p1,p2, ncol = 2)

##############################################
############ Robust ##########################
##############################################
library(DescTools)
sum(Winsorize(df_detrend$z, minval = quantile(df_detrend$z, 0.025),
              maxval = quantile(df_detrend$z, 0.975)) != df_detrend$z)

win_df = df_detrend
win_df$z = Winsorize(df_detrend$z, minval = quantile(df_detrend$z, 0.025),
                     maxval = quantile(df_detrend$z, 0.975))

par(mfrow = c(1,2))
hist(df_detrend$z, breaks = 20, freq = F, col = "cyan", border = "cyan",
     main = "Histogram of Rainfall", xlab = "Rainfall")
hist(win_df$z, col = 'purple', density = 20, freq = F, add = T)
legend("topright", pch = 15, legend = c("Not Winsorized","Winsorized"),
       col = c("cyan","purple"), cex = 0.8, box.lty = 0)
H = data.frame(z = c(df_detrend$z, win_df$z),
               Win = rep(c("Not Winsorized","Winsorized"), c(143,143)))
boxplot(z ~ Win, data = H, col = c('cyan','purple'))

#omit outliers
omit_df = df_detrend[which(win_df$z == df_detrend$z),]
par(mfrow = c(1,2))
hist(df_detrend$z, breaks = 20, freq = F, col = "yellow", border = "yellow",
     main = "Histogram of Rainfall", xlab = "Rainfall")
hist(omit_df$z, col = 'red', density = 20, freq = F, add = T)
legend("topright", pch = 15, legend = c("Complete Data","Omited Data"),
       col = c("yellow","red"), cex = 0.8, box.lty = 0)
H = data.frame(z = c(df_detrend$z, omit_df$z),
               Win = rep(c("Complete Data","Omited Data"), c(143,length(omit_df$z))))
boxplot(z ~ Win, data = H, col = c('yellow','red'))

#dev.off()


#Kriging winsorised data
library(automap)
vario_model_win = autofitVariogram(z~x+y, win_df)$var_model
plot(autofitVariogram(z~x+y, win_df))

pred = c()
for (i in 1:length(win_df$z)){
  uni_method = krige(
    z~x+y,
    win_df[-i,],
    win_df[i,],
    model = vario_model_win           
  )
  pred = append(pred, uni_method$var1.pred)
}
CVMSE_win_uni = sum((pred - win_df$z)^2)/length(win_df$z)
CVMSE_win_uni

#kriging omitted data
vario_model_omit = autofitVariogram(z~x+y, omit_df)$var_model
plot(autofitVariogram(z~x+y, omit_df))

pred = c()
for (i in 1:length(omit_df$z)){
  uni_method = krige(
    z~x+y,
    omit_df[-i,],
    omit_df[i,],
    model = vario_model_omit           
  )
  pred = append(pred, uni_method$var1.pred)
}
CVMSE_omit_uni = sum((pred - omit_df$z)^2)/length(omit_df$z)
CVMSE_omit_uni

omit_krig = krige(z ~ x+y,omit_df,grid, model = vario_model_omit)

omit_krig = data.frame(x = grid$x, y = grid$y,
                            predicts = omit_krig$var1.pred + predict(model3,grid2),
                            vars = omit_krig$var1.var
                       )
library(ggplot2)
p1 = ggplot(data = as.data.frame(parana$border), aes(x = east, y = north)) +
  geom_polygon(fill = "white") + 
  geom_point(data = as.data.frame(omit_krig),
             aes(x = x, y = y, color = predicts),size = 2) +
  scale_color_gradient(low = c('#caf0f8', '#ade8f4'), high = c("#0077b6",'#001845')) +
  labs(title = "Rainfall prediction for Parana State, Brasil",
       subtitle = "Method: Universal Kriging on Omitted Data") +
  theme_void()
p2 = ggplot(data = as.data.frame(parana$border), aes(x = east, y = north)) +
  geom_polygon(fill = "white") + 
  geom_point(data = as.data.frame(omit_krig),
             aes(x = x, y = y, color = vars),size = 2) +
  scale_color_gradient(low = c('#00296b', '#00509d'), high = c("#fdc500",'#ffd500')) +
  ggtitle("Rainfall prediction var") +
  theme_void()
library(gridExtra)
grid.arrange(p1,p2, ncol = 2)

aw = as.data.frame(omit_krig)
pred = cut(aw$predicts, breaks = 20,
           labels = hcl.colors(20, palette = "Viridis", rev = T), )
vars = cut(aw$predicts, breaks = 20,
           labels = hcl.colors(20, palette = "Viridis", rev = F), )

library(rgl)
open3d(windowRect=c(50,50,800,800))
plot3d(aw$x,aw$y,aw$predicts, col = pred, type = 'h')
play3d(spin3d())
# movie3d(spin3d(rpm = 10), duration = 5, fps = 10,
#         movie = "movie_pred", dir = "C://Workspace//")

plot3d(aw$x,aw$y,aw$vars, col = vars, type = 'h')
play3d(spin3d())
# movie3d(spin3d(rpm = 10), duration = 5, fps = 10,
#         movie = "movie_var", dir = "C://Workspace//")


CVMSE1 #ordinary with raw Data (WLS)
CVMSE2 #ordinary with detrended data (WLS)
CVMSE3 #Universal (WLS)
CVMSE4 #ordinary with raw Data (MLE)
CVMSE5 #ordinary with detrended data (MLE)
CVMSE6 #Universal (MLE)
CVMSE_omit_uni
CVMSE_win_uni

##############################
##### Indicator Kriging ######
##############################
ind_df = as.data.frame(df)
ind_df$ind = ifelse(ind_df$z >= mean(ind_df$z),1,0)
table(ind_df$ind)

library(sp)
coordinates(ind_df) =~x+y
head(ind_df)

library(automap)
vario_model_ind = autofitVariogram(ind~x+y, ind_df)$var_model
plot(autofitVariogram(ind~x+y, ind_df))

plot(ind_df$x, ind_df$y, col = ind_df$ind+16, pch = 19,
     ylim = c(0,500), xlim = c(100,800),
     xlab = "x", ylab = "y", cex = 2)
lines(parana$border)
legend("bottomright" , legend = c('rainfall > Mean','rainfall < Mean'),box.lty = 0,
       col = c(16,17),pch = 19, cex = 0.6, horiz = F,pt.cex = 1,
       xpd = F)

ind_method = krige(
  ind ~ x+y,
  ind_df,
  new_data,
  model = vario_model_ind,
  indicators = T
)
ind_method

ind_method_grided = krige(ind ~ x+y,ind_df,grid,
                          model = vario_model_ind,indicators = T)
plot(ind_method_grided, col = hcl.colors(15, palette = "Cividis", rev = T), main = "Indicator Kriging Predictions")
points(ind_df$x,ind_df$y,col= ind_df$ind,pch = 19,cex = 0.7)

pred = c()
for (i in 1:length(ind_df$ind)){
  ind_method = krige(
    ind~x+y,
    ind_df[-i,],
    ind_df[i,0],
    model = vario_model_ind
  )
  pred = append(pred, ind_method$var1.pred)
}
CVMSE_ind = sum((pred - ind_df$ind)^2)/length(ind_df$ind)
CVMSE_ind

##############################
##### probability kriging ####
##############################
dev.off()

library(MASS)
library(fitdistrplus)
fitn <- fitdist(df$z, "norm")
summary(fitn)
ks.test(df$z, 'pnorm', mean = fitn$estimate[1], sd = fitn$estimate[2] )

U = dnorm(df$z , fitn$estimate[1], fitn$estimate[2])
ind_df$U = U

g = gstat(NULL, id = "ind", form = ind ~ x+y, data = ind_df)
g = gstat(g, id = "U", form = U ~ x+y, data = ind_df)

v.cross = variogram(g)
plot(v.cross, pl=F)

g = gstat(g, id = "ind", model = vario_model_ind, fill.all=T)
g = fit.lmc(v.cross, g)
plot(variogram(g), model=g$model)

CK = predict(g, grid)
plot(CK)

CK_df = as.data.frame(CK)
head(CK_df)

library(ggplot2)
p1 = ggplot(data = CK_df, aes(x = x, y = y, color = ind.pred)) +
  geom_point(size = 5) +
  scale_color_gradient(low  = hcl.colors(15, palette = "Cividis", rev = T)) +
  theme_void()
p2 = ggplot(data = CK_df, aes(x = x, y = y, color = ind.var)) +
  geom_point(size = 5) +
  scale_color_gradient(low = c('black', 'red'), high = c("blue",'darkblue')) +
  theme_void()

library(gridExtra)
grid.arrange(p1,p2, ncol = 2)


#calculate CVMSE
models = vgm()[1]
models = as.character(models$short)
pred = c()
rm(CKcv)
defaultW = getOption("warn")
options(warn = -1)
for (i in seq(1,143)){
  gcv = gstat(NULL, id = "U", form = U ~ x+y, data = ind_df[-i,])
  gcv = gstat(gcv, id = "ind", form = ind ~ x+y, data = ind_df[-i,])
  vario_model_ind_cv = autofitVariogram(ind~1, ind_df[-i,])$var_model
  v.cross = variogram(gcv)
  gcv = gstat(gcv, id = "ind", model = vario_model_ind_cv, fill.all=T)
  gcv = fit.lmc(v.cross, gcv)
  
  t = try(expr = {CKcv = predict(gcv, ind_df[i,0])},silent=TRUE)
  
  if("try-error" %in% class(t)){
    j = 1
    t2 = t
    while ("try-error" %in% class(t2)) {
      j = j + 1
      gcv = gstat(NULL, id = "U", form = U ~ x+y, data = ind_df[-i,])
      gcv = gstat(gcv, id = "ind", form = ind ~ x+y, data = ind_df[-i,])
      vario_model_ind_cv = autofitVariogram(ind~x+y, ind_df[-i,],
                                            model = models[j])$var_model
      v.cross = variogram(gcv)
      gcv = gstat(gcv, id = "ind", model = vario_model_ind_cv, fill.all=T)
      gcv = fit.lmc(v.cross, gcv)
      t2 = try(expr = {CKcv = predict(gcv, ind_df[i,0])},silent=TRUE)
    }
  }
  if(is.na(CKcv$ind.pred)){print(i);break}
  pred = append(pred, CKcv$ind.pred)
}
options(warn = defaultW)
length(pred)

CVMSE_prob = sum((pred - ind_df$ind[-130])^2)/length(pred)

CVMSE_ind
CVMSE_prob



##############################
###### Baysian Kriging #######
##############################
grid = as.data.frame(df_detrend)[1:2]
grid[144,] = c(800, 550)
grid[145,] = c(100,50)
coordinates(grid) =~ x + y 
grid = makegrid(grid, cellsize = 10) 
colnames(grid) = c('x','y')
grid = locations.inside(grid, parana$border)
# grid = SpatialPixelsDataFrame(points = grid[c("x", "y")] ,data = grid)
grid2 = as.data.frame(grid)[1:2]
grid2$x2 = grid2$x ^ 2
grid2$y2 = grid2$y ^ 2
grid2$xy = grid2$x * grid2$y

#from previous analysis spherical is best model for this data set
parana.ml.sph
vario_model_mle
#be dalil inke trend as ghabl hazf shode: trend.d = "cte", trend.l = "cte"
#data parana normal ast pas lambda = 1
#khroji "bayes_krige" dar gozaresh avarde shode ast.
#be elat zaman bar boodan ravand ejra behtar ast ke in code ra ejra nakonid
bayes_cov_model = model.control(trend.d = "cte", trend.l = "cte",
                                cov.model = "spherical", lambda = 1)

bays_data = as.geodata(df_detrend)
#don't run takes too long
bayes_krige = krige.bayes(bays_data, coords = bays_data$coords,
                          data = bays_data$data,
                          locations = grid, borders = parana$border,
                          model = bayes_cov_model)

plot(bayes_krige)
  

bayes_krige$prior$beta$dist
bayes_krige$prior$sigmasq$dist
bayes_krige$prior$phi$dist
bayes_krige$prior$tausq.rel$dist

bayes_prior_model = post2prior(bayes_krige)
bayes_prior_model$beta.prior
bayes_prior_model$sigmasq.prior
bayes_prior_model$phi.prior
bayes_prior_model$tausq.rel.prior

bayes_krige_final = krige.bayes(bays_data, coords = bays_data$coords,
                                data = bays_data$data,
                                locations = grid, borders = parana$border,
                                model = bayes_cov_model,
                                prior = bayes_prior_model)

predictions = data.frame(x = grid$x, y = grid$y,
                         pred = as.vector(bayes_krige_final$predictive$mean.simulations)+ predict(model3,grid2),
                         var = as.vector(bayes_krige_final$predictive$variance.simulations))
head(predictions)
library(ggplot2)
p1 = ggplot(data = predictions, aes(x = x, y = y, color = pred)) +
  geom_point(size = 5, pch = 15) +
  scale_color_gradient(low = c('#caf0f8', '#ade8f4'), high = c("#0077b6",'#001845')) +
  theme_void()
p2 = ggplot(data = predictions, aes(x = x, y = y, color = var)) +
  geom_point(size = 5, pch = 15) +
  scale_color_gradient(low = c('black', 'red'), high = c("blue",'darkblue')) +
  theme_void()

library(gridExtra)
grid.arrange(p1,p2, ncol = 2)


#CVMSE
pred = c()
for (i in 1:143){
  loc = bays_data$coords[i,]
  bayes_krige_cv = krige.bayes(coords = bays_data$coords[-i,],
                               data = bays_data$data[-i],
                               locations = data.frame(x = loc[1], y = loc[2]),
                               model = bayes_cov_model,
                               prior = bayes_prior_model)
  
  pred = append(pred, bayes_krige_cv$predictive$mean.simulations)
}

CVMSE_bayes = sum((pred - bays_data$data)^2)/length(bays_data$data)
CVMSE_bayes
#637.2796

###########################################
###### Baysian Kriging Omitted data #######
###########################################
grid = as.data.frame(omit_df)[1:2]
grid[136,] = c(800, 550)
grid[137,] = c(100,50)
coordinates(grid) =~ x + y 
grid = makegrid(grid, cellsize = 10) 
colnames(grid) = c('x','y')
grid = locations.inside(grid, parana$border)
# grid = SpatialPixelsDataFrame(points = grid[c("x", "y")] ,data = grid)
grid2 = as.data.frame(grid)[1:2]
grid2$x2 = grid2$x ^ 2
grid2$y2 = grid2$y ^ 2
grid2$xy = grid2$x * grid2$y
#Matern with kappa 1.1 is best model for this data set
autofitVariogram(z~x+y, omit_df, 
                 model = c("Sph", "Exp", "Gau","Lin","Mat","Ste"))$var_model

model = lm(z~x+y, data = omit_df)
summary(model) 
#be dalil inke adam vojod trend: trend.d = "cte", trend.l = "cte"
#data parana normal ast pas lambda = 1
#khroji "bayes_krige" dar gozaresh avarde shode ast.
#be elat zaman bar boodan ravand ejra behtar ast ke in code ra ejra nakonid
bayes_cov_model = model.control(trend.d = "cte", trend.l = "cte",
                                cov.model = "matern",kappa = 1.1,
                                lambda = 1)

bays_data = as.geodata(omit_df)
#don't run takes too long
bayes_krige = krige.bayes(bays_data, coords = bays_data$coords,
                          data = bays_data$data,
                          locations = grid, borders = parana$border,
                          model = bayes_cov_model)
bayes_krige$model
bayes_krige$prior$beta$dist
bayes_krige$prior$sigmasq$dist
bayes_krige$prior$phi$dist
bayes_krige$prior$tausq.rel$dist
plot(bayes_krige)

bayes_prior_model = post2prior(bayes_krige)
bayes_prior_model$beta.prior
bayes_prior_model$sigmasq.prior
bayes_prior_model$phi.prior
bayes_prior_model$tausq.rel.prior

bayes_krige_final = krige.bayes(bays_data, coords = bays_data$coords,
                                data = bays_data$data,
                                locations = grid, borders = parana$border,
                                model = bayes_cov_model,
                                prior = bayes_prior_model)

plot(bayes_krige_final)

predictions = data.frame(x = grid$x, y = grid$y,
                         pred = as.vector(bayes_krige_final$predictive$mean.simulations)+ predict(model3,grid2),
                         var = as.vector(bayes_krige_final$predictive$variance.simulations))
head(predictions)
library(ggplot2)
p1 = ggplot(data = predictions, aes(x = x, y = y, color = pred)) +
  geom_point(size = 5, pch = 15) +
  scale_color_gradient(low = c('#caf0f8', '#ade8f4'), high = c("#0077b6",'#001845')) +
  theme_void()
p2 = ggplot(data = predictions, aes(x = x, y = y, color = var)) +
  geom_point(size = 5, pch = 15) +
  scale_color_gradient(low = c('black', 'red'), high = c("blue",'darkblue')) +
  theme_void()

library(gridExtra)
grid.arrange(p1,p2, ncol = 2)


#CVMSE
pred = c()
for (i in 1:length(bays_data$data)){
  loc = bays_data$coords[i,]
  bayes_krige_cv = krige.bayes(coords = bays_data$coords[-i,],
                               data = bays_data$data[-i],
                               locations = data.frame(x = loc[1], y = loc[2]),
                               model = bayes_cov_model,
                               prior = bayes_prior_model)
  
  pred = append(pred, bayes_krige_cv$predictive$mean.simulations)
}

CVMSE_bayes = sum((pred - bays_data$data)^2)/length(bays_data$data)
CVMSE_bayes #omitted
# 502.557
