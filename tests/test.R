#https://www.kaggle.com/datasets/joyshil0599/mlb-hitting-and-pitching-stats-through-the-years
dat<-read.csv("data/pitcher_data.csv")
dat$AVG<-as.numeric(dat$AVG)
dat1<-na.omit(dat)[1:500,]
mars_obj<-mars(Win~Games.played+Innings.pitched+Hit.Batsmen+base.on.balls+WHIP+AVG,dat1,mars.control(Mmax=8,d=3,trace=T))
summary(mars_obj)
plot(mars_obj)
pred_dat<-na.omit(dat)[500:600,]
predict(mars_obj,pred_dat)

#https://www.kaggle.com/datasets/neuromusic/avocado-prices
dat2<-read.csv("data/avocado.csv")[1:400,]
mars_obj2<-mars(AveragePrice~Total.Volume+X4046+X4225+X4770+Small.Bags+Large.Bags+XLarge.Bags+year,dat2)
print(mars_obj2)
anova(mars_obj2)

#https://www.kaggle.com/datasets/grubenm/austin-weather
dat3<-read.csv("data/austin_weather.csv")[1:400,]
dat3$PrecipitationSumInches<-as.numeric(dat3$PrecipitationSumInches)
dat3<-na.omit(dat3)
dat3<-dat3[-305,]
mars_obj3<-mars(PrecipitationSumInches~HumidityAvgPercent+TempAvgF+SeaLevelPressureAvgInches+WindAvgMPH,dat3)
summary(mars_obj3)
