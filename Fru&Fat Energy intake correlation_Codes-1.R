
# Clear the environment
rm(list=ls())

# Load the packages
library(plyr)
library(fields)
library(mgcv)
library(sp)


# Load the data
data<-read.csv("../data.csv") 

head(data)
str(data)

data$in.F<-as.factor(data$in.F) 

score_F<-sort(unique(data$in.F)) 
col_F<-c("black", "red","blue")
data$my_col<-col_F[match(data$in.F, score_F)]

par(mar=c(5,5,5,1))	
plot(data$in.Fru, data$energy.in, col=data$my_col, ylim=c(10, 70), xlim=c(-0.5, 8), pch=1, cex=0.4, xlab="Fructose (kJ/g)", ylab="Energy Intake (kJ/day)", cex.lab=1.8, cex.axis=1.2)

legend(6.03, 72, as.character(sort(unique(data$in.F))), col=col_F, pch=1, cex=1.4, title = "Fat (kJ/g)")

model<-gam(energy.in ~ s(in.Fru, by=in.F, k=4) + in.F, data=data)
summary(model)

x<-fitted(model, type="response")
y<-resid(model, type="pearson")
test<-gam(y ~s(x))
print(summary(test))

new.data<-data.frame(in.F=as.factor(1.43), in.Fru=seq(0, 7, 1))  
out<-predict(model, newdata=new.data, se.fit=T)


lines(new.data$in.Fru, out$fit, col=col_F[1], lwd=2)

lines(new.data$in.Fru, out$fit + out$se.fit, lty=2, col=col_F[1])
lines(new.data$in.Fru, out$fit - out$se.fit, lty=2, col=col_F[1])

new.data<-data.frame(in.F=as.factor(2.85), in.Fru=seq(0, 6, 1)) 
out<-predict(model, newdata=new.data, se.fit=T)

lines(new.data$in.Fru, out$fit, col=col_F[2], lwd=2)

lines(new.data$in.Fru, out$fit + out$se.fit, lty=2, col=col_F[2])
lines(new.data$in.Fru, out$fit - out$se.fit, lty=2, col=col_F[2])

new.data<-data.frame(in.F=as.factor(4.28), in.Fru=seq(0, 5, 1)) 
out<-predict(model, newdata=new.data, se.fit=T)

lines(new.data$in.Fru, out$fit, col=col_F[3], lwd=2)

lines(new.data$in.Fru, out$fit + out$se.fit, lty=2, col=col_F[3])
lines(new.data$in.Fru, out$fit - out$se.fit, lty=2, col=col_F[3])



