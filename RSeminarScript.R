### R code from vignette source 'RSeminar2.Rnw'

setwd("")

###################################################
### code chunk number 11: RSeminar2.Rnw:193-195
###################################################
5L
as.integer(5.0)


###################################################
### code chunk number 12: RSeminar2.Rnw:198-201
###################################################
5.25
pi
as.numeric(1L)


###################################################
### code chunk number 13: RSeminar2.Rnw:210-212
###################################################
"hello world"
as.character(1.0)


###################################################
### code chunk number 14: RSeminar2.Rnw:215-218
###################################################
TRUE
1 > 2
as.logical(1L)


###################################################
### code chunk number 15: RSeminar2.Rnw:227-229
###################################################
NA
is.na(NA)


###################################################
### code chunk number 16: RSeminar2.Rnw:232-234
###################################################
NULL
is.null(NULL)


###################################################
### code chunk number 17: RSeminar2.Rnw:241-242 
###################################################
rm(list=ls())


###################################################
### code chunk number 18: RSeminar2.Rnw:244-248
###################################################
x1= 1/3
x2 <- 1/6
1/2 -> x3
class(x1)


###################################################
### code chunk number 19: RSeminar2.Rnw:252-253
###################################################
(x1+x2-2*x3)*0.1


###################################################
### code chunk number 20: RSeminar2.Rnw:260-263
###################################################
v1=c(1,2,4,6)
v2=1:4
( v3=0:3*5+1 )


###################################################
### code chunk number 21: RSeminar2.Rnw:267-270
###################################################
( v1 = c(x1=1,x2=2,x3=4,x4=6) )
names(v2) = c("x1","x2","x3","x4"); print(v2)
names(v3) = paste("x",1:4,sep=""); print(v3)


###################################################
### code chunk number 22: RSeminar2.Rnw:278-283
###################################################
v1[4]
v1["x4"]
v1[3:4]
v1[c(1,2,4)]
v1[-c(1,2,4)]


###################################################
### code chunk number 23: RSeminar2.Rnw:290-294
###################################################
( v4 = c("Hello","World") )
class(v4)
( v5 = rep(TRUE,4) )
class(v5)


###################################################
### code chunk number 24: RSeminar2.Rnw:300-301
###################################################
( fac1= factor(rep(1:3,3),labels=paste("Drug",1:3,sep="")) )


###################################################
### code chunk number 25: RSeminar2.Rnw:304-306
###################################################
levels(fac1) = paste("Drug",3:1,sep="")
fac1


###################################################
### code chunk number 26: RSeminar2.Rnw:313-316
###################################################
( mat1 = matrix(0L,2,3) )
( mat2 = matrix(1:6,2,3) )
( mat3 = matrix(1:6,2,3,byrow=TRUE) )


###################################################
### code chunk number 27: RSeminar2.Rnw:322-325
###################################################
rownames(mat2)=paste("row",1:2)
colnames(mat2)=paste("column",1:3)
mat2


###################################################
### code chunk number 28: RSeminar2.Rnw:331-335
###################################################
mat2[2,3]
mat2["row 2","column 3"]
mat2[6]
mat2[,c(1,3)]


###################################################
### code chunk number 29: RSeminar2.Rnw:343-350
###################################################
require("Matrix")
mat4 = matrix(0,1000,50)
mat5 = Matrix(0,1000,50,sparse=TRUE)
class(mat4)
object.size(mat4)
class(mat5)
object.size(mat5)


###################################################
### code chunk number 30: RSeminar2.Rnw:383-392 
###################################################
# require("bigalgebra")
# mat6 = matrix(1.0,5000,5000)
# mat7 = as.big.matrix(mat6)
# object.size(mat6);object.size(mat7)
# #200MB vs. 0.6KB
# system.time(mat6%*%mat6)
# system.time(mat7%*%mat7)


###################################################
### code chunk number 31: RSeminar2.Rnw:398-399
###################################################
( arr1 = array(1:12,c(2,3,2)) )


###################################################
### code chunk number 32: RSeminar2.Rnw:405-407
###################################################
arr1[2,1,1]
arr1[1,,2]


###################################################
### code chunk number 33: RSeminar2.Rnw:413-414
###################################################
set.seed(1) 


###################################################
### code chunk number 34: RSeminar2.Rnw:416-421
###################################################
( df1 = data.frame(id=1:6,
          gender=factor(rep(c("M","F"),3)),
          treatment1=factor(rep(LETTERS[1:3],2)),
          treatment2=factor(rep(LETTERS[1:2],each=3)),
          response=rnorm(6)) )


###################################################
### code chunk number 35: RSeminar2.Rnw:429-432
###################################################
df1$treatment1
df1[,"treatment1"]
df1[,3]


###################################################
### code chunk number 36: RSeminar2.Rnw:439-444
###################################################
ls1 = list()
ls1[[1]] = "Hello World"
ls1[[2]] = pi
ls1[[3]] = matrix(1:4,2,2)
print(ls1)


###################################################
### code chunk number 37: RSeminar2.Rnw:451-456
###################################################
DRfun = function(dose){
  return( 1-exp(-2.18E-04*dose) )
}
DRfun(3.18E+03)
DRfun(1+03)


###################################################
### code chunk number 38: RSeminar2.Rnw:462-466
###################################################
DRfun = function(dose,K){
  return( 1-exp(-K*dose) )
}
DRfun(1e3,2e-4)


###################################################
### code chunk number 39: RSeminar2.Rnw:473-474
###################################################
df1[which(df1$response<0),]


###################################################
### code chunk number 40: RSeminar2.Rnw:480-481
###################################################
with(df1,response[which(treatment1=="C")])


###################################################
### code chunk number 41: RSeminar2.Rnw:485-486
###################################################
with(df1,by(response,gender,mean))


###################################################
### code chunk number 42: RSeminar2.Rnw:502-505
###################################################
for(i in 1:4){
  print(i^2)
}


###################################################
### code chunk number 43: RSeminar2.Rnw:511-515
###################################################
x = c(1,2,4,8,16)
for(i in x){
  print(log2(i))
}


###################################################
### code chunk number 44: RSeminar2.Rnw:521-526
###################################################
count=0
while(count<5){
  print(count)
  count = count + 1
}


###################################################
### code chunk number 45: RSeminar2.Rnw:531-537
###################################################
set.seed(1)
x=0
while(x<0.5){
  x = rnorm(1)
  print(x)
}


###################################################
### code chunk number 46: RSeminar2.Rnw:543-550
###################################################
set.seed(1)
uu = runif(1)
if(uu < 0.5){
  print("heads")
}else{
  print("tails")
}


###################################################
### code chunk number 47: RSeminar2.Rnw:556-560
###################################################
set.seed(1)
uu = runif(1)
result = ifelse(uu < 0.5, "heads", "tails")
print(result)


###################################################
### code chunk number 48: RSeminar2.Rnw:593-597 
###################################################
sleepstudy = read.csv("http://myweb.uiowa.edu/dksewell/sleepstudy.csv")
test = read.table("http://myweb.uiowa.edu/dksewell/sleepstudy.csv",
                       sep=",",header=TRUE)
all.equal(sleepstudy,test)


###################################################
### code chunk number 49: RSeminar2.Rnw:611-612
###################################################
head(sleepstudy)


###################################################
### code chunk number 50: RSeminar2.Rnw:619-620
###################################################
tail(sleepstudy)


###################################################
### code chunk number 51: RSeminar2.Rnw:627-628 
###################################################
View(sleepstudy)  #In RStudio only


###################################################
### code chunk number 52: RSeminar2.Rnw:635-636 
###################################################
print(sleepstudy)


###################################################
### code chunk number 53: RSeminar2.Rnw:643-644
###################################################
colnames(sleepstudy)


###################################################
### code chunk number 54: RSeminar2.Rnw:651-652
###################################################
str(sleepstudy)


###################################################
### code chunk number 55: RSeminar2.Rnw:659-660
###################################################
dim(sleepstudy)


###################################################
### code chunk number 56: RSeminar2.Rnw:666-667
###################################################
blowdown = alr3::blowdown


###################################################
### code chunk number 57: RSeminar2.Rnw:669-670 
###################################################
?blowdown


###################################################
### code chunk number 58: RSeminar2.Rnw:675-676
###################################################
summary(blowdown)


###################################################
### code chunk number 59: RSeminar2.Rnw:681-682
###################################################
summary(blowdown$D)


###################################################
### code chunk number 60: RSeminar2.Rnw:687-688
###################################################
with(blowdown,by(D,y,summary))


###################################################
### code chunk number 61: RSeminar2.Rnw:693-695
###################################################
colMeans(blowdown[,-4])
#rowMeans()


###################################################
### code chunk number 62: RSeminar2.Rnw:700-701
###################################################
apply(blowdown[,-4],2,mean)


###################################################
### code chunk number 63: RSeminar2.Rnw:707-708
###################################################
apply(blowdown[,-4],2,function(x) mean(abs(x-mean(x))))


###################################################
### code chunk number 64: RSeminar2.Rnw:713-716
###################################################
blowdownScaled = scale(blowdown[,-4])
apply(blowdownScaled,2,
      function(x) round(c(mean=mean(x),sd=sd(x)),10))


###################################################
### code chunk number 65: RSeminar2.Rnw:723-726
###################################################
require("reshape2")
###Long form:
head(sleepstudy)


###################################################
### code chunk number 66: RSeminar2.Rnw:733-738
###################################################
###Wide form:
sleepWide = acast(sleepstudy,Subject~factor(Days),
                  value.var="Reaction")
dim(sleepWide)
head(round(sleepWide,1))


###################################################
### code chunk number 67: RSeminar2.Rnw:745-749
###################################################
###Back to tall form:
sleepTall = melt(sleepWide,value.name="Reaction")
dim(sleepTall)
head(round(sleepTall,1))


###################################################
### code chunk number 68: RSeminar2.Rnw:754-756 
###################################################
write.csv(sleepWide,file="sleepWide.csv")
write.table(sleepWide,file="sleepWide.txt",sep="\t")


###################################################
### code chunk number 69: RSeminar2.Rnw:771-774
###################################################
X = model.matrix(~Days,data=sleepstudy)
class(X)
dim(X)


###################################################
### code chunk number 70: RSeminar2.Rnw:780-781
###################################################
qr(X)$rank


###################################################
### code chunk number 71: RSeminar2.Rnw:787-789
###################################################
XtX = t(X)%*%X
all.equal(XtX,crossprod(X,X))  #Uses less memory


###################################################
### code chunk number 72: RSeminar2.Rnw:795-798
###################################################
XtXgenInv = MASS::ginv(XtX)
XtXInv = chol2inv(chol(XtX))
XtXInv = solve(XtX)


###################################################
### code chunk number 73: RSeminar2.Rnw:807-809
###################################################
XtXInv%*%t(X)%*%sleepstudy$Reaction
lm(Reaction~Days,data=sleepstudy)$coef


###################################################
### code chunk number 74: RSeminar2.Rnw:832-833 
###################################################
?USJudgeRatings


###################################################
### code chunk number 75: RSeminar2.Rnw:840-848 
###################################################
judgesCent = scale(as.matrix(USJudgeRatings[,-1]),scale=FALSE)
Sigma = cov(judgesCent)
eigs = eigen(Sigma)
Scores = judgesCent %*% eigs$vectors[,1:2]
plot(Scores,pch=16,ylim=c(-1,1)*2.5,xlab="",ylab="")
arrows(0,0,eigs$vec[,1]*3,eigs$vec[,2]*3)
text(eigs$vec[,1:2]*3,labels=colnames(judgesCent),
     adj=c(1,0))



###################################################
### code chunk number 77: RSeminar2.Rnw:897-899
###################################################
yMeans = as.numeric(with(blowdown,by(y,SPP,mean)))
barplot(yMeans,names.arg=levels(blowdown$SPP))


###################################################
### code chunk number 78: RSeminar2.Rnw:904-911 
###################################################
par(mfrow=c(1,2))
pie(yMeans,labels=levels(blowdown$SPP),
    main="Survival rates")
plotrix::pie3D(yMeans,labels=levels(blowdown$SPP),
               explode=0.1,
               col=rgb(20/256,c(1:9*20)/256,120/256))
par(mfrow=c(1,1))


###################################################
### code chunk number 81: RSeminar2.Rnw:931-932
###################################################
boxplot(blowdown$D)


###################################################
### code chunk number 82: RSeminar2.Rnw:937-938
###################################################
boxplot(blowdown$D,main="Blowdown data",ylab="Diameter")


###################################################
### code chunk number 83: RSeminar2.Rnw:943-947
###################################################
par(mfrow=c(1,2))
with(blowdown,by(D,y,boxplot,main="Blowdown data",
                 ylab="Diameter"))
par(mfrow=c(1,1))


###################################################
### code chunk number 84: RSeminar2.Rnw:952-955
###################################################
boxplot(blowdown$D~blowdown$y,main="Blowdown data",
        ylab="Diameter",names=c("survive","died"),
        cex.lab=1.5,cex.axis=1.5)


###################################################
### code chunk number 85: RSeminar2.Rnw:960-962
###################################################
hist(blowdown$D,main="Blowdown data",xlab="",
        ylab="Diameter",cex.lab=1.5,cex.axis=1.5)


###################################################
### code chunk number 87: RSeminar2.Rnw:974-976
###################################################
plot(density(blowdown$D),ylab="Diameter",cex.lab=1.5,
     cex.axis=1.5,main="blowdown data")


###################################################
### code chunk number 89: RSeminar2.Rnw:989-992
###################################################
hist(blowdown$D,main="Blowdown data",xlab="",freq=FALSE,
        ylab="Diameter",cex.lab=1.5,cex.axis=1.5)
lines(density(blowdown$D,bw=2.5),col="red",lwd=2)


###################################################
### code chunk number 90: RSeminar2.Rnw:997-999
###################################################
plot(ecdf(blowdown$D),main="Emp. Cum. Distn Func.",
     xlab="",cex.lab=1.5)


###################################################
### code chunk number 91: RSeminar2.Rnw:1010-1013
###################################################
DRfun = function(dose){
  return( 1-exp(-2.18E-04*dose) )
}


###################################################
### code chunk number 92: RSeminar2.Rnw:1024-1027
###################################################
curve(DRfun(x),from=0,to=50000,xlab="Dose",
      ylab="Probability of Infection",cex.lab=1.5,
      cex.axis=1.5,lwd=2)


###################################################
### code chunk number 94: RSeminar2.Rnw:1046-1052
###################################################
set.seed(123)
x1 = rnorm(100)
x2 = rgamma(100,5,10)
par(mfrow=c(1,2))
qqnorm(x1);qqline(x1)
qqnorm(x2);qqline(x2)



###################################################
### code chunk number 96: RSeminar2.Rnw:1072-1077
###################################################
par(mfrow=c(1,2))
qqplot(qgamma(ppoints(500),5,10),x1)
qqline(x1,distribution=function(p)qgamma(p,5,10))
qqplot(qgamma(ppoints(500),5,10),x2)
qqline(x2,distribution=function(p)qgamma(p,5,10))
par(mfrow=c(1,1))


###################################################
### code chunk number 97: RSeminar2.Rnw:1087-1090 
###################################################
class(ldeaths) #From datasets package
plot(ldeaths,ylab="# Deaths")
plot(stl(ldeaths,s.window="periodic"))


###################################################
### code chunk number 100: RSeminar2.Rnw:1128-1129
###################################################
pairs(airquality)


###################################################
### code chunk number 101: RSeminar2.Rnw:1136-1138
###################################################
plot(Ozone~Solar.R,data=airquality,xlab="Solar Radiation",
     ylab="Ozone",pch=16)


###################################################
### code chunk number 102: RSeminar2.Rnw:1145-1153 
###################################################
###################################################
### code chunk number 103: RSeminar2.Rnw:1161-1169
###################################################
CEXs = seq(0.5,5,length.out=500)
CEXs = with(airquality,
            CEXs[length(CEXs)*(max(Wind)+1-Wind)/
                   (max(Wind)+1)])
with(airquality,
     plot(Ozone~Solar.R,xlab="Solar Radiation",
     ylab="Ozone",pch=16,cex=CEXs,col=Month,
     cex.lab=1.5,cex.axis=1.5))


###################################################
### code chunk number 105: RSeminar2.Rnw:1189-1199
###################################################
with(airquality,sapply(unique(Month),
      function(x){ind=which(Month==x);
                  abline(lm(Ozone[ind]~Solar.R[ind]),
                         col=x,lwd=3)}))
with(airquality,legend("topleft",lwd=rep(4,5),cex=1.5,
       col=unique(Month),legend=unique(Month)))


###################################################
### code chunk number 106: RSeminar2.Rnw:1206-1213 
###################################################
plot(sleepTall$Reaction~sleepTall$Var2,pch=16,
     xlab="Days of sleep deprivation",
     ylab="Reaction time",cex.lab=1.5,cex.axis=1.5)
for(i in 1:nrow(sleepWide)){
  lines(sleepWide[i,]~c(0:9))
}
lines(colMeans(sleepWide)~c(0:9),col="blue",lwd=4,lty=2)


###################################################
### code chunk number 108: RSeminar2.Rnw:1234-1240 
###################################################
curve(dnorm(x,-1.5),-4.5,3.5)
xseq = seq(0,3.5,length.out=500)
yseq = dnorm(xseq,-1.5)
polygon(x=c(xseq,xseq[length(xseq):1],0),
        y=c(yseq,rep(0,length(xseq)+1)),
        col=rgb(0.25,0.75,1,alpha=0.5))


###################################################
### code chunk number 110: RSeminar2.Rnw:1260-1264 
###################################################
?text
?segments
?arrows
?symbols


###################################################
### code chunk number 111: RSeminar2.Rnw:1271-1276 
###################################################
?jpeg
?pdf
?png
?bmp
?tiff


###################################################
### code chunk number 112: RSeminar2.Rnw:1280-1283 
###################################################
## jpeg("foo.jpg",height=800,width=800)#in pixels
## plot( ... )
## dev.off()


###################################################
### code chunk number 114: RSeminar2.Rnw:1316-1323
###################################################
require("ggplot2")
blowdown$y <- 
  factor(blowdown$y,labels=c("Surv","Died"))
ggplot(blowdown,aes(factor(y),D))+
  geom_boxplot()+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=0,colour="white"))


###################################################
### code chunk number 116: RSeminar2.Rnw:1342-1346
###################################################
ggplot(blowdown,aes(D,fill=SPP))+
  geom_histogram()+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=0,colour="white"))


###################################################
### code chunk number 118: RSeminar2.Rnw:1365-1369
###################################################
ggplot(blowdown,aes(D,fill=SPP,colour=SPP))+
  geom_density(alpha=0.1)+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=0,colour="white"))


###################################################
### code chunk number 120: RSeminar2.Rnw:1388-1392
###################################################
ggplot(blowdown,aes(D,fill=SPP,colour=SPP))+
  geom_density(position="stack")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=0,colour="white"))


###################################################
### code chunk number 122: RSeminar2.Rnw:1411-1416
###################################################
ggplot(airquality,aes(Solar.R,Ozone))+
  geom_point(size=6,aes(color=Temp))+
  geom_smooth(method="lm")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=0,colour="white"))


###################################################
### code chunk number 113: RSeminar2.Rnw:1435-1449 
###################################################
require("emdbook")
require("plot3D")

normMix = function(x,y){
  mixProbs = c(1/3,1/2,1/6)
  ret = dnorm(x,-2)*dnorm(y,-2) +
    dnorm(x,0)*dnorm(y,2) +
    dnorm(x,2)*dnorm(y,0)
  return(ret)
}
curve3d(normMix(x,y),from=c(-5,-5),to=c(5,5),sys3d="contour",
        xlab="",ylab="",labcex=1.5,nlevels=20)
curve3d(normMix(x,y),from=c(-5,-5),to=c(5,5),sys3d="persp",theta=-15,
        xlab="",ylab="",zlab="")
curve3d(normMix(x,y),from=c(-5,-5),to=c(5,5),sys3d="rgl",
        xlab="",ylab="",zlab="",
        col = rgb(20/256,60/256,120/256,0.5))


###################################################
### code chunk number 116: RSeminar2.Rnw:1490-1496 
###################################################
require("plot3D");require("alr3")
pollution = sniffer[,c("TankTemp","TankPres")]
pollution = data.frame(pollution,logy=log(sniffer$Y))
with(pollution,scatter3D(TankTemp,TankPres,logy,theta=45,
                         phi=20,xlab="Temp",ylab="Press",
                         zlab="Pollution",pch=16))


###################################################
### code chunk number 118: RSeminar2.Rnw:1517-1521 
###################################################
with(pollution,scatter3D(TankTemp,TankPres,logy,theta=45,
                         phi=20,xlab="Temp",ylab="Press",
                         zlab="Pollution",pch=16,
                         type="h"))

###################################################
### code chunk number 120: RSeminar2.Rnw:1538-1552 
###################################################
fit = lm(logy~TankTemp+TankPres,data=pollution)
xgrid=with(pollution,
       seq(min(TankTemp),max(TankTemp),length.out=15))
ygrid=with(pollution,
       seq(min(TankPres),max(TankPres),length.out=15))
xygrid = expand.grid(TankTemp=xgrid,TankPres=ygrid)
logyPred = matrix(predict(fit,xygrid),15,15)
fitPts = predict(fit)
with(pollution,
     scatter3D(TankTemp,TankPres,logy,theta=45,
               phi=20,xlab="Temp",ylab="Press",
               zlab="Pollution",pch=16,
               surf=list(x=xgrid,y=ygrid,z=logyPred,
               fit=fitPts,facets=NA)))


###################################################
### code chunk number 123: RSeminar2.Rnw:1594-1600
###################################################
require("plot3Drgl")
BDTab = table(cut(blowdown$D,seq(5,85,by=10),
                  include.lowest=TRUE),
              blowdown$SPP)
hist3D(z=BDTab,col = rgb(20/256,60/256,120/256),
       border = "black",shade = 0.4,space = 0.15,
       xlab="Diameter",ylab="Species",zlab="")
plotrgl()


###################################################
### code chunk number 124: RSeminar2.Rnw:1640-1648
###################################################
RFun = function(n,start=1){
  ret=0
  for(i in start:n){
    ret = ret + i
  }
  return(ret)
}
RFun(100) == 100*101/2


###################################################
### code chunk number 125: RSeminar2.Rnw:1655-1659
###################################################
require("compiler")
RFunComp = cmpfun(RFun)
system.time(RFun(1e7))
system.time(RFunComp(1e7))


###################################################
### code chunk number 126: RSeminar2.Rnw:1666-1680
###################################################
require("foreach")
require("doParallel")
registerDoParallel(cl=2)
system.time({
  temp = foreach(i=c(5e6,1e7)) %dopar% {
    if(i==5e6){
      RFun(i)
    }else{
      RFun(i,5e6+1)
    }
  }
  print(sum(unlist(temp)) == 1e7*(1e7+1)/2)
})
stopImplicitCluster()


###################################################
### code chunk number 127: RSeminar2.Rnw:1687-1698
###################################################
require("Rcpp");require("inline")
cppFun = cxxfunction(signature(nFromR="integer"),
                     body='
int n = Rcpp::as<int>(nFromR);
int ret = 0;
for(int i=1;i<n+1;i++){
  ret += i;
}
return wrap(ret);
',plugin="Rcpp")
cppFun(100)


###################################################
### code chunk number 128: RSeminar2.Rnw:1703-1705
###################################################
system.time(RFun(1e7))
system.time(cppFun(1e7))


###################################################
### code chunk number 129: RSeminar2.Rnw:1734-1737
###################################################
attach(blowdown)
dim(blowdown)
head(blowdown)


###################################################
### code chunk number 130: RSeminar2.Rnw:1741-1742
###################################################
rm(list=setdiff(ls(),"blowdown"))


###################################################
### code chunk number 131: RSeminar2.Rnw:1746-1748
###################################################
y = factor(y,labels=c("Survived","Died"))
( tab = table(SPP,y) )


###################################################
### code chunk number 132: RSeminar2.Rnw:1754-1756
###################################################
chisq.test(tab)
prop.table(tab,1)


###################################################
### code chunk number 134: RSeminar2.Rnw:1775-1784
###################################################
image(t(prop.table(tab,1)),xaxt="n",yaxt="n",
      main="Tree Species and Survival")
box()
axis(1,at=0:1,labels=levels(y),cex.axis=1.5)
axis(2,at=seq(1,0,length.out=length(levels(SPP))),
     labels=levels(SPP),las=1,cex.axis=1.5)
text(cbind(rep(0:1,each=length(levels(SPP))),
          rep(seq(0,1,length.out=length(levels(SPP))),2)),
     labels=round(c(prop.table(tab,1)),2),cex=1.5)


###################################################
### code chunk number 135: RSeminar2.Rnw:1790-1792
###################################################
logD = log(D)
t.test(logD~y,var.equal=FALSE)


###################################################
### code chunk number 137: RSeminar2.Rnw:1806-1809
###################################################
par(mfrow=c(1,2))
by(logD,y,qqnorm)
par(mfrow=c(1,1))


###################################################
### code chunk number 138: RSeminar2.Rnw:1816-1817
###################################################
kruskal.test(logD,y)


###################################################
### code chunk number 140: RSeminar2.Rnw:1834-1839
###################################################
ind=which(y=="Died")
plot(density(x=logD[ind]),ylim=c(0,0.9),
     main="log(Diameter)")
lines(density(x=logD[-ind]),lty=2)
legend("topright",legend=c("Died","Survived"),lty=1:2)


###################################################
### code chunk number 141: RSeminar2.Rnw:1846-1847
###################################################
boxplot(logD~SPP)


###################################################
### code chunk number 142: RSeminar2.Rnw:1853-1855
###################################################
lmMod= lm(logD~SPP)
anova(lmMod)


###################################################
### code chunk number 143: RSeminar2.Rnw:1861-1862
###################################################
summary(lmMod)


###################################################
### code chunk number 145: RSeminar2.Rnw:1876-1878
###################################################
qqnorm(resid(lmMod));qqline(resid(lmMod))
plot(resid(lmMod)~fitted(lmMod));abline(h=0)


###################################################
### code chunk number 148: RSeminar2.Rnw:1899-1902
###################################################
plot(density(x=S),ylim=c(0,2.25),main="Storm Severity")
lines(density(x=S[-ind]),lty=2)
legend("topright",legend=c("Died","Survived"),lty=1:2)


###################################################
### code chunk number 149: RSeminar2.Rnw:1909-1910
###################################################
( logRegMod = glm(y~logD + SPP + S, family=binomial) )


###################################################
### code chunk number 150: RSeminar2.Rnw:1915-1916
###################################################
round(summary(logRegMod)$coef,4)


###################################################
### code chunk number 151: RSeminar2.Rnw:1924-1926
###################################################
( probRegMod = glm(y~logD + SPP + S, 
                   family=binomial(link="probit")) )


###################################################
### code chunk number 152: RSeminar2.Rnw:1932-1933
###################################################
summary(probRegMod)$coef


###################################################
### code chunk number 154: RSeminar2.Rnw:1958-1969
###################################################
logitFun = function(x,y){
  eta = logRegMod$coef["(Intercept)"] +
    logRegMod$coef["SPPBS"]+logRegMod$coef["logD"]*x+
    logRegMod$coef["S"]*y
  return(1/(1+exp(-eta)))
}
with(blowdown,
     curve3d(logitFun,from=c(min(logD),min(S)),
             to=c(max(logD),max(S)),sys3d="persp",
             theta=-20,xlab="Diameter",ylab="Storm",
             zlab="Probability of Dying"))


###################################################
### code chunk number 156: RSeminar2.Rnw:1995-1999
###################################################
yr <- floor(tt <- time(mdeaths))
plot(ldeaths,ylab="",
     xy.labels=paste(month.abb[12*(tt-yr)],
                     yr-1900,sep="'"))


###################################################
### code chunk number 157: RSeminar2.Rnw:2037-2053
###################################################
require("MCMCpack")
fun1= function(y,nsims=1000,stepsAhead=100){
  TT = length(y)
  BB = sum(y[-TT]^2)
  bb = sum(y[-1]*y[-TT])/BB
  Qb = sum((y[-1]-bb*y[-TT])^2)
  
  s2 = rinvgamma(nsims,shape=0.5*(TT-3),scale=0.5*Qb)
  phi = rnorm(nsims,mean=bb,sd=sqrt(s2/BB))
  ypred = matrix(0.0,nsims,stepsAhead)
  ypred[,1]=rnorm(nsims,phi*y[TT],sd=sqrt(s2))
  for(tt in 2:stepsAhead){
    ypred[,tt] = rnorm(nsims,phi*ypred[,tt-1],sd=sqrt(s2))
  }
  return(list(s2=s2,phi=phi,ypred=ypred))
}


###################################################
### code chunk number 158: RSeminar2.Rnw:2061-2077
###################################################
set.seed(1)
M=1000
phiTrue = 0.75
s2True = 0.5
TT=500
Y = matrix(NA,M,TT+1)
Y[,1]=rnorm(M,sd=sqrt(s2True/(1-phiTrue^2)))
phiHat = s2Hat = numeric(M)
for(tt in 1+1:TT){
  Y[,tt] = 0.75*Y[,tt-1] + rnorm(M,sd=sqrt(s2True))
}
for(iter in 1:M){
  temp = fun1(Y[iter,],stepsAhead=2)
  s2Hat[iter] = mean(temp$s2)
  phiHat[iter] = mean(temp$phi)
}


###################################################
### code chunk number 159: RSeminar2.Rnw:2083-2087 
###################################################
hist(phiHat)
abline(v=phiTrue,col="red",lwd=2)
hist(s2Hat)  
abline(v=s2True,col="red",lwd=2)


###################################################
### code chunk number 162: RSeminar2.Rnw:2107-2120 
###################################################
ldeaths1 = ldeaths-mean(ldeaths)
fit1 = fun1(ldeaths1,nsims=1e5,stepsAhead=12)
ypreds = colMeans(fit1$ypred)+mean(ldeaths)
yBounds = t(apply(fit1$ypred,2,quantile,
                probs=c(0.025,0.975)))+mean(ldeaths)
YL = range(c(ldeaths,yBounds))
plot(c(ldeaths),ylab="",xlab="",type="l",
     xlim=c(1,(length(ldeaths)+12)),ylim=YL,lwd=2)
lines(ypreds~c(length(ldeaths)+1:12),col="red",lwd=2)
lines(yBounds[,1]~c(length(ldeaths)+1:12),col="blue",
      lwd=2,lty=2)
lines(yBounds[,2]~c(length(ldeaths)+1:12),col="blue",
      lwd=2,lty=2)


###################################################
### code chunk number 165: RSeminar2.Rnw:2159-2165
###################################################
TT = length(ldeaths1);  BB = sum(ldeaths1[-TT]^2);  
bb = sum(ldeaths1[-1]*ldeaths1[-TT])/BB;  
Qb = sum((ldeaths1[-1]-bb*ldeaths1[-TT])^2)
curve(dgamma(x,shape=0.5*(TT-3),scale=0.5*Qb),
      lwd=2,xlab=expression(sigma^2),from=0,to=5e8,
      ylab="",cex.axis=1.5,cex.lab=1.5)


###################################################
### code chunk number 166: RSeminar2.Rnw:2170-2173
###################################################
curve(dnorm(x,bb,sd=sqrt(mean(fit1$s2)/BB)),
      lwd=2,xlab=expression(phi),from=0,to=1,
      ylab="",cex.axis=1.5,cex.lab=1.5)


###################################################
### code chunk number 167: RSeminar2.Rnw:2180-2182
###################################################
require("compiler")
fun1Comp = cmpfun(fun1)


###################################################
### code chunk number 168: RSeminar2.Rnw:2198-2242
###################################################
require("RcppArmadillo");require("inline")
funcpp = cxxfunction(
  signature(YY="numeric",NSIMS="integer",STEPSAHEAD="integer"),
  body='

//Environment stats("package:stats");
//Function rnorm = stats["rnorm"];
Environment MCMCpack("package:MCMCpack");
Function rinvgamma = MCMCpack["rinvgamma"];

arma::colvec yy = Rcpp::as<arma::colvec>(YY);
int nsims = Rcpp::as<int>(NSIMS);
int stepsAhead = Rcpp::as<int>(STEPSAHEAD);
int TT = yy.n_rows;
double BB=0, bb=0, Qb=0, sig=0;

for(int tt=1;tt<TT;tt++){
  BB += pow(yy(tt-1),2);
  bb += yy(tt)*yy(tt-1);
}
bb = bb/BB;

for(int tt=1;tt<TT;tt++){
  Qb += pow( yy(tt)-bb*yy(tt-1) ,2);
}

arma::colvec s2 = Rcpp::as<arma::colvec>(rinvgamma(nsims,0.5*(TT-3),0.5*Qb));
arma::colvec phi = arma::randn(nsims);
for(int i=0;i<nsims;i++){
  phi(i) = bb + sqrt(s2(i)/BB)*phi(i);
}

arma::mat ypreds = arma::randn(nsims,stepsAhead);
for(int i=0;i<nsims;i++){
sig = sqrt(s2(i));
ypreds(i,0) = phi(i)*yy(TT-1)+sig*ypreds(i,0);
for(int tt=1;tt<stepsAhead;tt++){
ypreds(i,tt) = phi(i)*ypreds(i,tt-1) + sig*ypreds(i,tt);
}
}

return Rcpp::List::create(wrap(ypreds),wrap(phi));

  ',plugin="RcppArmadillo")


###################################################
### code chunk number 169: RSeminar2.Rnw:2248-2260
###################################################
system.time({for(it in 1:100){
  Rsims = fun1(Y[1,])
  }
})
system.time({for(it in 1:100){
  RsimsComp = fun1Comp(Y[1,])
  }
})
system.time({for(it in 1:100){
  Rcppsims = funcpp(Y,NSIMS=1000,STEPSAHEAD=100)
  }
})


