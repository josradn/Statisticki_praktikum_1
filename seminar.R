#ucitavanje podataka:
x=zad53r$x
x=x[-1]
x=as.numeric(x)
x
y=zad53r$X.
y=y[-1]
y
n=length(x)
n

#(a)
lny=log(y)
lny
lnx=log(x)
lnx
krozx=1/x
krozx

par(mfrow=c(3,2))
plot(x,y,xlab="težina tijela",ylab="odreðeno svojstvo metabolizma",main="originalni podaci: (x,y)")
#(1)
plot(x,lny,xlab="težina tijela",ylab="odreðeno svojstvo metabolizma",ylim=c(4,6),main="transformirani podaci: (x,lny)")
#(2)
plot(lnx,lny,xlab="težina tijela",ylab="odreðeno svojstvo metabolizma",xlim=c(4.5,6.5),ylim=c(4,6),main="transformirani podaci: (lnx,lny)")
#(3)
plot(lnx,y,xlab="težina tijela",ylab="odreðeno svojstvo metabolizma",xlim=c(4.5,6.5),main="transformirani podaci: (lnx,y)")
#(4)
plot(krozx,lny,xlab="težina tijela",ylab="odreðeno svojstvo metabolizma",ylim=c(4,6),main="transformirani podaci: (1/x,lny)")
#(5)
plot(krozx,y,xlab="težina tijela",ylab="odreðeno svojstvo metabolizma",main="transformirani podaci: (1/x,y)")

#Pearsonov koeficijent korelacije
r=(sum(x*y)-n*mean(x)*mean(y))/(sqrt(sum(x^2)-n*mean(x)^2)*sqrt(sum(y^2)-n*mean(y)^2))
r
r1=(sum(x*lny)-n*mean(x)*mean(lny))/(sqrt(sum(x^2)-n*mean(x)^2)*sqrt(sum(lny^2)-n*mean(lny)^2))
r1
r2=(sum(lnx*lny)-n*mean(lnx)*mean(lny))/(sqrt(sum(lnx^2)-n*mean(lnx)^2)*sqrt(sum(lny^2)-n*mean(lny)^2))
r2
r3=(sum(lnx*y)-n*mean(lnx)*mean(y))/(sqrt(sum(lnx^2)-n*mean(lnx)^2)*sqrt(sum(y^2)-n*mean(y)^2))
r3
r4=(sum(krozx*lny)-n*mean(krozx)*mean(lny))/(sqrt(sum(krozx^2)-n*mean(krozx)^2)*sqrt(sum(lny^2)-n*mean(lny)^2))
r4
r5=(sum(krozx*y)-n*mean(krozx)*mean(y))/(sqrt(sum(krozx^2)-n*mean(krozx)^2)*sqrt(sum(y^2)-n*mean(y)^2))
r5

#vidimo da je apsolutna vrijednost Pearsonovog koeficijenta korelacije najveci kod modela (2),(3) i (5)

#(b)

#uzeti æemo modele (2) i (5)
par(mfrow=c(1,1))

#prvo æemo naæi parametre za model (2)
X2=matrix(0,nrow=n,ncol=2)
X2[,1]=1
X2[,2]=lnx
X2

thetakapa_model2=solve(t(X2)%*%X2)%*%t(X2)%*%lny
thetakapa_model2

thetakapa0_model2=thetakapa_model2[1,]
thetakapa0_model2
thetakapa1_model2=thetakapa_model2[2,]
thetakapa1_model2

plot(lnx,lny,xlab="ln(x)",ylab="ln(y)",xlim=c(4.5,6.5),ylim=c(4,6),main="(lnx,lny) i pravac y'=7.363538-0.4405617x'")
abline(a=thetakapa0_model2,b=thetakapa1_model2,col="red")

#sada æemo naæi parametre za model (5)
X5=matrix(0,nrow=n,ncol=2)
X5[,1]=1
X5[,2]=krozx
X5

thetakapa_model5=solve(t(X5)%*%X5)%*%t(X5)%*%y
thetakapa_model5
thetakapa0_model5=thetakapa_model5[1,]
thetakapa0_model5
thetakapa1_model5=thetakapa_model5[2,]
thetakapa1_model5

plot(krozx,y,xlab="1/x",ylab="y",main="(1/x,y) i pravac y'=78.79752+13813.99841x'")
abline(a=thetakapa0_model5,b=thetakapa1_model5,col="red")

#(b)

#raèunamo statistiku R^2

#prvo nam trebaju procjenitelji za Yi

ykapica_model2=thetakapa0_model2+thetakapa1_model2*lnx
ykapica_model2

ykapica_model5=thetakapa0_model5+thetakapa1_model5*krozx
ykapica_model5

SSE2=sum((lny-ykapica_model2)^2)
SSE5=sum((y-ykapica_model5)^2)
SSE2
SSE5

SYY2=sum((lny-mean(lny))^2)
SYY5=sum((y-mean(y))^2)
SYY2
SYY5

R_2=1-SSE2/SYY2
R_2
R_5=1-SSE5/SYY5
R_5

#TEST ADEKVATNOSTI

#za model (2)
z_lnx=table(lnx)
z_lnx
n_lnx=c(3,3,4,4)
n_lnx
l=length(z_lnx)
l

j=1
z_lnx=numeric(l)
for( i in 1:l){
  z_lnx[i]=lnx[j]
  j=j+n_lnx[i]
}
z_lnx

yij=matrix(0,l,l)
j=1
for(i in 1:l)
{
  for(k in 1:n_lnx[i])
  {
    yij[i,k]=lny[j]
    j=j+1
  }
}

arit_yni=numeric(l)
for(i in 1:l)
{
  suma=0
  for(j in 1:n_lnx[i])
  {
    suma=suma+yij[i,j]
  }
  arit_yni[i]=suma/n_lnx[i]
}
arit_yni

SSEpe2=0
for(i in 1:l)
{
  for(j in 1:n_lnx[i])
  {
    SSEpe2=SSEpe2+(yij[i,j]-arit_yni[i])^2
  }
}
SSEpe2
k=1
F2=(SSE2-SSEpe2)*(n-l)/(SSEpe2*(l-k-1))
F2
pvalue=1-pf(F2,df1=l-k-1,df2=n-l)
pvalue
#hipotezu H_0 ne odbacujemo za alpha<=0.4758953

#za model (5)
z_krozx=table(krozx)
z_krozx
n_krozx=c(3,3,4,4)
l=length(z_krozx)
l

j=1
z_krozx=numeric(l)
for( i in 1:l){
  z_krozx[i]=krozx[j]
  j=j+n_krozx[i]
}
z_krozx


yijnovi=matrix(0,l,l)
j=1
for(i in 1:l)
{
  for(k in 1:n_krozx[i])
  {
    yijnovi[i,k]=y[j]
    j=j+1
  }
}
yijnovi

arit_yninovi=numeric(l)
for(i in 1:l)
{
  suma=0
  for(j in 1:n_krozx[i])
  {
    suma=suma+yijnovi[i,j]
  }
  arit_yninovi[i]=suma/n_krozx[i]
}
arit_yninovi

SSEpe5=0
for(i in 1:l)
{
  for(j in 1:n_krozx[i])
  {
    SSEpe5=SSEpe5+(yijnovi[i,j]-arit_yninovi[i])^2
  }
}
SSEpe5
k=1
F5=(SSE5-SSEpe5)*(n-l)/(SSEpe5*(l-k-1))
F5
pvalue=1-pf(F5,df1=l-k-1,df2=n-l)
pvalue

#hipotezu H0 ne odbacujemo za alpha<=0.5773824

#model (5) je bolji




#(c)
#--> izabrali smo model (5) za (c)
reziduali=y-ykapica_model5
reziduali
plot(ykapica_model5,reziduali,xlab="ykapica",ylab="reziduali",ylim=c(-32.5,35.2),xlim=c(100,205),main="graf reziduala za model (5)")
abline(0,0,col="blue")
k=1
sigmakapa_2=SSE5/(n-k-1)
sigmakapa=sqrt(sigmakapa_2)
sigmakapa
H=X5%*%solve(t(X5)%*%X5)%*%t(X5)
hii=diag(H)
hii
stand_rezid=reziduali/(sigmakapa*sqrt(1-hii))
stand_rezid
plot(ykapica_model5,stand_rezid,xlab="ykapica",ylab="standardizirani reziduali",xlim=c(100,205),ylim=c(-2,2),main="graf standardiziranih reziduala za model (5)")
abline(0,0,col="blue")

i=1:n
qi=qnorm((i-3/8)/(n+1/4))
qi
plot(qi,sort(stand_rezid),main="normalni vjerojatnosni graf za standardizirane reziduale",ylim=c(-2,2),ylab="standardizirani reziduali")
abline(a=0,b=1,col="red")

dn=0
sortstrez=sort(stand_rezid)
for(i in 1:n)
{
  prvi=abs((i-1)/n-pnorm(sortstrez[i]))
  drugi=abs(i/n-pnorm(sortstrez[i]))
  maxi=max(prvi,drugi)
  if(maxi>dn)
  {
    dn=maxi
  }
}
dn
ks.test(sortstrez,"pnorm",0,1)
#Za p<=0.8745=87.45% æemo odbaciti hipotezu
#da standardizirani reziduali dolaze iz
#jediniène normalne distribucije
#tj. za nivo znaèajnosti za koji vrijedi alpha<=0.8745 ne odbacujemo hipotezu

#(d)
alpha=0.05
t=qt(1-alpha/2,n-k-1)
t
C=solve(t(X5)%*%X5)
cjj=diag(C)
cjj
donja_granica=thetakapa_model5-t*sigmakapa*sqrt(cjj)
gornja_granica=thetakapa_model5+t*sigmakapa*sqrt(cjj)
procjena=matrix(0,nrow=2,ncol=2)
procjena[,1]=donja_granica
procjena[,2]=gornja_granica
procjena

#POUZDANO PODRUÈJE
f=qf(1-alpha,k+1,n-k-1)
f
B=t(X5)%*%X5
B
koef=f*(k+1)*sigmakapa_2
koef
a=B%*%thetakapa_model5
a
b=t(a)
b
c=t(thetakapa_model5)%*%a
c
#F*(k+1)*sigmakapa^2-koef<=0
elipsa=function(theta0,theta1)
{
  rbind(c(theta0,theta1))%*%B%*%t(rbind(c(theta0,theta1)))-rbind(c(theta0,theta1))%*%B%*%thetakapa_model5-t(thetakapa_model5)%*%B%*%t(rbind(c(theta0,theta1)))+t(thetakapa_model5)%*%B%*%thetakapa_model5-koef
  
}

library("car")
confidenceEllipse(lm(y~krozx),xlab="",ylab="",main="pouzdano podruèje", col="red",center.pch=F) 
#1.put ovo gore

#2.put ovo dolje
A=B/((k+1)*sigmakapa^2)
A
D=diag(eigen(A)$values)
D
V=eigen(A)$vectors
V
B_=-B%*%thetakapa_model5/(sigmakapa^2)
B_
C=t(thetakapa_model5)%*%B%*%thetakapa_model5/((k+1)*sigmakapa^2)-f
C
Bc=t(V)%*%B_
Bc

sx=-Bc[1]/(2*D[1,1])
sx
sy=-Bc[2]/(2*D[2,2])
sy
Cnovi=C-D[1,1]*sx^2-D[2,2]*sy^2
Cnovi
rx=sqrt(-Cnovi/D[1,1])
rx
ry=sqrt(-Cnovi/D[2,2])
ry
phi=seq(0,2*pi,length=50)
phi
xe=sx+rx%*%cos(phi)
xe
ye=sy+ry%*%sin(phi)
ye
xp=V[1,1]*xe+V[1,2]*ye
yp=V[2,1]*xe+V[2,2]*ye
plot(xp,yp, type="l",main="Pouzdano podruèje",col="red")

#3.put dolje
M=t(X5)%*%X5/((k+1)*sigmakapa^2*f)
M
install.packages("expm")
library(expm)
N = sqrtm(M)
N
xcrta=cos(phi)
xcrta
ycrta=sin(phi)
ycrta


#(e)
xkoord=seq(100,510,1)
xkoord
l=length(xkoord)
l
ykoord=thetakapa0_model5+thetakapa1_model5*(1/xkoord)
plot(x,y,xlab="težina tijela",ylab="odreðeno svojstvo metabolizma",ylim=c(90,200),xlim=c(100,550),main="95% p.i. za srednju vrijednost od Y (uz dano x)")
lines(xkoord,ykoord, col="red")

alpha=0.05
y_donja=numeric(l)
y_gornja=numeric(l)
X=matrix(0,nrow=n,ncol=2)
X[,1]=1
X[,2]=1/x
for(i in 1:l)
{
  z=c(1,1/xkoord[i])
  y_donja[i]=thetakapa0_model5+thetakapa1_model5*(1/xkoord[i])-t*sigmakapa*sqrt(z%*%solve(t(X)%*%X)%*%z)
  y_gornja[i]=thetakapa0_model5+thetakapa1_model5*(1/xkoord[i])+t*sigmakapa*sqrt(z%*%solve(t(X)%*%X)%*%z)
}
points(xkoord,y_donja,col="blue",type="l")
points(xkoord,y_gornja,col="blue",type="l")

plot(x,y,xlab="težina tijela",ylab="odreðeno svojstvo metabolizma",ylim=c(50,250),xlim=c(100,550),main="95% p.i. za vrijednost od Y (uz dano x)")
lines(xkoord,ykoord, col="red")
y_donjanovi=numeric(l)
y_gornjanovi=numeric(l)
for(i in 1:l)
{
  z=c(1,1/xkoord[i])
  y_donjanovi[i]=thetakapa0_model5+thetakapa1_model5*(1/xkoord[i])-t*sigmakapa*sqrt(1+z%*%solve(t(X)%*%X)%*%z)
  y_gornjanovi[i]=thetakapa0_model5+thetakapa1_model5*(1/xkoord[i])+t*sigmakapa*sqrt(1+z%*%solve(t(X)%*%X)%*%z)
}
points(xkoord,y_donjanovi,col="blue",type="l")
points(xkoord,y_gornjanovi,col="blue",type="l")


ykapica=thetakapa0_model5+thetakapa1_model5*krozx
plot(y,ykapica,main="(y,ykapica) s pravcem y=x")
abline(0,1,col="red")
#Model je valjda? dobar jer su podaci otprilike simetricno oko dijagonale
