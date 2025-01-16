#zadanie 1
x = matrix(c(4,0,3,17),nrow=2,dimnames = list(Uszkodzenie=c("Tak","Nie"),Temperatura=c("<65", ">65")))
x
ftable(x)
fisher.test(ftable(x))
fisher.test(x[,1], x[,2]) #wektory 0,10,101,10,
prop.test(x[,1],rowSums(x), alternative = "t")


#zadanie2
reakcja <- read.csv2(file="Reakcja.csv",header=TRUE)
ftable(reakcja)
#2a
table2a = ftable(reakcja,col.vars = "Reakcja", row.vars="Dawka")
fisher.test(table2a)
table2a
#2b #prop test
table2b = ftable(reakcja,col.vars = "Reakcja", row.vars="Rodzaj")
table2b
fisher.test(table2b)
prop.test(table2b)
#2c
table2c = ftable(reakcja,col.vars = "Reakcja", row.vars="Miejsce")
table2c
fisher.test(table2c)
prop.test(table2c)


#zadanie3
tabela1 = matrix(c(32, 44, 60, 70, 22, 38 ,104 ,125, 13,48 ,61 ,113 , 3 ,18 ,54 ,96),nrow=4,byrow=TRUE)
dimnames(tabela1) <- list(Wynagrodzenie=c("ponizej 6000 ", "6000-15000 ", "15000-25000" , "powyzej 25000"),Zadowoleniezpracy=c("b. niezadow.", "niezadow.", "zadow." ,"b. zadow."))

chisq.test(tabela1)
chisq.test(tabela1,correct=FALSE)
chisq.test(tabela1,simulate.p.value = TRUE)


library(vcd)
assoc(tabela1,shade=TRUE)



#zadanie4
test_iloraz_wiarygodnosci = function(dane){
  g2=1
  for (i in 1:dim(dane)[1]){
    for (j in 1:dim(dane)[2]){
      a = ((sum(dane[i,])*sum(dane[,j]))/(sum(tabela1)*dane[i,j]))^dane[i,j]
      g2 = a*g2
}
  }
  G2 =  -2*log(g2)
  p = 1-pchisq(G2,prod((dim(dane)-1)))
  return (p)
}
library(vcd)
assocstats(tabela1)
test_iloraz_wiarygodnosci(tabela1)


#zadanie5
library(DescTools)
reakcja <- read.csv2(file="Reakcja.csv",header=TRUE)
ftable(reakcja)

GoodmanKruskalGamma(ftable(reakcja,col.vars="Reakcja",row.vars = "Dawka"),direction = "column")
ftable(reakcja,col.vars="Dawka",row.vars = "Reakcja")

GoodmanKruskalTau(ftable(reakcja,col.vars="Reakcja",row.vars = "Miejsce"),direction = "column")
ftable(reakcja,col.vars="Miejsce",row.vars = "Reakcja")



#zadanie 6
GoodmanKruskalGamma(ftable(tabela1,col.vars = "Wynagrodzenie",row.vars = "Zadowoleniezpracy"),direction="columnn")
GoodmanKruskalGamma(tabela1)
#zadanie7
leki = matrix(c( 35, 0, 0, 22, 22, 0, 15, 15, 15, 0, 40, 10, 18, 3, 5 ),nrow=5,byrow=TRUE)
dimnames(leki) <- list(Lek=c("Ibuprom","Apap","Paracetamol","Ibuprofen","Panadol"),Wiek=c("do lat 35" ,"od 36 do 55" ,"powyzej 55"))
leki
ftable(addmargins(leki))
     
GoodmanKruskalTau(ftable(leki,col.vars = "Lek",row.vars = "Wiek"),direction="column")
GoodmanKruskalTau(leki,direction="column")
#zadanie7

#analiza korespondencji  
leki[,1]
leki/sum(leki)
P = function(dane) {
  return (dane/sum(dane))
}
r = function(dane){
  return ((addmargins(dane)[nrow(dane)+1,]/sum(dane))[-ncol(dane)-1])
}
c = function(dane){
  return ((addmargins(dane)[,ncol(dane)+1]/sum(dane))[-nrow(dane)-1])
}

R = function(dane){
  return (inv(diag(r(dane)))*P(dane))
}

C = function(dane){
  return (P(dane)%*%inv(diag(c(dane))))
}

A = function(dane){
  return ((solve (diag(r(dane))^(1/2)) )%*%(t(P(dane))-r(dane) %*% t(c(dane))) %*% (solve (diag(c(dane))^(1/2)) ))
}

FG = function(dane){
  a = A(dane)
  u = svd(dane)$u
  v = svd(dane)$v
  GAM = svd(dane)$d
  F = t((solve(diag(r(dane))^(1/2)) )%*%t(u))%*%GAM
  G = (solve (diag(c(dane))^(1/2)) )%*%t(v)%*%GAM
  return (list("F"=F,"G"=G))
}
FG(leki)
A(leki)
svd(dane)


expandMatrix <- function(X, nrow, ncol) {
  X <- cbind(X, matrix(0, nrow = nrow(X), ncol = ncol - ncol(X)))
  X <- rbind(X, matrix(0, nrow = nrow - nrow(X), ncol = ncol(X)))
  X
}
ca(dane)


library(ca)
plot(ca(dane))

FG = function(dane){
  c = t(((addmargins(dane)[nrow(dane)+1,]/sum(dane))[-ncol(dane)-1]))
  r = t(((addmargins(dane)[,ncol(dane)+1]/sum(dane))[-nrow(dane)-1]))
  Dr = diag(c(r))
  Dc = diag(c(c))
  P = dane/sum(dane)
  A = solve(Dr^(1/2)) %*% (P-t(r) %*% c) %*% solve(Dc^(1/2))
  A_rozklad = svd(A)
  U = A_rozklad$u
  V = A_rozklad$v
  D = diag(A_rozklad$d)
  F = solve(Dr^(1/2)) %*% U %*% D
  G = solve(Dc^(1/2)) %*% V %*% D
  inercja = sum(D^2)
  inercje_procent = A_rozklad$d^2/inercja
  return (list("F"=F,"G"=G,"inercja"=inercja,"procent"=inercje_procent))
}
wsF = FG(leki)$F
wsG = FG(leki)$G
inercja = FG(leki)$inercja
inercja_wymiarow= FG(leki)$procent
xlimit = c(min(min(wsF[,1]),min(wsG[,1])), max(max(wsF[,1]),max(wsG[,1])))
ylimit = xlim = c(min(min(wsF[,2]),min(wsG[,2])), max(max(wsF[,2]),max(wsG[,2])))
xlabel = paste("Wymiar 1-inercja",toString(round(inercja_wymiarow[1]*100,1)),"%")
ylabel = paste("Wymiar 1-inercja",toString(round(inercja_wymiarow[2]*100,1)),"%")
plot(wsF, col=c("blue"),pch=16,xlim=xlimit,ylim=ylimit,xlab=xlabel,ylab=ylabel)
abline(h=0,v=0,lty=2, lwd=1)
points(wsG, col=c("red"),pch=15)
leki_nazwy = c("Ibuprom","Apap","Paracetamol","Ibuprofen","Panadol")
wiek_nazwy = c("do lat 35" ,"od 36 do 55" ,"powyzej 55")
for (i in 1:nrow(wsF)){
  text(x=wsF[i,1]+0.07,y=wsF[i,2]+0.02,leki_nazwy[i])
}
for (i in 1:nrow(wsG)){
  text(x=wsG[i,1]+0.07,y=wsG[i,2]+0.02,wiek_nazwy[i])
}
title("Analiza korespondencji dla zadania 7")


wsF = FG(tabela1)$F
wsG = FG(tabela1)$G
inercja = FG(tabela1)$inercja
inercja_wymiarow= FG(tabela1)$procent
xlimit = c(min(min(wsF[,1]),min(wsG[,1])), max(max(wsF[,1]),max(wsG[,1])))
ylimit = xlim = c(min(min(wsF[,2]),min(wsG[,2])), max(max(wsF[,2]),max(wsG[,2])))
xlabel = paste("Wymiar 1-inercja",toString(round(inercja_wymiarow[1]*100,1)),"%")
ylabel = paste("Wymiar 1-inercja",toString(round(inercja_wymiarow[2]*100,1)),"%")
plot(wsF, col=c("blue"),pch=16,xlim=xlimit,ylim=ylimit,xlab=xlabel,ylab=ylabel)
abline(h=0,v=0,lty=2, lwd=1)
points(wsG, col=c("red"),pch=15)
zarobki_nazwy = c("ponizej 6000 ", "6000-15000 ", "15000-25000" , "powyzej 25000")
zadow_nazwy = c("b. niezadow.", "niezadow.", "zadow." ,"b. zadow.")
for (i in 1:nrow(wsF)){
  text(x=wsF[i,1]+0.01,y=wsF[i,2]+0.01,zarobki_nazwy[i])
}
for (i in 1:nrow(wsG)){
  text(x=wsG[i,1]+0.01,y=wsG[i,2]+0.01,zadow_nazwy[i])
}
inercja


plot(ca(leki))
