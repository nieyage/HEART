setwd("D:/ArchR/fibroblast/")
library(ggplot2)
FB1<-c(4.18,5.80,9.97 ,5.33 ,2.85 ,5.41 ,12.08 ,7.64 ,19.17 ,14.90 ,3.79 ,4.75 ,5.67)
FB2<-c(19.47 ,5.98 ,5.60 ,4.46 ,16.71 ,7.36 ,18.29 ,4.39 ,6.82 ,10.63 ,25.65 ,5.23 ,11.47 ,7.55 ,7.65 ,4.38 ,5.85 ,8.27 ,4.87,7.01 )
FB3<-c(10.28 ,11.88 ,19.00 ,8.13 ,9.69 ,7.83 ,18.17 ,10.99 ,21.64 ,8.93 ,10.23 ,13.75 ,7.46 ,5.37 ,6.52 ,7.12 ,20.75 ,8.58 ,7.51 )
FB4<-c(11.14 ,7.96 ,4.41 ,6.94 ,16.62 ,13.76 ,6.10 ,7.40 ,7.99 ,8.43 ,6.84 ,2.62 ,5.79 ,6.67 ,8.00 ,22.45 ,8.30 ,5.59 ,5.09 ,7.41 ,9.92 ,4.77 ,15.40 ,7.69 ,6.58)

FB1<-ecdf(FB1)
FB2<-ecdf(FB2)
FB3<-ecdf(FB3)
FB4<-ecdf(FB4)
plot(FB1,col="#A20056B2")
lines(FB2, col = "#3B4992B2")
lines(FB3, col = "#F48639")
lines(FB4, col = "#D51F26")
lines(x, col = "black")
ks.test(FB3,x)
##random select some numbers for test###
set.seed(3)
x=runif(n=20,min = 0,max = 25)

####density plot####
plot(density(FB4), col = "#D51F26")

lines(density(FB1),col="#A20056B2")
lines(density(FB2), col = "#3B4992B2")
lines(density(FB3), col = "#F48639")
?plot


str(FB4)
subtype<-c(rep("FB1",14),rep("FB2",19),rep("FB3",20),rep("FB4",25))
value<-c(FB1,FB2,FB3,FB4)
den<-data.frame(Group=subtype,value=value)
den
p<-ggplot(den, aes(x = value))
p + geom_density(aes(color = Group))+scale_fill_manual(values = c("#A20056B2", "#3B4992B2","#F48639","#D51F26"))
