
arv<-read.csv2('../dados/arvores.csv',stringsAsFactors=F);
View(arv);
arv$dap<-arv$cap/pi;
arv$cldap<-cut(arv$dap,seq(5,33,4),right=F);

graphics.off(); x11()
hist(arv$dap,breaks=seq(5,33,4))

tf<-table(arv$cldap);
tf<-data.frame(cldap=names(tf),fo=as.vector(tf));
tf$fo<-round(tf$fo*10000/1800,2);
tf$vc<-seq(7,31,4);
tf<-tf[,c(1,3,2)];

tf$lnfo<-log(tf$fo);
View(tf)

aj_lin<-lm('log(fo)~vc',tf);
coefs_lin<-as.vector(coef(aj_lin));

numerador<-sum(tf$vc*tf$lnfo)-sum(tf$vc)*sum(tf$lnfo)/nrow(tf)
denominador<-sum(tf$vc^2)-sum(tf$vc)^2/nrow(tf)
b1<-numerador/denominador
b0<-mean(tf$lnfo)-b1*mean(tf$vc)

aj_nlin<-nls('fo~b0*exp(b1*vc)',tf,start=list(b0=exp(7.54),b1=-0.18));
coefs_nlin<-as.vector(coef(aj_nlin))

log(coefs_nlin[1])

#tf$fe<-predict(aj_nlin);

tf$fe<-exp(predict(aj_lin));

nc<-nrow(tf)
#Quociente de De Liocourt
tf$fe[1:(nc-1)]/tf$fe[2:nc]

with(tf,plot(vc,fo,type='l',col=2,xlab='dap(cm)',ylab='frequência'));
#with(tf,lines(vc,fe,col=3));

fpred<-function(x, ps) exp(ps[1]+ps[2]*x)
curve(
  expr = fpred(x, coefs_lin),
  from = min(tf$vc),
  to = max(tf$vc),
  add=T,
  col=3
)
legend(
  'topright', 
  legend=c('F. observada', 'F. estimada'), 
  text.col=c(2,3),
  horiz = T,
  bty = 'n'
)

b1
q<-tf$fe[1]/tf$fe[2]
b1<-log(q)/-4

View(tf)

tf$gi<-pi*tf$vc^2/40000

tf$Gi<-tf$gi*tf$fo

Gorig<-sum(tf$Gi) #m²/ha

Grem<-Gorig*0.5  #m²/ha

b0_original<-exp(coefs_lin[1])
b1_original<-coefs_lin[2]

b1<-b1_original;

b0<-with(tf,(40000*Grem)/(pi*sum(vc^2*exp(b1*vc))))
log(b0);

tf$freman<-b0*exp(b1*tf$vc);

with(tf,plot(vc,fo,type='l',col=2,xlab='dap(cm)',ylab='frequência'));
#with(tf,lines(vc,freman,col=3));
curve(
  expr = fpred(x, c(log(b0),b1)),
  from = min(tf$vc),
  to = max(tf$vc),
  add=T,
  col=3
)
legend(
  'topright', 
  legend=c('F. observada', 'F. remanescente'), 
  text.col=c(2,3),
  horiz = T,
  bty = 'n'
)

tf$fremov<-with(tf,fo-freman)


##Novo plano de manejo
#Novo q: 20% superior ao q original
#G remanescente: 50% do G original
q_original<-tf$fe[1]/tf$fe[2]

nq<-q_original*1.2
b1<-log(nq)/-4

Grem<-Gorig*0.5
b0<-with(tf,(40000*Grem)/(pi*sum(vc^2*exp(b1*vc))))
log(b0);

tf$freman<-b0*exp(b1*tf$vc);
tf$fremov<-with(tf,fo-freman)

with(tf,plot(vc,fo,type='l',col=2,xlab='dap(cm)',ylab='frequência'));
#with(tf,lines(vc,freman,col=3));
curve(
  expr = fpred(x, c(log(b0),b1)),
  from = min(tf$vc),
  to = max(tf$vc),
  add=T,
  col=3
)
legend(
  'topright', 
  legend=c('F. observada', 'F. remanescente'), 
  text.col=c(2,3),
  horiz = T,
  bty = 'n'
)

View(tf)