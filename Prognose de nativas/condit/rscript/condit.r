arv<-read.csv2('../dados/fustes_med1_med2_xylopia.csv',fileEncoding = 'latin1')

arv$datamed1<-as.Date(arv$datamed1,format = "%d/%m/%Y")
arv$datamed2<-as.Date(arv$datamed2,format = "%d/%m/%Y")

#dmin=min(arv$dap1)
dapmin<-5
anos_prog<-10;

#cmrExam::disjoin(arv,vars=c('arv','dap1','dap2','datamed1','datamed2'))

condit<-function(arv,dap1,dap2,datamed1,datamed2,anos_prog, dapmin=min(dap1[dap1>0]), opcao=1, fact=5){
  arv<-data.frame(arv=arv,dap1=dap1,dap2=dap2, dt=round(as.vector(datamed2-datamed1)/365.25,2))
  arv<-arv[arv$dap1>=dmin & arv$dap2>=arv$dap1,]
  arv$g<-with(arv,(log(dap2)-log(dap1))/dt)

  if(opcao==1){
    coef_start<-as.vector(coef(lm(log(g)~log(dap1),data=arv,subset=g>0)))
    ajn<-nls('g~a*exp(b*log(dap1))',data=arv,start=list(a=coef_start[1],b=coef_start[2]))
    
    arv$ge<-predict(ajn)
    arv$ga<-as.vector(arv$g+abs(residuals(ajn)))
    aja<-nls('ga~a*exp(b*log(dap1))',data=arv,start=list(a=coef_start[1],b=coef_start[2]))
    arv$gae<-predict(aja)
    
    ft<-function(dap,a,b,m) -((1/(a*b))*exp(-b*log(dap)))+m
    
    cn<-as.vector(coef(ajn))
    mn<--ft(dapmin,cn[1],cn[2],0)
    arv$idade_prog_n<-ft(arv$dap2,cn[1],cn[2],mn)+anos_prog
    arv$dap_prog_n<-exp(log(cn[1]*cn[2]*(mn-arv$idade_prog_n))/-cn[2])

    ca<-as.vector(coef(aja))
    ma<--ft(dapmin,ca[1],ca[2],0)
    arv$idade_prog_a<-ft(arv$dap2,ca[1],ca[2],ma)+anos_prog
    arv$dap_prog_a<-exp(log(ca[1]*ca[2]*(ma-arv$idade_prog_a))/-ca[2])

  }else{
    coef_start<-as.vector(coef(lm(log(g)~dap1,data=arv,subset=g>0)))
    ajn<-nls('g~a*exp(b*dap1)',data=arv,start=list(a=coef_start[1],b=coef_start[2]))
    
    arv$ge<-predict(ajn)
    arv$ga<-as.vector(arv$g+abs(residuals(ajn)))
    aja<-nls('ga~a*exp(b*dap1)',data=arv,start=list(a=coef_start[1],b=coef_start[2]))
    arv$gae<-predict(aja)

    ftp<-function(dap,b,fact) sum(Reduce('c',lapply(1:fact,function(x)((-b*dap)^x)/(x*factorial(x)))))
    ft<-function(dap,a,b,m,fact) (1/a)*(log(dap)+ftp(dap,b,fact))-m
    
    #ft<-function(dap,a,b,m) (1/a)*(log(dap)+(-b*dap)+((-b*dap)^2)/(2*factorial(2))+((-b*dap)^3)/(3*factorial(3))-m)
    
    fo<-function(dap,a,b,m,idade,fact){abs(idade-ft(dap,a,b,m,fact))};
    
    fp<-function(x,g='g',ge='ge',cs,m,fact,anos_prog){
      idp<-ft(x$dap2,cs[1],cs[2],m,fact)+anos_prog
      g1<-x[[ge]]-abs(x[[ge]]-x[[g]])
      if(g1<0) g1<-0
      g2<-x[[ge]]+abs(x[[ge]]-x[[g]])
      ld1<-exp(g1*idp+log(x$dap2))
      ld2<-exp(g2*idp+log(x$dap2))
      ld1<-ifelse(ld1==ld2,ld1*0.9,ld1)
      #print(x$arv)
      return(optimize(fo,c(ld1,ld2),cs[1],cs[2],m,anos_prog,fact,maximum=F)$minimum)
    }
    cn<-as.vector(coef(ajn))
    mn<-ft(dapmin,cn[1],cn[2],0,fact)
    arv$idade_prog_n<-Reduce('c',lapply(arv$dap2,ft,cn[1],cn[2],mn,fact))+anos_prog
    arv$dap_prog_n<-Reduce('c',lapply(split(arv,1:nrow(arv)),fp,g='g',ge='ge',cs=cn,m=mn,fact,anos_prog=anos_prog))
    
    ca<-as.vector(coef(aja))
    ma<-ft(dapmin,ca[1],ca[2],0,fact)
    arv$idade_prog_a<-Reduce('c',lapply(arv$dap2,ft,ca[1],ca[2],ma,fact))+anos_prog
    arv$dap_prog_a<-Reduce('c',lapply(split(arv,1:nrow(arv)),fp,g='ga',ge='gae',cs=ca,m=ma,fact,anos_prog=anos_prog))
  }
  return(list(arv=arv,ajn=ajn,aja=aja))
}    

# opcao=1: Primeira modificação
prog1<-condit(arv$arv,arv$dap1,arv$dap2,arv$datamed1,arv$datamed2,anos_prog, dapmin, opcao=1, fact=5)
prog1$arv
prog1$ajn
prog1$aja

# opcao=2: Segunda modificação
prog2<-condit(arv$arv,arv$dap1,arv$dap2,arv$datamed1,arv$datamed2,anos_prog, dapmin, opcao=2, fact=5)
prog2$arv
prog2$ajn
prog2$aja

