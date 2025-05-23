fst<-read.csv2('../dados/fustes_med1_med2_spp.csv',fileEncoding = 'latin1')

#dmin=min(fst$dap1[fst$dap1>0])
dmin<-5
#dmax<-max(c(fst$dap1,fst$dap2))
dmax<-45
#nclasses<-nclass.Sturges(fst$dap2)
nclasses<-8

n_proj<-1 #Prognose a partir da segunda medição para n_proj * intervalo de tempo entre medições

fator_extrapolacao<-10000/1200; #(Área do povoamento)/(Área da amostra)

# Métodos: 
#  1: Razão de movimentação
#  2: Wahlenberg
#  3: Matriz de transição

prog_nativas<-function(dap1,dap2,dapmin=min(dap1[dap1>0]),dapmax=max(c(dap1,dap2)),
                       nclasses=nclass.Sturges(dap2), fator_extrapolacao=1, 
                       metodo=3, n_proj=1, breaks=NULL){
  if(is.null(breaks)){
    breaks=seq(dmin,dmax,len=nclasses+1)
  }else{
    breaks<-sort(breaks)
    dmin<-min(breaks)
    if(max(breaks)>=max(c(dap1,dap2))){
      dmax<-max(breaks)
    }else{
      dmax<-max(c(dap1,dap2))
      breaks[length(breaks)]<-dmax
    }
    nclasses<-length(breaks)-1
  } 
  fst<-data.frame(dap1=dap1,dap2=dap2)
  fst<-fst[(fst$dap1==0 | fst$dap1>=dmin) & (fst$dap2==0 | fst$dap2>=dmin),]
  fst$cld1<-cut(fst$dap1,breaks=breaks,include.lowest=T,right=F,labels=F)
  
  tfreq<-function(x,ndigits=2,...){
    hx<-hist(sort(x),...);
    if(!is.null(ndigits)) hx$mids<-round(hx$mids,ndigits)
    return(data.frame(vc=hx$mids,fo=hx$counts))
  }
  fo_recrut<-tfreq(fst$dap2[fst$dap1==0],breaks=breaks,include.lowest=T,right=F, plot=F)
  fo_mortas<-tfreq(fst$dap1[fst$dap2==0],breaks=breaks,include.lowest=T,right=F, plot=F)
  fo_med2<-tfreq(fst$dap2[fst$dap2>0],breaks=breaks,include.lowest=T,right=F,plot=F)
  
  fg<-function(cln,fst,breaks,metodo){
    sfst<-fst[fst$cld1==cln,]
    if(nrow(sfst)>0 & sum(sfst$dap2)>0){
      if(metodo==2 | metodo==3){
        prob=hist(sfst$dap2[sfst$dap2>0],breaks=breaks,include.lowest=T,right=F,plot=F)$counts/nrow(sfst)
      }else{
        rm<-mean(sfst$ip)/(breaks[cln+1]-breaks[cln])
        prob<-rep(0,length(breaks)-1)
        cls<-cln+trunc(rm);cls<-c(cls,cls+1)
        prob[cls[1]]<-1-(rm-trunc(rm))
        prob[cls[2]]<-1-prob[cls[1]]
      }
    }else{
      prob<-rep(0,length(breaks)-1)
      prob[cln]<-1
    }
    return(matrix(prob))
  }
  
  if(metodo==1 | metodo==2){
    fst$ip<-fst$dap2-fst$dap1
    g<-Reduce('cbind',lapply(1:nclasses,fg,fst[fst$dap1>0 & fst$dap2>0,],breaks,metodo))
    dimnames(g)<-list(fo_med2$vc,fo_med2$vc)
    
    fest<-fo_med2$fo
    for(i in 1:n_proj) fest<-g %*% matrix(fest) + matrix(fo_recrut$fo) - matrix(fo_mortas$fo)
    fest<-as.vector(fest)
    fest[fest<0]<-0
    prog<-data.frame(
      classe=cut(fo_med2$vc,breaks=breaks,include.lowest=T,right=F),
      fo_med2=fo_med2$fo*fator_extrapolacao,
      fo_recrutamento=fo_recrut$fo*fator_extrapolacao,
      fo_mortas=fo_mortas$fo*fator_extrapolacao,
      fo_esperada=fest*fator_extrapolacao
    )
  }else if(metodo==3){
    g<-Reduce('cbind',lapply(1:nclasses,fg,fst[fst$dap1>0,],breaks,metodo))
    dimnames(g)<-list(fo_med2$vc,fo_med2$vc)
    
    fest<-matrix(fo_med2$fo)
    for(i in 1:n_proj) fest<-g %*% fest + matrix(fo_recrut$fo)
    fest<-as.vector(fest)
    prog<-data.frame(
      classe=cut(fo_med2$vc,breaks=breaks,include.lowest=T,right=F),
      fo_med2=fo_med2$fo*fator_extrapolacao,
      fo_recrutamento=fo_recrut$fo*fator_extrapolacao,
      fo_esperada=fest*fator_extrapolacao
    )
  }
  return(list(prog=prog,g=g))
}    

# 1: Razão de movimentação
prog<-prog_nativas(fst$dap1,fst$dap2,dapmin,dapmax,nclasses,
                   fator_extrapolacao, metodo=1, n_proj)
prog$g
prog$prog

# 2: Wahlenberg
prog<-prog_nativas(fst$dap1,fst$dap2,dapmin,dapmax,nclasses,
                   fator_extrapolacao, metodo=2, n_proj)
prog$g
prog$prog

# 3: Matriz de transição
prog<-prog_nativas(fst$dap1,fst$dap2,dapmin,dapmax,nclasses,
                   fator_extrapolacao, metodo=3, n_proj)
prog$g
prog$prog

