#***************************************************
#   Function to calculate partial linear model
#***************************************************
.First.lib <- function(lib, pkg) {
    library.dynam("twslm", pkg, lib)
}

library(splines)                            #load bs package

#****************************************************************
#             Constant variance model                           *
#****************************************************************
non.robust.twslm<-function(sld,geneid,rt,intn,df=12,degree=3,norm.only=TRUE,
   tol=1e-5){

     cat("Calculation is running..."," ",date(),"\n")
     EPS<-tol
     id<-as.integer(unclass(factor(geneid)))     #make id from 1,2,....J
     slide<-as.integer(unclass(factor(sld)))     #make slide 1,2,...,n
     ii<-order(slide,id)

     id<-id[ii]
     slide<-slide[ii]
     ratio<-rt[ii]
     intensity<-intn[ii]

     uniqueid<-unique(geneid[ii])             #unique geneid vector
     uniqueslide<-sld[ii]

     geneno<-length(uniqueid)        #total gene number
     obsno<-length(id)               #number of row
     slidno<-length(unique(slide))   #slide number
     Nj<-as.vector(table(id))        #number of each clones in the experiment
     Ni<-as.vector(table(slide))     #number of clones in each slide
   #initial beta as the mean of ratio
    beta0<-beta<-.C("ymean",
        y=as.double(ratio),
      gno=as.integer(geneno),
       id=as.integer(id),
      obs=as.integer(obsno),
       nj=as.integer(Nj),
    ymean=double(geneno),PACKAGE="twslm")$ymean
    beta0[1]<--sum(beta0[2:length(beta0)])

    Bfit<-NULL
    fittedvalue<-rep(0,obsno)
    residi<-resid0i<-sqrt(sum(ratio^2))     #initialize resid0 and resid
    eps<-iter<-1

while(eps>EPS){
  Bfit<-NULL
  #B-spline by slide
  for(i in unique(slide)){
     betaresid<-.Call("resid", as.double(ratio[slide==i]),
                               as.integer(id[slide==i]),
                               as.double(beta),PACKAGE="twslm")
      b<-intensity[slide==i]
      tmp<-fitted(lm(betaresid~bs(b,df=df,degree=degree)))
      Bfit<-c(Bfit,as.vector(tmp))
   }

   fittedvalue<-.Call("fittedvalue",as.double(Bfit),as.double(beta),as.integer(id),PACKAGE="twslm")
   residi<-sqrt(sum((ratio-fittedvalue)^2))
   eps<-abs(residi-resid0i)/resid0i
   resid0i<-residi
   iter<-iter+1
   beta0<-beta

   #calculate beta value
   beta<-.Call("UAnovamean",y=as.double(ratio-Bfit),
                          nj=as.double(Nj*1.0),
                          gno=as.integer(geneno),
                          id=as.integer(id),PACKAGE="twslm")
 }
#****  end of parameter estimation  *****
  beta0<-as.matrix(beta0)
  rownames(beta0)<-uniqueid
  colnames(beta0)<-"beta"
  betamean<-.Call("Wymean",as.double(ratio-Bfit),
                          as.integer(id),
                          as.integer(geneno),
                          as.double(Nj*1.0),PACKAGE="twslm");
  betamean<-as.matrix(betamean)
  rownames(betamean)<-uniqueid
  colnames(betamean)<-"ymean"

  parm<-list(name=uniqueid,beta=beta0,ymean=betamean,bvar=NULL,
             fittedvalue=fittedvalue, bfit=Bfit, slide=NULL,id=geneid[ii],ratio=NULL,
             intensity=NULL,scale=NULL,rscale=NULL)

  if(norm.only==TRUE){
    parm$slide<-uniqueslide
    parm$ratio<-ratio
    parm$intensity<-intensity
    cat("Calculation is end. ",date(),"\n")
    return (parm)
  }
  else{

#****Begin to calculate covariance matrix of betahat********
#construct B spline matrix

  B<-Bmean<-NULL
  tmprow<-0;count<-0
  for(i in unique(slide)){
       b<-intensity[slide==i]
       count<-count+1
       tmp<-bs(b,df=df,degree=degree)[1:length(b),1:df]
       tmp<-cbind(1,tmp)
       tmprow<-tmprow+length(b)
        if(count==1)   B<-cbind(B,rbind(tmp,matrix(0,nrow=obsno-length(b),ncol=df+1)))
        else if(count==slidno) B<-cbind(B,rbind(matrix(0,nrow=obsno-length(b),ncol=df+1),tmp))
        else B<-cbind(B,rbind(matrix(0,nrow=tmprow-length(b),ncol=df+1),
                  tmp,matrix(0,nrow=obsno-tmprow,ncol=df+1)))
    }
   tmpm<-.C("UAnovabmean",b=as.double(B),
             nr=as.integer(nrow(B)),
             nc=as.integer(ncol(B)),
            gid=as.integer(id),
            gno=as.integer(geneno),
             nj=as.double(Nj*1.0),
          bmean=double(geneno*ncol(B)),PACKAGE="twslm")$bmean
     Bmean<-matrix(tmpm,nrow=geneno,ncol=ncol(B))
     rm(tmpm)
  Bresid<-.Call("Bresid",b=B,bmean=Bmean,id=as.integer(id),PACKAGE="twslm")
  rm(B)
  Finverse<-solve(t(Bresid)%*%Bresid)                     #F^(-1)
  rm(Bresid)
  Bmean.Finverse<-Bmean[-1,]%*%Finverse                   #(Z'Z)^(-1)Z'BF^(-1)
  rm(Finverse)
  ss2<-residi^2/(obsno-slidno*(df+1)-geneno+1)  #S^2
   var.beta<-.Call("betavar",bf=Bmean.Finverse,
                  b=t(Bmean[-1,]),
                  nj=as.double(Nj*1.0),
                  gno=as.integer(geneno),PACKAGE="twslm")*ss2  #variance of betahat
    rm(Bmean.Finverse,Bmean)
    parm$bvar<-var.beta
    parm$slide<-uniqueslide
    parm$ratio<-ratio
    parm$intensity<-intensity
    parm$scale<-ss2
    rm(var.beta,slide,ratio,intensity,uniqueslide,Nj)

    cat("Calculation is end. ",date(),"\n")

    return(parm)
   }
}

#******************************************************************************
#                          Robust Model                                       *
#******************************************************************************
robust.twslm<-function(sld,geneid,rt,intn,df=12,degree=3, norm.only=TRUE,
      robust=TRUE,robust.name="Tukey",scale.constant=2.5,weight.constant=4.685,ibeta=NULL,iscale=NULL,tol=1e-5){

     cat("Calculation is running..."," ",date(),"\n")
     EPS<-tol
     id<-as.integer(unclass(factor(geneid)))     #make id from 1,2,....J
     slide<-as.integer(unclass(factor(sld)))     #make slide 1,2,...,n
     ii<-order(slide,id)

     id<-id[ii]
     slide<-slide[ii]
     ratio<-rt[ii]
     intensity<-intn[ii]
     uniqueid<-unique(geneid[ii])    #unique geneid vector
     uniqueslide<-sld[ii]

     geneno<-length(uniqueid)        #total gene number
     obsno<-length(id)               #number of row
     slidno<-length(unique(slide))   #slide number
     Nj<-as.vector(table(id))        #number of each clones in the experiment
     Ni<-as.vector(table(slide))     #number of clones in each slid
     beta<-beta0<-NULL                               #initial beta as the mean of ratio
    if(is.null(ibeta)==TRUE){
        beta0<-.C("ymean",
              as.double(ratio),
              as.integer(geneno),
              as.integer(id),
              as.integer(obsno),
              as.integer(Nj),
              double(geneno),PACKAGE="twslm")[[6]]
    beta0[1]<--sum(beta0[2:length(beta0)])}
    else beta0<-ibeta  

    beta<-beta0
    Bfit<-NULL
    fittedvalue<-rep(0,obsno)
    residi<-resid0i<-sqrt(sum(ratio^2))     #initialize resid0 and resid
    resido0<-residi
     w0<-w<-rep(1,obsno)                                    #initial weights as 1s
    eps<-1;iter<-1;iter1<-1;eps1<-1;
    sumw0<-sumw<-NULL
    isigma<-NULL 
    if(is.null(iscale)==TRUE) 
         isigma<-1.482602*median(abs(ratio-median(ratio)))               #initialize sigma
    else isigma<-iscale
    fsigma<-isigma
while(eps>EPS){

    Bfit<-NULL
    for(i in unique(slide)){
        betaresid<-.Call("resid", as.double(ratio[slide==i]),
                                  as.integer(id[slide==i]),
                                  as.double(beta),PACKAGE="twslm")
               b<-intensity[slide==i]
             tmp<-fitted.values(lm(betaresid~bs(b,df=df,degree=degree),weights=w[slide==i]))
            Bfit<-c(Bfit,as.vector(tmp))
      }
      fittedvalue<-.Call("fittedvalue",as.double(Bfit),as.double(beta),as.integer(id),PACKAGE="twslm")
      beta0<-beta
      residual<-ratio-fittedvalue

      if(robust.name=="Huber")
        residi=.Call("Robustobj",
                      as.double(residual),
                      as.double(isigma),
                      as.double(weight.constant),
                  as.integer(1),PACKAGE="twslm")                              #Huber's objective function
      else if(robust.name=="Tukey")
        residi<-.Call("Robustobj",
                  as.double(residual),
                  as.double(isigma),
                  as.double(weight.constant),
                  as.integer(2),PACKAGE="twslm")
      eps<-abs(residi-resid0i)/resid0i
      resid0i<-residi
      iter<-iter+1
      isigma<-fsigma               #assign fisigma to isigma
      sumw0<-sumw
      w0<-w
    if (is.na(isigma)) isigma=1
    #Robust regression
      if(robust==TRUE){
        if(robust.name=="Tukey"){
          fsigma<-.Call("Tukeyscale",
                    as.double(residual),
                        as.integer((df+1)*slidno+geneno-1),
                        as.double(isigma),
                        as.double(scale.constant),PACKAGE="twslm")
      w<-.Call("TukeyWeight",as.double(residual),as.double(isigma),as.double(weight.constant),PACKAGE="twslm")
        }
        else if (robust.name=="Huber"){
            fsigma<-.Call("Huberscale",
                         as.double(residual),
                         as.integer((df+1)*slidno+geneno-1),
                         as.double(isigma),
                         as.double(scale.constant),PACKAGE="twslm")
             fsigma<-sqrt(fsigma)

      w<-.C("Huber",
          as.double(residual),
              double(obsno),
              as.integer(obsno),
              as.double(isigma),
              as.double(weight.constant),PACKAGE="twslm")[[2]]
        }
       }
       sumw<-.Call("sumweight",
                 as.double(w),
                 as.integer(id),
                 as.integer(geneno),PACKAGE="twslm")  #sum of weight for each gene

      beta<-.Call("UAnovamean",
                 as.double(w*(ratio-Bfit)),
                 as.double(sumw),
                 as.integer(geneno),
                 as.integer(id),PACKAGE="twslm")
    }

#save the results in
  beta0<-as.matrix(beta0)
  rownames(beta0)<-uniqueid
  colnames(beta0)<-"beta"

  betamean<-.Call("Wymean",as.double(w0*(ratio-Bfit)),
                          as.integer(id),
                          as.integer(geneno),
                          as.double(sumw0),PACKAGE="twslm")
  betamean<-as.matrix(betamean)
  rownames(betamean)<-uniqueid
  colnames(betamean)<-"ymean"

  parm<-list(name=uniqueid,beta=beta0,ymean=betamean,bvar=NULL,
             fittedvalue=fittedvalue,bfit=Bfit,slide=NULL,id=geneid[ii],ratio=NULL,
             intensity=NULL,scale=isigma^2,rscale=NULL)
  if(norm.only==TRUE){
    parm$slide<-uniqueslide
    parm$ratio<-ratio
    parm$intensity<-intensity
    cat("Calculation is end. ",date(),"\n")
    return (parm)
  }
else{
#*****calculate the robust variance estimator of beta*****
   #construct B spline matrix****
    B<-Bmean<-NULL
    tmprow<-0;count<-0
    for(i in unique(slide)){
      count<-count+1
      b<-intensity[slide==i]
      tmp<-bs(b,df=df,degree=degree)[1:length(b),1:df]
      tmp<-cbind(1,tmp)
      tmpm<-.C("UAnovabmean",
              as.double(tmp),
              as.integer(nrow(tmp)),
              as.integer(ncol(tmp)),
              as.integer(id[slide==i]),
              as.integer(geneno),
              as.double(Nj*1.0),
              double(geneno*ncol(tmp)),PACKAGE="twslm")[[7]]
       Bmean<-cbind(Bmean,matrix(tmpm,nrow=geneno,ncol=ncol(tmp)))      ##Bmeans
       tmprow<-tmprow+length(b)
     if(count==1)    B<-cbind(B,rbind(tmp,matrix(0,nrow=obsno-length(b),ncol=df+1)))
     else if(count==slidno) B<-cbind(B,rbind(matrix(0,nrow=obsno-length(b),ncol=df+1),tmp))
     else B<-cbind(B,rbind(matrix(0,nrow=tmprow-length(b),ncol=df+1),
                        tmp,matrix(0,nrow=obsno-tmprow,ncol=df+1)))
    }
  #call bmean function to calculate B(1)
    Bresid<-.Call("Bresid",b=B,bmean=Bmean,id=as.integer(id),PACKAGE="twslm")
   rm(B)
    Finverse<-solve(t(Bresid)%*%Bresid)               #F^(-1)
    Bmean.Finverse<-Bmean[-1,]%*%Finverse             #(Z'Z)^(-1)Z'BF^(-1)
    ss2<-NULL               #calculate Robust variance
    if(substr(robust.name,1,5)=="Huber") {
      ss2<-.Call("Robustvar",
                as.double(ratio-parm$fittedvalue),
                as.integer(slidno*(df+1)-geneno+1),
            as.double(weight.constant),
            as.double(isigma),
            as.integer(1),PACKAGE="twslm")
    }
    else if(substr(robust.name,1,5)=="Tukey")  {
      ss2<-.Call("Robustvar",
                 as.double(ratio-parm$fittedvalue),
                 as.integer(slidno*(df+1)-geneno+1),
             as.double(weight.constant),
             as.double(isigma),
             as.integer(2),PACKAGE="twslm")           #calculate Robust variance
    }
    var.beta<-.Call("betavar",
                     Bmean.Finverse,
                     t(Bmean[-1,]),
                 as.double(Nj*1.0),
                 as.integer(geneno),PACKAGE="twslm")*ss2*isigma^2  #variance of betahat
                                        #cat("ss2=",ss2,"\n")

    parm$bvar<-var.beta
    parm$slide<-uniqueslide;
    parm$ratio<-ratio;
    parm$intensity<-intensity
    parm$rscale<-ss2
    ww2<-.Call("w2",as.double(w0),as.integer(id),as.integer(geneno),PACKAGE="twslm")
    rm(ratio,slide,intensity,id,var.beta,Ni,Nj,uniqueslide,w0,sumw0)

    cat("Calculation is end. ",date(),"\n")

   return(parm)
 }
}

#************************************
#   Blockwise Normalization Function*
#************************************

BlockByBlock<-function(sld,blk,geneid,rt,intn,df=12,degree=3,norm.only=TRUE,
  robust=TRUE,robust.name="Tukey",scale.constant=2.5, weight.constant=4.685,ibeta=NULL,iscale=NULL,tol=1e-5){

 beta<-NULL;bvar<-NULL;ymean<-NULL; slide<-NULL; block<-NULL;ID<-NULL
 ratio<-NULL;intensity<-NULL; bfit<-NULL; fittedvalue<-NULL;scale<-NULL;name<-NULL;
 rscale<-NULL
 cat("Block Number=")
 blockid<-unique(cbind(blk,geneid))
 for(i in unique(blk)){
  cat(i," ")
   as<-sld[blk==i]
   ai<-geneid[blk==i]
   ar<-rt[blk==i]
   ain<-intn[blk==i]
   bi<-rep(i,length(as))
   a<-NULL
   ii<-(blockid[,1]==i)
   if(robust==TRUE)
     a<-robust.twslm(sld=as,geneid=ai,rt=ar,intn=ain,df=df,degree=degree,
     norm.only=norm.only,robust=robust,robust.name=robust.name, scale.constant=scale.constant,
     weight.constant=weight.constant,ibeta=ibeta[ii==TRUE],iscale=iscale[i],tol=tol)
  else a<-non.robust.twslm(sld=as,geneid=ai,rt=ar,intn=ain,df=df,degree=degree,
          norm.only=norm.only,tol=tol)
   beta<-c(beta,a$beta)
   bvar<-c(bvar,a$bvar)
   ymean<-c(ymean,a$ymean)
   slide<-c(slide,a$slide)
   block<-c(block,bi)
   ID<-c(ID,a$id)
   ratio<-c(ratio,a$ratio)
   intensity<-c(intensity,a$intensity)
   bfit<-c(bfit,a$bfit)
   fittedvalue<-c(fittedvalue,a$fittedvalue)
   scale<-c(scale,a$scale)
   rscale<-c(rscale,a$rscale)
   name<-c(name,a$name)
 }
   beta<-as.matrix(beta)
   rownames(beta)<-name
   colnames(beta)<-"beta"
   ymean<-as.matrix(ymean)
   rownames(ymean)<-name
   colnames(ymean)<-"ymean"

 result<-list(name=name,beta=beta,bvar=bvar,ymean=ymean,bfit=bfit,
             fittedvalue=fittedvalue,ratio=ratio,intensity=intensity,
             slide=slide,id=ID,block=block,scale=scale,rscale=rscale)
 return(result)
}


#***************************************************
# Semiparametric model normalization main function *
#***************************************************

twslm<-function(sld,blk,geneid,rt,intn,df=12,degree=3,norm.only=TRUE, block.norm=FALSE,
    robust=TRUE,robust.name="Tukey", scale.constant=2.5,weight.constant=4.685,ibeta=NULL,iscale=NULL,tol=1e-5)
{
  param<-NULL
  if(block.norm==FALSE){
    if(robust==FALSE)
        param<-non.robust.twslm(sld=sld,geneid=geneid,rt=rt,intn=intn,df=df,degree=degree,
        norm.only=norm.only,tol=tol)
    else param<-robust.twslm(sld=sld,geneid=geneid,rt=rt,intn=intn,df=df,degree=degree,
    norm.only=norm.only,robust=robust,robust.name=robust.name,
    scale.constant=scale.constant,weight.constant=weight.constant,ibeta=ibeta,iscale=iscale,tol=tol)
  }
  else param<-BlockByBlock(sld=sld,blk=blk,geneid=geneid,rt=rt,intn=intn,df=df,degree=degree,
       norm.only=norm.only,robust=robust,robust.name=robust.name,scale.constant=scale.constant,
       weight.constant=weight.constant,ibeta=ibeta,iscale=iscale,tol=tol)

  return(param)
}
