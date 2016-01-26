#Grid creation

gridPcp <- function(filled,points,sts,inidate,enddate,parallel=TRUE,ncpu=2){

  #matrix of distances
  ddsts=rbind(sts,points)
  distanc=dist(cbind(ddsts$X,ddsts$Y))/1000; distanc=as.matrix(distanc); colnames(distanc)=ddsts$ID; rownames(distanc)=ddsts$ID
  w= which(colnames(distanc)%in%as.character(sts$ID))
  distanc=distanc[,-c(w)]
  w= which(rownames(distanc)%in%as.character(points$ID))
  distanc=distanc[-c(w),]

  #vector of dates
  datess=seq.Date(inidate,enddate,by='day')

  predday=function(x,datess,filled,distanc=distanc,points=points,sts=sts){

    dw=which(x==datess)
    d=filled[dw,]
    pred=data.frame(matrix(NA,ncol=2,nrow=length(d))); pred[,1] <- names(d)
    names(pred)=c('ID','obs')
    pred$obs=as.numeric(d)
    predpoint=data.frame(matrix(NA,ncol=2,nrow=nrow(points))); predpoint[,1] <- points$ID
    names(predpoint)=c('ID','pred')
    #if all data is zero...
    if(max(pred$obs,na.rm=T)==0){
      predpoint$pred=0
    } else{#else...

      #point by point
      for(h in 1:nrow(predpoint)){
        kk=data.frame(ID=rownames(distanc),D=distanc[,which(predpoint$ID[h]==colnames(distanc))],
                      obs=pred$obs[which(rownames(distanc)%in%pred$ID)],
                      stringsAsFactors=F)
        kk=kk[order(kk$D),]
        wna=which(is.na(kk$obs))#Should not be missing data, it would introduce inhomogeneities
        if(length(wna)>0) kk=kk[-c(wna),]#just in case
        if(nrow(kk)<10) kk=kk else kk=kk[1:10,]

        if(max(kk$obs,na.rm=T)==0){predpoint$pred[h]=0} else{

          sts_can <- points[h,]
          sts_nns <- sts[which(sts$ID%in%kk$ID),]

          #binomial
          b <- kk$obs; b[b>0]=1
          DF <- data.frame(y=b,alt=sts_nns$ALT,lat=sts_nns$Y,lon=sts_nns$X)
          fmtb <- suppressWarnings(glm(y~alt+lat+lon, data=DF, family=binomial()))
          newdata=data.frame(alt=sts_can$ALT,lat=sts_can$Y,lon=sts_can$X)
          pb <- predict(fmtb,newdata=newdata,type='response')
          if(pb<0.001 & pb>0) pb <- 0.001
          pb <- round(pb,3)

          #data
          mini=min(kk$obs)/2
          maxi=max(kk$obs)+(max(kk$obs)-min(kk$obs))
          yr=as.numeric((kk$obs-mini)/maxi)
          DF$y=yr
          fmt <- suppressWarnings(glm(y~alt+lat+lon, data=DF, family=quasibinomial()))
          p <- predict(fmt,newdata=newdata,type='response')
          if(p<0.001 & p>0) p <- 0.001
          p <- round((p*maxi)+mini,3)

          #introducing data
          if(pb<0.5) predpoint$pred[h]=0
          if(pb>=0.5 & p<1) predpoint$pred[h]=1
          predpoint$pred[h]=p
        }#next pred calculation
      }#next station
    }
    dir.create('./gridded/',showWarnings = F)
    write.table(predpoint,paste('./gridded/',x,'.txt',sep=''),quote=F,row.names=F,sep='\t',na='')
  }

  #RUN gridPcp
  if(parallel){
    sfInit(parallel=TRUE,cpus=ncpu)
  }
  if(parallel){
    print('Creating daily files of grid')
    sfLapply(datess,fun=predday,filled=filled,datess=datess,distanc=distanc,points=points,sts=sts)
  } else {
    lapply(datess,FUN=predday,filled=filled,datess=datess,distanc=distanc,points=points,sts=sts)
  }
  gc()
  sfStop()

}
