# UCM data processing
  #' Usage: load(subjects=c("Paula","Angelo"),wd=c("/Users/angelobj/Box Sync/ASA old paper/Prueba/"),Task="ASA")
  #' @subjects: Name of subjets or folder and file name structure to search for file and analyze date,
  #' @wd: Working directory where Folders contain each folder. Use Finder or File explorer to set it (use slash to separate folders), and there HAS TO BE a slash(/) at the END
  #' @Task: Name of conditions or task evaluated or file name structure to search for file,
  #' @timeR: Time (s) of Ramp task for Regression Analysis for Enslaving Matrix calculation, by default between 5.5 and 7.5 s,
  #' @timeS: Begining and End of search time (in seconds) for ASAs. It should be specified by concatenate function (i.e c(1,20)). By default starts at steady state up to the end c(6,20),
  #' @freq: Sampling frequency of data acquisition. By default is set at 1000Hz,
  #' @cutfreq: Cut off frequency. By default is set at 20Hz,
  #' @order: Filter order for data processing.
  #' @J: Jacobian Matrix, by default is [1 1 1 1]
  
UCM<-function(subjects,wd=NULL,Task=NULL,timeR=c(5.5,7.5),timeS=c(6,20),freq=1000,cutfreq=0.02,order=2,J=matrix(c(1,1,1,1),ncol=1)){
  ## NEEDS to be checked: timeS 
if (missing(subjects)) 
    stop("There are no subjetcs specified to analize")
  if (missing(wd)) 
    warning("Working Directory was set by default as setwd()")
  if (missing(Task)) 
    stop("There are no Taks specified to analize")
  if (missing(timeS)) 
    warning("There are no specified time intervals to search for ASAs")
  if (missing(freq)) 
    warning("Sampling frequency by default is set at 1kHz")
  if (missing(cutfreq)) 
    warning("Cut off frequency for filtering by default is set at 20Hz")
  if (missing(order)) 
    warning("Butterworth filter is set as 2nd order by default")
  
  if (is.null(wd)) {wd<-setwd()}else{wd<-wd}# Check for Working Directory
  
  #### Function to load starts here!. First stablish WD, then list files and finally load.############
  Ramp<-lapply(subjects,function(x){
    setwd(paste(wd,as.character(x),sep=""));getwd() # Set working directory
    Files<-list.files(pattern=paste("Ramp",sep=""))
    
    # Confirm files and directory if necessary with print(getwd()) print(ramp)
    # Loading of Ramp files corresponds to fingers 1,3,4,2
    lapply(Files[c(1,3,4,2)],function(i){read.csv(i, header=FALSE, sep="\t",col.names=c("I","M","R","L")  )})})
  ###################### End Here ##############
  # Naming lists inside lists
  names(Ramp)<-subjects;# Lists names
  Ramp<-lapply(Ramp,function(x){names(x)<-c("I","M","R","L")
  x})# Nested lists names. Remember to retreive element to save names correctly
  
  #### Estimates Total force and derives Individual Forces
  Ramp_tot<-sapply(names(Ramp),simplify=FALSE,USE.NAMES=T,
                   function(x){sapply(names(Ramp[[x]]),simplify=FALSE,USE.NAMES=T,function(i)
                   {data.frame(rbind(apply(Ramp[[x]][[i]][-c(1:10),],diff,MARGIN=2),1),Ftot=rowSums(Ramp[[x]][[i]][-c(1:10),]))}
                   )}) # F_tot computation. I also remove first datum
  
  #### Data filtering ####
  library("signal")# install.packages("signal") if needed, load otherwise
  
  if(!is.null(cutfreq)){W<-cutfreq/freq}else{W=cutfreq}
  bf2<-butter(2, W=W, type = "low") # Filter characteristics
  
  Ramp_tot_f<- sapply(names(Ramp_tot),simplify=FALSE,USE.NAMES=T,
                      function(x){sapply(names(Ramp_tot[[x]]),simplify=FALSE,USE.NAMES=T,function(i)
                      {apply(Ramp_tot[[x]][[i]],filter,filt=bf2,MARGIN=2,names=i)}
                      )}) 
  ########################
  
  ###### Extract ascending or selected part of Ramp
  (if (is.null(timeR)) {
    init.r<-0 # After ramp begins
    final.r<-as.numeric(min(unlist(sapply(names(Ramp_tot_f),USE.NAMES=T,simplify=FALSE,
                                          function(x){sapply(names(Ramp_tot_f[[x]]),USE.NAMES=TRUE,simplify=FALSE,
                                                             function(i){dim(Ramp_tot_f[[x]][[i]])[1]})}))))} # before ramp ends
  else{
    init.r<-timeR[1]*freq;
    final.r<-timeR[2]*freq})
  Ramp_tot_f<-sapply(names(Ramp_tot_f),USE.NAMES=T,simplify=FALSE,
                     function(x){sapply(names(Ramp_tot_f[[x]]),USE.NAMES=TRUE,simplify=FALSE,
                                        function(i){Ramp_tot_f[[x]][[i]][init.r:final.r,]})})
  
  ########################
  # Check if regression is to be calculated based on derivatives of forces
  E<-sapply(names(Ramp_tot_f),USE.NAMES=T,simplify=FALSE,
            function(x){sapply(names(Ramp_tot_f[[x]]),USE.NAMES=TRUE,
                               function(i){lm(Ftot~I+M+R+L-1,data=as.data.frame(Ramp_tot_f[[x]][[i]]))$coefficients})})
  
  #### Plot every subject Individual Force
  par(mfrow=c(2,2))# Plotting data
  plot1<-sapply(names(Ramp),USE.NAMES=T,simplify=FALSE,
                function(x){sapply(names(Ramp[[x]]),USE.NAMES=TRUE,
                                   function(i){apply(Ramp[[x]][[i]][,-5],FUN=plot,type="l",
                                                     MARGIN=2,main="Individual finger forces during Ramp",ylab=paste("Master Finger",i),xlab=paste("Subject",x))})});
  #############Ramp Analisys ends here ######

#### Loading Trials
  asa<-sapply(subjects,USE.NAMES=TRUE,simplify=FALSE,function(x){
    setwd(paste(wd,as.character(x),sep=""));getwd() # Set working directory
    Files<-list.files(pattern=paste(x,Task,sep="")) # Check if necessary to append subjects name paste(Task,sep="")
    # Confirm files and directory if necessary with print(getwd()) print(ramp)
    sapply(Files,USE.NAMES=T,simplify=FALSE,function(i){
    read.csv(i, header=FALSE, sep="\t",col.names=c("I","M","R","L"))
    })})
################################
#### Derivatives and F tot calculation
  asa_df<-sapply(names(asa),simplify=FALSE,USE.NAMES=T,
    function(x){sapply(names(asa[[x]]),simplify=FALSE,USE.NAMES=T,function(i)
       {data.frame(
          rbind(apply(asa[[x]][[i]],function(j)
            {diff(filter(bf2,j))},MARGIN=2),0),Ftot=(filter(bf2,rowSums(asa[[x]][[i]]))))}
                )})   # Check for filtering
#################################
#### Time of ASA initation
  level=5 # Treshold to find tASA
  tsh<-function(x){which(diff(x)*100/max(diff(x))>=level)[1]} # Function to derive, normalize and find treshold

    init.s<-6000;
    final.s<-14000

  tASA<-sapply(names(asa_df),USE.NAMES=T,simplify=FALSE,
               function(x){sapply(names(asa_df[[x]]),USE.NAMES=TRUE,simplify=FALSE,
                  function(i){(tsh(asa_df[[x]][[i]][init.s:final.s,5]))+init.s})})
#################################
#### Plot aligned Forces
colors<-c("antiquewhite3","antiquewhite4","aquamarine4","chartreuse4","chocolate","azure4",
               "coral4","cornflowerblue","cyan3","darkblue","brown","brown1","darkgoldenrod","darkgoldenrod1",
               "deepskyblue4","firebrick","darkorchid4","forestgreen","darkseagreen","darkslategray","darkslategrey",
               "lightcoral","hotpink3","khaki")
colors<-sapply(names(asa_df),USE.NAMES=TRUE,simplify=FALSE,
               function(x){sapply(names(asa_df[[x]]),USE.NAMES=TRUE,simplify=FALSE,
                                        function(i){sample(colors,1,replace=FALSE)})})
pch<-sapply(names(asa_df),USE.NAMES=TRUE,simplify=FALSE,
            function(x){sapply(names(asa_df[[x]]),USE.NAMES=TRUE,simplify=FALSE,
                               function(i){sample(seq(1:24),1,replace=FALSE)})})
  t=seq(from=-2000, to=1000,by=1) # To plot appropriate time

  par(mfrow=c(1,1))
  sapply(names(asa_df),USE.NAMES=T,simplify=TRUE,
         function(x){
           plot(t,asa_df[[x]][[1]][((tASA[[x]][[1]]-2000):(tASA[[x]][[1]]+1000)),5], type="l",
                main=paste("Force Profile for subject",x),ylab="Force (N)", xlab="Time (ms)",ylim=c(0,5))
           sapply(names(asa_df[[x]]),function(i){
           lines(t,asa_df[[x]][[i]][((tASA[[x]][[i]]-2000):(tASA[[x]][[i]]+1000)),5],
                 col=as.character(colors[[x]][[i]]))
             abline(v=0,lty=2,col="red")
             legend(-2000,5,legend=paste("Trial",substr(names(asa_df[[x]]),18,20)),lty=1,col=as.character(colors[[x]]),cex=0.7)
             })})
##############################
##### From forces to Modes
  library("MASS") # install.packages("MASS) if necessary
  # Null space for Finger Force J=matrix(1,1,1,1,byrow=F)
  if (!is.null(J)) {J=J}else{J=matrix(c(1,1,1,1),ncol=1)}
  
asa_modes<-sapply(names(asa_df),USE.NAMES=TRUE,simplify=FALSE,
  function(x){sapply(names(asa_df[[x]]),USE.NAMES=TRUE,simplify=FALSE,
  function(i){as.data.frame(apply(asa_df[[x]][[i]][((tASA[[x]][[i]]-2000):(tASA[[x]][[i]]+500)),1:4],
    function(j){solve(E[[x]])%*%as.matrix(j)},
                  MARGIN=1))})}) # Inverse Matrix with row-wise multiplication
  
V<-sapply(asa_modes,USE.NAMES=TRUE,simplify=FALSE,function(x){
  (apply(array(as.numeric(unlist(x)), dim=c( 4,2501,length(x) )),function(i){list(i)},MARGIN=2))}) # Organized as [[t]][F,T]. Joined two functions because I'm too lazy to shorten the code

e=Null(J);e # Gives similar results to Matlab svd(t(J),nu = nrow(t(J)))$u[,1]. t(J)%*%e gives an all-zero Matrix
P=e%*%t(e) # Projection matrix, i.e P=P^2, P%*%J=t(P)%*%J

# Check code to verify results with Matlab
UCM_null<-svd(J, nu = ncol(J))$u[,1]

v_null<-sapply(V,USE.NAMES=TRUE,simplify=F,function(x)
  {sapply(x,USE.NAMES=TRUE,simplify=F,function(i)
    {sapply(i,USE.NAMES=TRUE,simplify=F,function(j)
    {apply(j,function(k){(P%*%k)},MARGIN=2)})})})# Seems to be working... confirm!

v_null<-sapply(v_null,USE.NAMES=TRUE,simplify=F,function(x){names(x)<-paste("Time",seq(1:length(x)),sep="")
x})# Nested lists names. Remember to retreive element to save names correctly


V<-sapply(V,USE.NAMES=TRUE,simplify=F,function(x){names(x)<-paste("Time",seq(1:length(x)),sep="")
x})# Nested lists names. Remember to retreive element to save names correctly


v_ort<-sapply(names(V),USE.NAMES=T,simplify=FALSE,function(x)
        sapply(x,USE.NAMES=T,simplify=FALSE,function(i){mapply(function(z,y) matrix(unlist(z)-unlist(y),ncol=dim(V[[x]][[1]][[1]])[2]), 
                              z=V[[x]],y=v_null[[x]],SIMPLIFY =FALSE)}))

# unlist(apply(v[,1:2,],function(x){P%*%x},MARGIN=2)#,dim=c(4,2,length(ASA_tot_f)))
Vucm<-lapply(v_null,function(x){lapply(x,function(i){lapply(i,function(j){sum(apply(j,1,var))})})})
Vort<-lapply(v_ort,function(x){lapply(x,function(i){lapply(i,function(j){sum(apply(j,1,var))})})})

library("Matrix") # Install if needed: install.packages("Matrix")
degUCM<-rankMatrix(v_null[[1]][[1]][[1]])[1] # Rango de la Matrix
degORT<-rankMatrix(v_ort[[1]][[1]][[1]])[1]

VucmNorm = sapply(Vucm,function(x){unlist(x)/degUCM},simplify=FALSE)
VortNorm = sapply(Vort,function(x){unlist(x)/degORT},simplify=FALSE)

Vtot = mapply(function(x,y) unlist(x)+unlist(y), x=Vucm,y=Vort,SIMPLIFY =FALSE)
a= mapply(function(x,y,a,b) (unlist(x)/unlist(a))-(unlist(y)/unlist(b)), x=Vucm,y=Vort,a=degUCM,b=degORT,SIMPLIFY =FALSE)
b=mapply(function(x,a,b) unlist(x)/(unlist(a)+unlist(b)),  x=Vtot, a=degUCM, b=degORT, SIMPLIFY =FALSE)

dV = mapply(function(x,y) unlist(a)/unlist(b), x=a, y=b, SIMPLIFY=FALSE)

  return(list("Ramp"=Ramp,"Ramp tot f"=Ramp_tot_f,"ASA"=asa,"ASA dF"=asa_df,"tASA"=tASA,
              "F-modes"=asa_modes,"V matrix"=V,"v null"=v_null,"v ort"=v_ort,"degUCM"=degUCM,"degORT"=degORT,
              "V UCM Norm"=VucmNorm,"V ORT Norm"=VortNorm,"Vtot"= Vtot, "a"=a,"b"=b,"dV"=dV,
              "E"=E))
}
ASA<-UCM(subjects=c("Paula","Angelo"),wd=c("/Users/angelobj/Box Sync/ASA old paper/Prueba/"),Task="ASA",cutfreq = 20)
rm(list = ls())
bf2<-butter(2, W=0.02, type = "low") 
cutfreq=20
freq=1000
timeR=c(5.5,7.5)
timeS=c(6,14)
wd=c("/Users/angelobj/Box Sync/ASA old paper/Prueba/")
Task="ASA"
subjects<-c("Paula","Angelo")
# Reads files names (ramp) and stores (Ramp) data
rm("cutfreq","freq","timeR","timeS","Task","E","init.r","init.s","final.r","final.s","asa","asa_tot",
   "asa_tot_f","asa_f","ASA","tASA","plot1", "Ramp","Ramp_tot","Ramp_tot_f","t","tASA","W","level")
