multiSAR <-
function(modelList,data,nBoot=999,crit="Info",norTest="lillie",alpha=0.05,verb=FALSE) {

######## multiSAR : model selection and averaging function 
#
# modelList a vector of model names : ex : c("power","expo","negexpo","monod","logist","ratio","lomolino","weibull")
# data : an mmSAR data object : a list (run data(data.galap)) for an example
# nBoot : the number of bootstrap resamples
# crit : "Info" or "Bayes"
# norTest : "lillie" or "shapiro" or "null"
# alpha : level of tests for normality and homoscedasticity
# verb : FALSE or TRUE
#
##########################################################

#will check if the data is OK
nPoints <- length(data$data[[1]])
nlig <- length(modelList)

############Criteria must be "Info" for AIC/AICc or "Bayes" for BIC
#choosing an IC criterion (AIC or AICc or BIC)
if(crit == "Info") {
	if ( (nPoints / 3) < 40 ) { IC = "AICc" } else { IC = "AIC"}
	} else {
	if(crit == "Bayes") { IC = "BIC" } else { stop("Criteria must be 'Info' for AIC/AICc or 'Bayes' for BIC")}
	}

library(numDeriv)


#GO
if(verb) cat("##### multiSAR #####\n#")
if(verb) cat("# Choosen criterion is ",IC,"\n")

#This should eventually be an option
#Test on data points (if one richness == 0 then the data point is deleted)
#isNull = which(data$data[[2]]==0)

#if (length(isNull)!=0) {
#	if(verb) cat("Dataset contained ",length(isNull)," zero abundance point(s) that was(were) deleted for analysis\n")
#	data$data = data$data[-isNull,]
#}#end of if isNull

#matrix of optimisation results
vars <- c("p1","p2","p3","AICc","D.AICc","AICcW","AIC","D.AIC","AICW","BIC","D.BIC","BICW","RSS","R2","Norm Stat","Norm p.val","Pearson","Pea p.val")
optimResult = matrix(0,nlig,length(vars))
colnames(optimResult) = vars
rownames(optimResult) <- modelList

#List of Jacobian and Hat matrix
matList = list()

#matrix of calculated values and residuals and transformed residuals
pointsNames <- paste("S",c(1:nPoints))
calculated <- residuals <- transResiduals <- matrix(0,nlig,length(pointsNames))
colnames(calculated) <- colnames(residuals) <- colnames(transResiduals) <- pointsNames
rownames(calculated) <- rownames(residuals) <- rownames(transResiduals) <- modelList

#vector of final values (model averaging)
finalVect <- vector("numeric",length(pointsNames))
names(finalVect) <- pointsNames

#getting optimization results for all models
for (i in 1:nlig){

	optimres = rssoptim(eval(parse(text=as.character(modelList[i]))),data,norTest,verb)

	for (j in 1:eval(parse(text=as.character(modelList[i])))$paramnumber) {optimResult[i,paste("p",j,sep="")] <- optimres$par[j]}
	optimResult[i,"AIC"] <- optimres$AIC
	optimResult[i,"AICc"] <- optimres$AICc
	optimResult[i,"BIC"] <- optimres$BIC
	optimResult[i,"RSS"] <- optimres$value
	optimResult[i,"R2"] <- optimres$R2
	optimResult[i,"Norm Stat"] <- optimres$normaStat
	optimResult[i,"Norm p.val"] <- optimres$normaPval
	optimResult[i,"Pearson"] <- optimres$pearson
	optimResult[i,"Pea p.val"] <- optimres$pearpval

	calculated[i,] <- optimres$calculated
	residuals[i,] <- optimres$calculated - data$data[,2]
	
	#jacobian and Hat Matrix
	#first data Point
	jacob = jacobian( eval(parse(text=as.character(modelList[i])))$rssfun,optimres$par,data=data$data[1,],opt=FALSE)

	for (k in 2:nPoints) {
	jacob = rbind(jacob,jacobian( eval(parse(text=as.character(modelList[i])))$rssfun,optimres$par,data=data$data[k,],opt=FALSE))
	}
	
	jacobbis <- t(jacob)%*%jacob
	s <- svd(jacobbis)
	jacobbismun = s$v%*%(diag(1/s$d))%*%(t(s$u))
	hatMat = jacob%*%jacobbismun%*%t(jacob)
	matList[[i]] <- list(jacob=jacob,hatMat=hatMat)

	#Residuals transformation from Davidson and Hinkley, 1997 "Bootstrap methods and their applications" p 259 eq (6.9)
	diagHatMat = diag(hatMat)
	transResiduals[i,] <- residuals[i,] - mean(residuals[i,])
	transResiduals[i,] <- transResiduals[i,] / sqrt( 1 - diagHatMat )
	
}#end of for

names(matList) = modelList


#Fit validation
flags <- vector("numeric",nlig)

if(norTest!="null"){
	for (i in 1:nlig) { if (optimResult[i,"Norm p.val"]<alpha || optimResult[i,"Pea p.val"]<alpha) {flags[i]<-"KO"} else {flags[i]<-"OK"}  }
}else{
	flags <- rep("OK",nlig)
}#eo ifelse

filtOptimResult <- subset(optimResult,flags=="OK")
filtCalculated <- subset(calculated,flags=="OK")
filtMatList <- matList[flags=="OK"]
filtModelList <- modelList[flags=="OK"]

if(verb) cat("# Valid models : ", paste(rownames(filtOptimResult),collapse=", "),"\n#\n")

#Models comparison

DeltaICvect <- vector()
akaikeweightvect <- vector()

filtNlig <- dim(filtOptimResult)[1]

for (i in 1:filtNlig){

	#Delta IC = ICi - ICmin 
	DeltaIC <- filtOptimResult[i,IC] - min(filtOptimResult[,IC])
	DeltaICvect <- c(DeltaICvect,DeltaIC)
}

for (i in 1:filtNlig){
	#Akaike Weigths
	akaikesum <- sum(exp( -0.5*(DeltaICvect)))
	akaikeweight <- exp(-0.5*DeltaICvect[i]) / akaikesum
	akaikeweightvect <- c(akaikeweightvect,akaikeweight)
}

	
columnDelta = paste("D.",IC,sep="")
filtOptimResult[,columnDelta] <- DeltaICvect
columnW = paste(IC,"W",sep="")
filtOptimResult[,columnW] <- akaikeweightvect

if(verb) cat("# Akaike weigths : ", paste(round(filtOptimResult[,columnW],4),collapse=", "),"\n#\n")


#Averaging
for (i in 1:nPoints) {
	finalVect[i] <- sum(akaikeweightvect*filtCalculated[,i])
}


#Averaging validation
avResiduals = data$data[[2]] - finalVect
shapRes= shapiro.test(avResiduals)
if(verb) cat("# Averaging residuals normality (p.value) : ",shapRes$p.value,"\n#\n")

cor <- cor.test(avResiduals,data$data[[1]])
if(verb) cat("# Averaging residuals/X values correlation (method: ",cor$method,") (Value,p.value) : ",cor$estimate,",",cor$p.value,"\n#\n")

################################################################################
#Bootstrapping residuals and model averaging
################################################################################

#Matrix of boot Samples
bootMatrix=matrix(0,nBoot,nPoints)

#array of optimisation results
optimBootResult = array(0,c(nlig,length(vars),nBoot),dimnames=list(modelList,vars,seq(1,nBoot)))

#array of calculated values
bootCalculated <- array(0,c(nlig,length(pointsNames),nBoot),dimnames=list(modelList,pointsNames,seq(1,nBoot)))

#flags for fitting validation
flags <- matrix(0,nlig,nBoot)

if(verb) cat("# Bootstrap Samples creation\n#\n")

#vector of choosen models
choosenModels = vector()

#test variable
nGoodBoot = 1 

while (nGoodBoot < nBoot+1) {
    
	test <- 1

	chousModel = filtModelList[rmultinom(1, 1, akaikeweightvect)==1]
	if(verb) cat("Boot Sample : ",nGoodBoot," Choose model : ",chousModel,"\n")

	choosenModels[nGoodBoot] = chousModel

   	while (test != 0 ) {

		for (l in 1:nPoints) {
			positives = transResiduals[chousModel,][transResiduals[chousModel,] > 0]
			negatives = transResiduals[chousModel,][transResiduals[chousModel,] < 0]

			if (calculated[chousModel,l] > 0 ) {
				vtci = negatives[abs(negatives) <= calculated[chousModel,l] ]
				vtci = c(vtci,positives)
				value = sample(vtci,1)
			} else {
				
				vtci = positives[positives >= abs(calculated[chousModel,l]) ]
				value = sample(vtci,1)
			}

			bootMatrix[nGoodBoot,l] <- calculated[chousModel,l] + value
		}#end of for


		#test if one species richness is negative
		test=length( which(bootMatrix[nGoodBoot,]<0) )

		if(verb) cat("BootSample : ",bootMatrix[nGoodBoot,],"\n")
		

    	}#end of while
	
	###########################################################
	#Do the model averaging for each bootstrap sample
	###########################################################

	if(verb) cat("# model optimization on bootstrap resamples\n#\n")

	for (k in 1:nlig){
		badBoot = FALSE

		optimres = tryCatch(rssoptim(eval(parse(text=as.character(modelList[k]))),data=list(name="bootSample",data=data.frame(a=data$data[[1]],s=bootMatrix[nGoodBoot,])),norTest,verb),error = function(e) { list(convergence=999) } )

		if (optimres$convergence != 0) {
			badBoot=TRUE
		} else { 
			if (sum(optimres$calculated)==0) {badBoot=TRUE} else { 
				for (j in 1:eval(parse(text=as.character(modelList[k])))$paramnumber){optimBootResult[k,paste("p",j,sep=""),nGoodBoot] <- optimres$par[j]}
				optimBootResult[k,"AIC",nGoodBoot] <- optimres$AIC
				optimBootResult[k,"AICc",nGoodBoot] <- optimres$AICc
				optimBootResult[k,"BIC",nGoodBoot] <- optimres$BIC
				optimBootResult[k,"RSS",nGoodBoot] <- optimres$value
				optimBootResult[k,"R2",nGoodBoot] <- optimres$R2
				optimBootResult[k,"Norm Stat",nGoodBoot] <- optimres$normaStat
				optimBootResult[k,"Norm p.val",nGoodBoot] <- optimres$normaPval
				optimBootResult[k,"Pearson",nGoodBoot] <- optimres$pearson
				optimBootResult[k,"Pea p.val",nGoodBoot] <- optimres$pearpval
				bootCalculated[k,,nGoodBoot] <- optimres$calculated

				#Fitting validation
				if(norTest!="null"){
					if (optimBootResult[k,"Norm p.val",nGoodBoot]<alpha || optimBootResult[k,"Pea p.val",nGoodBoot]<alpha || length(which(bootCalculated[k,,nGoodBoot]<0)) !=0 ) { flags[k,nGoodBoot]<-"KO"
					} else {
						flags[k,nGoodBoot]<-"OK"
					}#end of if/else on Shap and Corr
				}else{
						flags[k,nGoodBoot]<-"OK"
				}#eo if/else norTest!=NULL
			}#end of if/else on convergence 2
		}#end of if/else on convergence 1
			
	}#end of for k

	if ( length(which(flags[,nGoodBoot]!="KO")) == 0 ) badBoot=TRUE
	if ( length(which(flags[,nGoodBoot]==0)) != 0 ) badBoot=TRUE

	if (badBoot == FALSE) { 
			    #write the bootSample to a file
			    #bootFileName = paste("bootSamp_",data$name,".txt",sep="")
			    #bootText = paste("BootSample",nGoodBoot,"\n",sep="")
			    #cat(bootText,file = bootFileName,append=TRUE)
			    #write(bootMatrix[nGoodBoot,], file = bootFileName,ncolumns= nPoints, append = TRUE, sep = " ")
			    nGoodBoot <- nGoodBoot + 1
	}#end of if badBoot
}#end of while



#Applying the filter (flags)
#transform 3D table to list
filtOptimBootResult=vector("list", nBoot)
for (i in 1:nBoot) filtOptimBootResult[[i]] <- optimBootResult[,,i]
for (i in 1:nBoot) filtOptimBootResult[[i]] <- subset(filtOptimBootResult[[i]],flags[,i]=="OK")
filtBootCalculated = vector("list", nBoot)
for (i in 1:nBoot) filtBootCalculated[[i]] <- bootCalculated[,,i]
for (i in 1:nBoot) filtBootCalculated[[i]] <- subset(filtBootCalculated[[i]],flags[,i]=="OK")

bootHat = matrix(0,nBoot,nPoints)
nBadBoot = 0
f=0

	if(verb) cat("# model averaging on bootstrap resamples\n#\n")

for (k in 1:nBoot) {

	DeltaICvect <- vector()
	akaikeweightvect <- vector()
	filtNlig <- dim(filtOptimBootResult[[k]])[1]

	if (filtNlig != 0) {

	for (i in 1:filtNlig){

		#Delta IC = ICi - ICmin 
		DeltaIC <- filtOptimBootResult[[k]][i,IC] - min(filtOptimBootResult[[k]][,IC])
		DeltaICvect <- c(DeltaICvect,DeltaIC)
	}#end of for i

	for (i in 1:filtNlig){
		#Akaike Weigths
		akaikesum <- sum(exp( -0.5*(DeltaICvect)))
		akaikeweight <- exp(-0.5*DeltaICvect[i]) / akaikesum
		akaikeweightvect <- c(akaikeweightvect,akaikeweight)
	}#end of for i

	columnDelta = paste("D.",IC,sep="")
	columnW = paste(IC,"W",sep="")

	filtOptimBootResult[[k]][,columnDelta] <- DeltaICvect
	filtOptimBootResult[[k]][,columnW] <- akaikeweightvect

	#Averaging
	for (i in 1:nPoints) {
		bootHat[k,i] <- sum(akaikeweightvect*filtBootCalculated[[k]][,i])
	}#end of for i

	} else { bootHat[k,] <- rep(0,nPoints) }

} #end of for k

bootSort=apply(bootHat,2,sort)

res=list(data=data,models=modelList,optimRes=optimResult,filtOptimRes=filtOptimResult,calculated=calculated,filtCalculated=filtCalculated,averaged=finalVect,DeltaIC=DeltaICvect,akaikeweight=akaikeweightvect,avResiduals=avResiduals,shapAvRes=shapRes,corAvRes=cor,bootMatrix=bootMatrix,optimBootResult=optimBootResult,bootCalculated=bootCalculated,flags=flags,filtOptimBootResult=filtOptimBootResult,filtBootCalculated=filtBootCalculated,bootSort=bootSort,bootHat=bootHat,bootMatrix=bootMatrix,IC=IC) 

if(verb) cat("###################\n")

invisible(res)

} #end of multiSAR

