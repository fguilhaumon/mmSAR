nlsModel2 <- function (form, data, start, wts, upper = NULL) 
{
   thisEnv <- environment()
   env <- new.env(hash = TRUE, parent = environment(form))
   for (i in names(data)) assign(i, data[[i]], envir = env)
   ind <- as.list(start)
   parLength <- 0
   for (i in names(ind)) {
      temp <- start[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp, envir = env)
      ind[[i]] <- parLength + seq_along(start[[i]])
      parLength <- parLength + length(start[[i]])
   }
   getPars.noVarying <- function() unlist(setNames(lapply(names(ind), 
                                                          get, envir = env), names(ind)))
   getPars <- getPars.noVarying
   internalPars <- getPars()
   if (!is.null(upper)) upper <- rep_len(upper, parLength)
   useParams <- rep(TRUE, parLength)
   lhs <- eval(form[[2L]], envir = env)
   rhs <- eval(form[[3L]], envir = env)
   .swts <- if (!missing(wts) && length(wts)) sqrt(wts) else rep_len(1, length(rhs))
   assign(".swts", .swts, envir = env)
   resid <- .swts * (lhs - rhs)
   dev <- sum(resid^2)
   if (is.null(attr(rhs, "gradient"))) {
      getRHS.noVarying <- function() {
         if (is.null(upper)) 
            numericDeriv(form[[3L]], names(ind), env)
         else numericDeriv(form[[3L]], names(ind), env, ifelse(internalPars < 
                                                                  upper, 1, -1))
      }
      getRHS <- getRHS.noVarying
      rhs <- getRHS()
   }
   else {
      getRHS.noVarying <- function() eval(form[[3L]], envir = env)
      getRHS <- getRHS.noVarying
   }
   dimGrad <- dim(attr(rhs, "gradient"))
   marg <- length(dimGrad)
   if (marg > 0L) {
      gradSetArgs <- vector("list", marg + 1L)
      for (i in 2L:marg) gradSetArgs[[i]] <- rep(TRUE, dimGrad[i - 
                                                                  1])
      useParams <- rep(TRUE, dimGrad[marg])
   }
   else {
      gradSetArgs <- vector("list", 2L)
      useParams <- rep(TRUE, length(attr(rhs, "gradient")))
   }
   npar <- length(useParams)
   gradSetArgs[[1L]] <- (~attr(ans, "gradient"))[[2L]]
   gradCall <- switch(length(gradSetArgs) - 1L, call("[", gradSetArgs[[1L]], 
                                                     gradSetArgs[[2L]], drop = FALSE), call("[", gradSetArgs[[1L]], 
                                                                                            gradSetArgs[[2L]], gradSetArgs[[2L]], drop = FALSE), 
                      call("[", gradSetArgs[[1L]], gradSetArgs[[2L]], gradSetArgs[[2L]], 
                           gradSetArgs[[3L]], drop = FALSE), call("[", gradSetArgs[[1L]], 
                                                                  gradSetArgs[[2L]], gradSetArgs[[2L]], gradSetArgs[[3L]], 
                                                                  gradSetArgs[[4L]], drop = FALSE))
   getRHS.varying <- function() {
      ans <- getRHS.noVarying()
      attr(ans, "gradient") <- eval(gradCall)
      ans
   }
   QR <- qr(.swts * attr(rhs, "gradient"))
   qrDim <- min(dim(QR$qr))
   if (QR$rank < qrDim) stop("singular gradient matrix at initial parameter estimates")
   
   getPars.varying <- function() unlist(setNames(lapply(names(ind), 
                                                        get, envir = env), names(ind)))[useParams]
   setPars.noVarying <- function(newPars) {
      assign("internalPars", newPars, envir = thisEnv)
      for (i in names(ind)) assign(i, unname(newPars[ind[[i]]]), 
                                   envir = env)
   }
   setPars.varying <- function(newPars) {
      internalPars[useParams] <- newPars
      for (i in names(ind)) assign(i, unname(internalPars[ind[[i]]]), 
                                   envir = env)
   }
   setPars <- setPars.noVarying
   on.exit(remove(i, data, parLength, start, temp, m))
   m <- list(resid = function() resid, fitted = function() rhs, 
             formula = function() form, deviance = function() dev, 
             lhs = function() lhs, gradient = function() .swts * attr(rhs, 
                                                                      "gradient"), conv = function() {
                                                                         if (npar == 0) return(0)
                                                                         rr <- qr.qty(QR, resid)
                                                                         sqrt(sum(rr[1L:npar]^2)/sum(rr[-(1L:npar)]^2))
                                                                      }, incr = function() qr.coef(QR, resid), setVarying = function(vary = rep(TRUE, 
                                                                                                                                                length(useParams))) {
                                                                         assign("useParams", if (is.character(vary)) {
                                                                            temp <- logical(length(useParams))
                                                                            temp[unlist(ind[vary])] <- TRUE
                                                                            temp
                                                                         } else if (is.logical(vary) && length(vary) != length(useParams)) stop("setVarying : 'vary' length must match length of parameters") else {
                                                                            vary
                                                                         }, envir = thisEnv)
                                                                         gradCall[[length(gradCall) - 1L]] <<- useParams
                                                                         if (all(useParams)) {
                                                                            assign("setPars", setPars.noVarying, envir = thisEnv)
                                                                            assign("getPars", getPars.noVarying, envir = thisEnv)
                                                                            assign("getRHS", getRHS.noVarying, envir = thisEnv)
                                                                            assign("npar", length(useParams), envir = thisEnv)
                                                                         } else {
                                                                            assign("setPars", setPars.varying, envir = thisEnv)
                                                                            assign("getPars", getPars.varying, envir = thisEnv)
                                                                            assign("getRHS", getRHS.varying, envir = thisEnv)
                                                                            assign("npar", length(seq_along(useParams)[useParams]), 
                                                                                   envir = thisEnv)
                                                                         }
                                                                      }, setPars = function(newPars) {
                                                                         setPars(newPars)
                                                                         assign("resid", .swts * (lhs - assign("rhs", getRHS(), 
                                                                                                               envir = thisEnv)), envir = thisEnv)
                                                                         assign("dev", sum(resid^2), envir = thisEnv)
                                                                         assign("QR", qr(.swts * attr(rhs, "gradient")), envir = thisEnv)
                                                                         return(QR$rank < min(dim(QR$qr)))
                                                                      }, getPars = function() getPars(), getAllPars = function() getPars(), 
             getEnv = function() env, trace = function() {
                cat(format(dev), ": ", format(getPars()))
                cat("\n")
             }, Rmat = function() qr.R(QR), predict = function(newdata = list(), 
                                                               qr = FALSE) eval(form[[3L]], as.list(newdata), env))
   class(m) <- "nlsModel"
   m
}




rssoptim <-
function(model,data,norTest="lillie",verb=TRUE){

####################################################
#                  INPUTS                          #   
####################################################
#                                                  #
# model: the model -list- (ex:  power)             #
# data: the data -data.frame-                      #
# norTest : shapiro OR lillie                      #
# verb : print info in console?                    #
####################################################

	if (norTest == "lillie") library(nortest)


	data.name=data$name
	data=data$data[,1:2]

	if (verb){
 	   cat("**********************************\n")
 	   cat("----------------------------------\n")
  	   cat("-FITTING: Residual Sum of Squares-\n")
 	   cat("< MODEL: ",model$name,">\n")
  	   cat("<  DATA: ",data.name,">\n")
  	   cat("----------------------------------\n")
	}

  	  l <- data[[2]]
  	  a <- data[[1]]

	if (verb){
  	  cat("--------------DATAS---------------\n")
   	 cat("A:",a,"\n")
   	 cat("S:",l,"\n")
   	 cat("----------------------------------\n")
	}

	#paramters bounds
	parLim = model$parLim

	#Transformed initial paramters
	start <- model$init(data)
	if(verb){cat("start :", start,"\n")}

	for (i in 1:length(start)) {
		if(parLim[i]!="R"){
			if(start[i]<=0){
				start[i]=0.1
			}
		} 
		if(parLim[i]=="unif"){
			if(start[i]>1){
				start[i]=0.9
			}
		}		

	}	

	startMod = transLink(start,parLim)
    
	#RSS function
	rssfun <- model$rssfun

	
	if (verb){
	    cat("------INITIAL VALUES--------------\n")
	    cat(start,"\n")
	    cat("----------------------------------\n")
	    cat("-transformed INITIAL VALUES-------\n")
	    cat(startMod,"\n")
	    cat("----------------------------------\n")

	}


	res1=optim(startMod,rssfun,hessian=FALSE,data=data,opt=TRUE, method="Nelder-Mead", control=list(maxit=50000)) # CG SANN Nelder-Mead BFGS

	#Backtransformation of parameters values

	res1$par = backLink(res1$par,parLim)

	#if (res1$par[2] <= 0) res1$par[2] <- start[2]
	names(res1$par)=model$paramnames

	l.calc=NULL
	l.calc = as.vector(model$fun(res1$par,data))
	residu = as.vector(l - l.calc)

	res2 = list(startvalues=start,data=data,model=model,calculated=l.calc,residuals=residu)

	#Residuals normality test

	if(norTest=="lillie"){

		if(length(l)<5) {

			cat("WARNING : the Lilliefors test cannot be used with less than 5 data points -> switching to the Shapiro test \n")

			norTest="shapiro"

		}#eo if length
	}#eo if lillie

	if(norTest=="shapiro"){

		if(length(l)<3) {

			cat("WARNING : the Shapiro test cannot be used with less than 3 data points -> residuals normality will not be checked\n")

			norTest="null"

		}#eo if length
	}#eo if shapiro

	normaTest = switch(norTest, "shapiro" = shapiro.test(residu) , "lillie" = lillie.test(residu) , "null"=list(statistic=NA,p.value=NA) )

	#Homogeneity of variance
	cor <- cor.test(residu^2,data[[1]])

	#Calcul des criteres R2a, AIC, AICc, BIC

	#variables communes
	n = length(a)
	P = model$paramnumber + 1  # + 1 pour la variance estimee

	#R2 (this definition of the R2 alloaws the comparison of all models)
	R2 <-  1 - ( (res1$value) /  sum((l - mean(l))^2) )

	#AIC
	AIC = n * log(res1$value / n) + 2 * P

	#AICc
	AICc = n * log(res1$value / n) + 2*P*(n / (n - P - 1))

	#BIC
	BIC = n *log(res1$value / n) + log(n) * P

	res3 = list(AIC=AIC, AICc=AICc, BIC=BIC, R2=R2)

	#estimates signifiance and confidence interval (95%)

	#constructing a nlsModel object
	nMod <- nlsModel2(model$form,data,res1$par)

	#number of parameters
	p <- model$paramnumber

	#residuals degrees of freedom
	rdf <- n - p

	#residuals variance
	resvar <- res1$value / rdf

	#calculating the inverse of the upper triangular factor
	#of the gradient array at estimated parameter values
	XtXinv <- chol2inv(nMod$Rmat())
	dimnames(XtXinv) <- list(model$paramnames, model$paramnames)

	#formating the table of estimates, standard eroor, t value and significance of parameters
	se <- sqrt(diag(XtXinv) * resvar)
    	tval <- res1$par/se
    	param <- cbind(res1$par, se, tval, 2 * pt(abs(tval), rdf, lower.tail = FALSE))
	dimnames(param) <- list(model$paramnames, c("Estimate", "Std. Error", 
        "t value", "Pr(>|t|)"))

	#95% confidence interval
	conf <- matrix(c(param[,"Estimate"] - 2 * param[,"Std. Error"], param[,"Estimate"] + 2 * param[,"Std. Error"]),p,2)
	colnames(conf) <- c("2.5%","97.5%")
	
	sigConf <- cbind(param,conf)

	if(verb){
		cat("----------FINAL VALUES-----------\n")
		print(sigConf)
		cat("\n")
		cat("----------------------------------\n")
		cat("RSS.value:",res1$value,"\n")
		cat("----------------------------------\n")
		cat("------RESIDUALS NORMALITY --------\n")
		if (norTest == "shapiro") {
			cat("Shapiro Test, W = ",normaTest$statistic,"\n")
			cat("Shapiro Test, p.value = ",normaTest$p.value,"\n")
		} else {
			cat("Lilliefors Test, D = ",normaTest$statistic,"\n")
			cat("Lilliefors Test, p.value = ",normaTest$p.value,"\n")
		}
	
        	cat("------HOMOGENEITY OF VARIANCE ---------\n")
		cat("Pearson coef. = " ,cor$estimate,"\n")
		cat("cor. Test p.value  = " ,cor$p.value,"\n")
	
		cat("--------- CRITERIA -------------------\n")
		cat("AIC : ",AIC,"\n")
		cat("AICc : ",AICc,"\n")
		cat("BIC : ",BIC,"\n")
		cat("R2 : ",R2,"\n")
		cat("**********************************\n")
	}#end of if verb

    res = c(res1,list(sigConf=sigConf,pearson=cor$estimate,pearpval=cor$p.value,normaTest=norTest,normaStat=normaTest$statistic,normaPval=normaTest$p.value),res2,res3,list(data.name=data.name))
    
    invisible(res)

#END OF RSSOPTIM
}

