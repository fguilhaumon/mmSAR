##########################################
# mmSAR plot functions
##########################################

#these functions are very basic and MUST be enhanced


#basic data plot
plot.mmSAR.data <- function(x,...) {

	plot(x$data,main=x$name,xlab="Area",ylab="Richness",...)
}


#basic model plot
plot.mmSAR.model <- function(model,pars=NULL,xs=1:100,...){

	if(is.null(pars)) error("Please supply parameter values")

	xs <- list(xs,c(1,2,3,4,5,6,7,8,9))

	xlims <- range(xs)[[1]]

	ys <- model$fun(pars,xs)

	main.text <- paste(model$name," (",paste(pars,collapse=" , "),")",sep="")

	plot(ys~xs[[1]],xlab="Area",ylab="Richness",main=main.text,...)


}
