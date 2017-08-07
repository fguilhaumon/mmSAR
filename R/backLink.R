backLink <-
function(x,boundType) {

	# fonction reciproque de la fonction logit
	invlogit=function(x) {
		1/(1+exp(-x))
	}#end of invlogit

if (length(x) ==1) {
	res = switch(boundType,
		R = x,
		Rplus = exp(x),
		unif = invlogit(x)
	)#end of switch
} else {

res=vector()

	for (i in 1:length(x)) {

		res[i] = switch(boundType[i],
				R = x[i],
				Rplus = exp(x[i]),
				unif = invlogit(x[i])
			 )#end of switch
	}#end of for
res
}#end of if/else
}#end of backLink

