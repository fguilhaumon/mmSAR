transLink <-
function(x,boundType) {

	# fonction logit
	logit=function(y) {
		log(y/(1-y))
	}#end of logit

if (length(x) ==1) {
	res = switch(boundType,
		R = x,
		Rplus = log(x),
		unif = logit(x)
	)#end of switch
} else {

res=vector()

	for (i in 1:length(x)) {

		res[i] = switch(boundType[i],
				R = x[i],
				Rplus = log(x[i]),
				unif = logit(x[i])
			)#end of switch
	}#end of for
res
}#end of if/else
}#end of transLink

