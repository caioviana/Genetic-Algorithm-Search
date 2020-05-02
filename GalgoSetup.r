#load data
library(galgo)
library(randomForest)

data <- readRDS("./illumina.prebatch.confirmedsplits.RDS")
data <- rbind(data$train,data$val,data$test)

#data <- read.csv("../data.csv",row.names=1)
#ref <- read.csv("reference.csv")



cls.use <- factor(data$sample.type)
data.use <- t(as.matrix(data[,-which(colnames(data) == "sample.type")]))
cls.tab <- table(cls.use)
cls.tab <- (1-(cls.tab/sum(cls.tab)))
cls.weight <- as.numeric(cls.tab)
names(cls.weight) <- names(cls.tab)

#options
maxSolutions <- 250	
maxBigBangs <- 250
chromosomeSize <- 30
goalFitness <- 0.95
saveFrequency <- 5

tit <- paste(levels(cls.use),collapse="-")
save.File <- gsub("\\?",".",paste(tit,".RData",sep=""))

rF.class <- function (chr, parent, tr, te, result) 
{
    require(randomForest)
	avg <- parent$data$avg
	train <- parent$data$data[tr, as.numeric(chr),drop=F]
	test <- parent$data$data[te, as.numeric(chr),drop=F]
	train.cl <- parent$classes[tr]
	test.cl <- parent$classes[te]
	trfd <- randomForest(train,train.cl,test,test.cl,ntree=parent$data$ntree,keep.forest=F,type="prob",mtry=parent$data$mtry,classwt=parent$data$cls.weight)
	#res <<- list(chr=chr,tr=tr,te=te,res=result,train=train,test=test,train.cl=train.cl,test.cl=test.cl)
	
	if(result == 1){	
		if(avg){
			classes <- levels(parent$classes)
			vals <- sapply(sapply(classes,function(x) which(test.cl == x)),function(y) sum(trfd$test$predicted[y] == test.cl[y])/length(y))
			xx <- mean(vals)
		}else{
			xx <- sum(trfd$test$predicted == test.cl)/length(test.cl)
		}
		if (is.numeric(xx)) 
			return(xx)
		else 0.01
	}else{
		return(trfd$test$predicted)
	}
}

galgo.search <- configBB.VarSel(
	data = data.use,
	classes = cls.use,
	main = tit,
	classification.method = "user",
	classification.train.error = "kfolds",
	classification.userFitnessFunc="rF.class",
	chromosomeSize = chromosomeSize,
	maxBigBangs = maxBigBangs,
	maxSolutions = maxSolutions,
	goalFitness = goalFitness,
	saveVariable = "galgo.search",
	saveFrequency = saveFrequency,
	saveFile = "galgo.search.parallel.Rdata",
	callBackFuncBB=NULL
	)
galgo.search$data$avg <- T
galgo.search$data$ntree <- 300
galgo.search$data$mtry <- 5
galgo.search$data$cls.weight <- cls.weight


saveObject(galgo.search, saveFile="galgo.search.parallel.Rdata")



