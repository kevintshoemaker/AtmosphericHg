


# Atmospheric Reactive Mercury project

#   Kevin Shoemaker
#   collaboration with Mae Gustin, Sarrah Dunham-Cheatham et al.


rm(list=ls())

# Load packages ----------

library(tidyverse)
library(readxl)
library(DataExplorer)
library(caret)
library(randomForest)
library(ranger)
library(faux)
library(mgcv)
library(hrbrthemes)

# Load functions ---------

sigfig <- function(vec, n=3){ 
  ### function to round values to N significant digits
  # input:   vec       vector of numeric
  #          n         integer is the required sigfig  
  # output:  outvec    vector of numeric rounded to N sigfig
  
  formatC(signif(vec,digits=n), digits=n,format="fg", flag="#") 
  
}      # end of function   sigfig

# Read in data ----------

mydf <- read_excel("GustinModelData2_KTS.xlsx",sheet=1,skip=2)


# Select variables for testing -----------

names(mydf)
tokeep <- c(
  "Duration_day",
  #"AirTUNR_C",    # exclude for now based on discussion
  "AHUNR_gm3",
  "CO_ppm",
  "NO2_ppb",     # run with and without this one
  "O3_ppm",
  "PM2.5_ugm3",
  "SO2_ppb",
  "MixDepthMax",
  "Eurasia",
  "EAsia",
  "POW",
  "POE",
  "Canada",
  "USA",
  "TRM_CEM",
  "TGOM_CEM",    # added per Mae
  "TPBM_PTFE",    # added per Mae
  "RMAreaO",
  "RMAreaBrCl",
  "RMAreaN",
  "RMAreaS",
  "RMAreaOrg"
)

mydf <- mydf %>% 
  select(all_of(tokeep))

# Define response and predictors ---------

rmas <- tokeep[grepl("RMA",tokeep)]
allresponses <- c("TRM_CEM","TGOM_CEM","TPBM_PTFE",rmas)  # total reactive mercury, gaseous oxidizded mercury, particle-bound mercury
duration_var <- "Duration_day"

predictors <- setdiff(tokeep,c(allresponses,duration_var))
length(predictors)  # 13 predictors

## divide response by duration
temp <- sapply(allresponses,function(t) mydf[[sprintf("%s_orig",t)]] <<- mydf[[t]]  )
temp <- sapply(allresponses,function(t) mydf[[t]] <<- mydf[[sprintf("%s_orig",t)]] / mydf[[duration_var]]  )

# Function: run analyses for single response variable -----------------------

response=allresponses[1]; no2=F
RunAnalysis <- function(response=allresponses[1], no2=F){
 
  # remove rows with NAs
  mydf_c <- mydf[complete.cases(mydf[,c(response,predictors)]),]
  # nrow(mydf_c)   # 34-35 observations
  # nrow(mydf) - nrow(mydf_c)   # 3-5 rows lost?
  
  # Correlation -------------
  
  mycor <- cor(mydf[,predictors],use="pairwise.complete.obs")
  # mycor
  
  torm <- caret::findCorrelation(mycor,cutoff=0.8)
  # predictors[torm]   # remove N02?
  
  torm_names <- predictors[torm]
  
  predictors2 <- predictors[setdiff(1:length(predictors),torm)]     # remove NO2
  
  if(no2) predictors2 <- unique(c(predictors2,"NO2_ppb"))  # add no2 back in
  
  length(predictors2)  # 12 predictors remain
  
  # PREV NOTE: this step removes air temp, which is the top variable if left in- but is highly correlated with many other variables.
  # NOTE: this step removes NO2, which is highly correlated with many other variables, especially CO.
  
  # Add random variables --------------
  
  mydf_c <- mydf_c %>% 
    mutate(
      rand1 = as.factor(sample(sample(letters[1:3], nrow(mydf_c), replace = TRUE))),
      rand2 = rnorm(nrow(mydf_c), mean = 10, sd = 2),
      corr1 = rnorm_pre(mydf_c[,response], mu = 10, sd = 2, r = 0.2, empirical = TRUE),
      corr2 = rnorm_pre(mydf_c[,response], mu = 10, sd = 2, r = 0.5, empirical = TRUE)
    )
  
  predictors3 <- c(predictors2,c("rand1","rand2","corr1","corr2"))
  
  # Transform data ---------
  
    ## scale all predictor variables
  mydf_c <- mydf_c %>%
    # Save categorical features as factors
    mutate(across(c("rand1"), as.factor)) %>%
    # Center and scale numeric features
    mutate(across(all_of(predictors3[sapply(mydf_c[,predictors3],is.numeric)]), ~as.vector(scale(.) ) )) 
  

  response_orig <- sprintf("%s_orig",response)
  
  
  # Visual exploration -----------
  
  mydf2 <- mydf_c[,c(response,response_orig,duration_var,predictors3)]
  
  # plot_intro(mydf2)
  # plot_bar(mydf2)
  # plot_correlation(mydf2)
  
  # histogram of response
  # resp <- mydf2 %>% 
  #   pull(all_of(response))
  # hist(resp)
  # hist(log(resp+0.1))  # log transformation may make sense here...
  
  mydf2[,response] <- log(mydf2[,response_orig]+0.1)  # do log transform
  
  # histograms of predictor variables
  # contpred <- setdiff(predictors3,"rand1")
  # temp <- lapply(contpred,function(x) { this <- ggplot(mydf2,aes(x=.data[[x]])) + geom_histogram(bins=10) + ggtitle(x); print(this)  })
  
  # RFE --------------
  
  control <- rfeControl(functions = rfFuncs, # random forest
                        method = "repeatedcv", # repeated cv
                        repeats = 10, # number of repeats
                        number = 10) # number of folds
  
  
  # Features
  x <- mydf2 %>%
    select(all_of(setdiff(predictors3,c(response,"rand1","rand2","corr1","corr2")))) %>%
    as.data.frame()
  
  # Target variable
  y <- mydf2[[response]]
  
  # run the RFE algorithm
  result_rfe1 <- rfe(x = x, 
                     y = y, 
                     sizes = c(1:length(setdiff(predictors3,c(response,"rand1","rand2","corr1","corr2")))),
                     rfeControl = control)
  
  # Print the results
  # result_rfe1
  
  # Print the selected features
  # predictors(result_rfe1)
  
  # result_rfe1
  
  # ggplot(data = result_rfe1, metric = "Rsquared") + theme_bw()
  no2t <- ifelse(no2,"NO2","noNO2")
  svg(sprintf("plots/rfe_%s_%s.svg",response,no2t),3,4)
    thisplot <- ggplot(data = result_rfe1, metric = "RMSE") + theme_bw() + 
      scale_x_continuous(limits=c(0,10),breaks = c(0, 2, 4, 6, 8, 10) )
    print(thisplot)
  dev.off()
  
  
  # select top variables
  set.size = 5 
  lm.vars <- result_rfe1$variables 
  lm.set <- lm.vars[lm.vars$Variables==set.size,  ] 
  lm.set <- aggregate(lm.set[, c("Overall")], list(lm.set$var), mean)
  lm.order <- order(lm.set[, c("x")], decreasing = TRUE)[1:set.size]
  tops <- lm.set[lm.order, ][["Group.1"]]
  
  ## run random forest with top variables -----------
  
  tops_plus <- c(tops,c("rand1","rand2","corr1","corr2"))
  
  formula <- as.formula(paste(response, "~", paste(tops_plus,collapse="+")))
  thismod <- ranger(formula,data=mydf2,importance = 'permutation')
  
  varimp <- importance(thismod)
  
  options(scipen=20)
  
  varimp <- varimp[order(varimp,decreasing=T)]
  
  svg(sprintf("plots/varimp_%s_%s.svg",response,no2t),4,4)
    par(mai=c(1,2,0.1,0.1))
    par(las=2)
    plot(length(varimp):1~varimp,yaxt="n",xlab="importance",ylab="",xlim=c(0,max(varimp)+0.1))
    axis(2,at=length(varimp):1,labels = names(varimp)  )
  dev.off()
  
  # effects plots -------------
  
  # bestvars <- predictors(result_rfe1)
  contvars <- setdiff(predictors3,"rand1")
  toplot2 <- names(varimp)[!names(varimp)%in%c("rand1","rand2","corr1","corr2")]
  toplot <- toplot2[1:min(4,length(toplot2))]
  
  svg(sprintf("plots/effectsRF_%s_%s.svg",response,no2t),6,5)
  par(mai=c(1,1,.1,.1))
  par(mfrow=c(2,2))
  p=1
  for(p in 1:length(toplot)){
    thisvar <- toplot[p]
    nd <- data.frame(x=seq(min(mydf2[[thisvar]]),max(mydf2[[thisvar]]),length=50))
    names(nd) <- thisvar
    
    othervars <- setdiff(contvars,thisvar)
    temp <- sapply(othervars,function(t) nd[[t]] <<- mean(mydf2[[t]])  )
    nd$rand1 <- factor("a",levels=levels(mydf2$rand1))
    #nd
    
    pred = predict(thismod,data=nd,type="response")$predictions
    
    plot(exp(pred)~nd[,1],type="l",xlab=thisvar,main="",ylab=response,xaxt="n")
    rug(mydf2[[thisvar]])
    axis(1,at=seq(min(mydf2[[thisvar]]),max(mydf2[[thisvar]]),length=5),
         labels=sigfig(pmax(0,seq(min(mydf2[[thisvar]]),max(mydf2[[thisvar]]),length=5) * sd(mydf[[thisvar]],na.rm=T) + mean(mydf[[thisvar]],na.rm=T)) ,2) )
  }
  dev.off()
  
  
  # interactions??  -----------------------
  
    # just test interactions among top 4 vars?
  # select top variables
  set.size = 4
  lm.vars <- result_rfe1$variables 
  lm.set <- lm.vars[lm.vars$Variables==set.size,  ] 
  lm.set <- aggregate(lm.set[, c("Overall")], list(lm.set$var), mean)
  lm.order <- order(lm.set[, c("x")], decreasing = TRUE)[1:set.size]
  tops <- lm.set[lm.order, ][["Group.1"]]
  
  # intdf <- expand.grid(tops,tops,stringsAsFactors = F)
  # intdf <- intdf[intdf$Var1!=intdf$Var2,]
  # intdf2 <- NULL
  # i=1
  # for(i in 1:nrow(intdf)){
  #   if(!paste0(intdf$Var2[i],intdf$Var1[i])%in%paste0(intdf2$Var1,intdf2$Var2)){
  #     intdf2 <- rbind(intdf2,intdf[i,])
  #   }
  # }
  # intdf2$IntStrength <- NA
  # i=1
  # for(i in 1:nrow(intdf2)){
  #   v1 <- intdf2$Var1[i]
  #   v2 <- intdf2$Var2[i]
  #   x1 <- as.numeric(quantile(mydf2[[v1]],seq(0,1,length=10)))
  #   x2 <- as.numeric(quantile(mydf2[[v2]],seq(0,1,length=10)))
  #   tempdf <- expand.grid(x1,x2)
  #   othervars <- setdiff(contvars,c(v1,v2))
  #   temp <- sapply(othervars,function(t) tempdf[[t]] <<- 0  )
  #   tempdf$rand1 <- factor("a",levels=levels(mydf2$rand1))
  #   #nd
  #   names(tempdf)[1:2] <- c(v1,v2)
  #   tempdf$pred = predict(thismod,data=tempdf,type="response")$predictions
  #   tempform <- as.formula(sprintf("pred ~ as.factor(%s)+as.factor(%s)",v1,v2))
  #   addmod <- lm(tempform,data=tempdf)
  #   tempdf$pred2 = predict(addmod)
  #   intdf2$IntStrength[i] <- sqrt(mean((tempdf$pred-tempdf$pred2)^2))
  # }
  # topint <- intdf2[which.max(intdf2$IntStrength),1:2]  
  # cor(mydf2[,toplot])  # double check if top variables are highly correlated
  
  # Conclusions so far --------------------
  
  #    if air temp (AirTUNR_C) is included, it tends to be the only important variable: resp is total Hg or oxygen species
  #    should test any final model against a model with only air temp- something tells me the air temp model will win out
  
  #    best model for RMAreaBrCl has:   
  
  par(mfrow=c(1,1))
  
  # FOR GENERAL MULTIPREDICTOR MODEL -----------------
  set.size = 5
  lm.vars <- result_rfe1$variables 
  lm.set <- lm.vars[lm.vars$Variables==set.size,  ] 
  lm.set <- aggregate(lm.set[, c("Overall")], list(lm.set$var), mean)
  lm.order <- order(lm.set[, c("x")], decreasing = TRUE)[1:set.size]
  tops <- lm.set[lm.order, ][["Group.1"]]
  
  thesepred <- predictors(result_rfe1)
  if(length(thesepred)>4) thesepred <- tops[1:4]
  if(length(thesepred)<2) thesepred <- tops[1]
  
  # thesepred <- unique(c(thesepred,as.vector(unlist(topint))))
  
  form2<-list()    # univariate formulas for visualization
  svg(sprintf("plots/scatterplots_%s_%s.svg",response,no2t),6,5)
    par(mai=c(1,1,.1,.1))
    par(mfrow=c(2,2))
    for(v in 1:min(4,length(thesepred))){
      form2[[v]] <- as.formula(paste(response,"~",thesepred[v]))
      form3 <- as.formula(paste("exp(",response,") ~",thesepred[v]))
      plot(form2[[v]],data=mydf2,ylab=sprintf("log %s",response),xlab=sprintf("%s_std",thesepred[v]))
      #plot(form3,data=mydf2)
    }
  dev.off()
  
  # run GAMs or GLMs! ----------------  
  form4 <- paste(response_orig,"~ offset(log(Duration_day))") 
  form5 <- form4
  # nonint <- setdiff(thesepred,as.vector(unlist(topint)))
  # if(thesepred[1]%in%nonint){
  #   nonint2 <- setdiff(thesepred,c(thesepred[1],as.vector(unlist(topint))))
  #   form4 <- paste(form4,"+ s(",thesepred[1],") + s(",topint[1],",",topint[2],")")
  #   form5 <- paste(form5,"+",thesepred[1]," + ",topint[1],"*",topint[2])
  #   if(length(thesepred)>2){
  #     for(v in 1:min(1,length(nonint2))){
  #       form4 <- paste(form4,"+ s(",nonint2[v],")")
  #       # form5 <- paste(form5,"+ ",thesepred[v],"")
  #     }
  #   }
  # }else{
  #   form4 <- paste(form4,"+ s(",topint[1],",",topint[2],")")
  #   form5 <- paste(form5,"+",topint[1],"*",topint[2])
  #   if(length(thesepred)>2){
  #     for(v in 1:min(2,length(nonint))){
  #       form4 <- paste(form4,"+ s(",nonint[v],")")
  #       form5 <- paste(form5,"+ ",nonint[v],"")
  #     }
  #   }
  # }
  
  form5 <- paste(form5, "+",paste(thesepred,collapse=" + "))
  
  # form5 <- paste(response_orig,"~ offset(log(Duration_day)) + ",thesepred[1],"")
  
  # form4 <- formula(form4)
  form5 <- formula(form5)
  # form4 <- as.formula(paste(response, "~", sprintf("(max*%s)/(K+%s)+ext",names(varimp)[1],names(varimp)[1])))
  # form5 <- as.formula(paste(response, "~", sprintf("Asym/(1 + exp((- %s)/scal)) + adj",names(varimp)[1]))) 
  # model1 <- nls(form5,data=mydf2,start=list(Asym=5,scal=1,adj=2))
  # model1 <- nls(form4,data=mydf2,start=list(max=5,K=1,ext=1))
  # model1 <- gam(form4,data=mydf2,family=Gamma(link="log"),select = T)
  # summary(model1)
  
  mydf2[[response_orig]] <- mydf2[[response_orig]]+0.1
  model2 <- glm(form5,data=mydf2,family=Gamma(link="log"))
  summary(model2)
  model3 <- step(model2)
  summary(model3)
  # plot(model3)
  
  # plot univariate effects plots  ----------------
  
  thesepred2 <- thesepred[thesepred%in% names(coef(model3))]
  # thesepred2 
  if(length(thesepred2)>2){
    svg(sprintf("plots/effectsplots_%s_%s.svg",response,no2t),6,5)
    par(mai=c(1,1,.1,.1))
    par(mfrow=c(2,2))
  }else{
    svg(sprintf("plots/effectsplots_%s_%s.svg",response,no2t),5,3)
    par(mai=c(1,1,.1,.1))
    par(mfrow=c(1,2))
  }
  
  v=1
  for(v in 1:length(thesepred2)){
    nd <- data.frame(
      temp = seq(range(mydf2[[thesepred2[v]]])[1],range(mydf2[[thesepred2[v]]])[2],length=100)
    )
    names(nd) = thesepred2[v]
    
    if(length(thesepred2)>1){
      othervars <- setdiff(thesepred2,thesepred2[v])
      for(v2 in 1:length(othervars)){
        nd[[othervars[v2]]] <- 0
      }
    }
    nd[[duration_var]] <- median(mydf2[[duration_var]])
    
    pred <- predict(model3,newdata = nd,se.fit=T,type="response") 
    
    # par(mfrow=c(1,1))
    plot(mydf2[[thesepred2[v]]],(mydf2[[response_orig]]/mydf2[[duration_var]])*median(mydf2[[duration_var]]),
         ylab=response,xlab=thesepred2[v],xaxt="n",pch=20,col=gray(0.6)) 
    lines(nd[[thesepred2[v]]],pred$fit,lwd=2)
    lines(nd[[thesepred2[v]]],pred$fit-2*pred$se.fit,lty=2)
    lines(nd[[thesepred2[v]]],pred$fit+2*pred$se.fit,lty=2)
    axis(1,at=seq(min(nd[[1]]),max(nd[[1]]),length=5),
         labels = sigfig(mean(mydf[[thesepred2[v]]],na.rm=T) + sd(mydf[[thesepred2[v]]],na.rm=T) * seq(min(nd[[1]]),max(nd[[1]]),length=5)  ,2)  )
    
  }
  
  dev.off()
  
  
  # visualize interaction --------------
  
  # if(any(grepl(":",names(coef(model3)) ))){
  #   svg(sprintf("plots/interaction_%s_%s.svg",response,no2t),5,5)
  #   thisint <- as.vector(unlist(topint))
  #   v1 <- thisint[1]
  #   v2 <- thisint[2]
  #   x <- seq(min(mydf2[[v1]]),max(mydf2[[v1]]),length=15)
  #   y <- seq(min(mydf2[[v2]]),max(mydf2[[v2]]),length=15)
  #   data <- expand.grid(X=x, Y=y)
  #   names(data) <- c(v1,v2)
  #   othervars <- setdiff(thesepred2,c(v1,v2))
  #   temp <- sapply(othervars,function(t) data[[t]] <<- 0)
  #   data[[duration_var]] <- median(mydf2[[duration_var]])
  #   data[[response]] <- predict(model3,newdata = data,type="response")
  #   
  #   # Heatmap 
  #   thisplot <- ggplot(data, aes(.data[[v1]], .data[[v2]], fill= .data[[response]])) + 
  #     geom_tile() +
  #     scale_fill_gradient(low="white", high="blue") +
  #     scale_x_continuous( #limits=c(max(0,min(x)),max(x)),
  #                        breaks = c(min(x),seq(min(x),max(x),length=4),max(x)),
  #                        labels = sigfig(pmax(0,c(min(x),seq(min(x),max(x),length=4),max(x)) * sd(mydf[[v1]],na.rm=T) + mean(mydf[[v1]],na.rm=T)),2) ) +
  #     scale_y_continuous( #limits=c(max(0,min(y)),max(y)),
  #                        breaks = c(min(y),seq(min(y),max(y),length=4),max(y)),
  #                        labels = sigfig(pmax(0,c(min(y),seq(min(y),max(y),length=4),max(y)) * sd(mydf[[v2]],na.rm=T) + mean(mydf[[v2]],na.rm=T)),2) )
  #   print(thisplot)
  #   dev.off()
  # }
  
  # test model
  theser <- DHARMa::simulateResiduals(model3)
  svg(sprintf("plots/GOF_%s_%s.svg",response,no2t),7,4)
  plot(theser)
  dev.off()
  DHARMa::testResiduals(theser)   # passes tests
  
  return(model3)
  
}  # end function "RunAnalysis"



# Run automated analysis -------------

allmodels_withno2 <- list()
allmodels_nono2 <- list()

r=1
for(r in 1:length(allresponses)){
  thisresponse <- allresponses[r]
  allmodels_nono2[[thisresponse]] <-  RunAnalysis(thisresponse,no2=F)
  allmodels_withno2[[thisresponse]] <-  RunAnalysis(thisresponse,no2=T)
}




## Test to see if these models (fitted to reno) work for predicting other sites...


# Read in new data


mydf_GUMO <- read_excel("GustinModelData_othersites.xlsx",sheet=2,skip=2)
mydf_NRGT <- read_excel("GustinModelData_othersites.xlsx",sheet=3,skip=2)
names(mydf)
names(mydf_GUMO)
names(mydf_NRGT)
mydf_GUMO$NO2_ppb   # note: no No2 info for GUMO
mydf_GUMO$O3_ppm
mydf_GUMO$RHUNR_percent
mydf_GUMO$RH_percent

r=3
for(r in 1:length(allresponses)){
  thisresponse <- allresponses[r]
  thismod <- allmodels_nono2[[thisresponse]]
  if("lm" %in% class(thismod)){
    thesevars <- names(thismod$coefficients)[-1]
    thesevars <- thesevars[!grepl(":",thesevars)]
    nd_GUMO <- mydf_GUMO[,c(thesevars,thisresponse,"Duration_day")]
    nd_NRGT <- mydf_NRGT[,c(thesevars,thisresponse,"Duration_day")]
    if(!any(sapply(nd_GUMO, function(t)  any(is.na(t))))){
      mydf_GUMO[[paste0(thisresponse,"pred")]] <- predict(thismod,nd_GUMO)
      plot(mydf_GUMO[[paste0(thisresponse,"pred")]],mydf_GUMO[[thisresponse]],
           xlab="predicted (from Reno model)",ylab="observed (GUMO)",main=thisresponse)
    }
    if(!any(sapply(nd_NRGT, function(t)  any(is.na(t))))){
      mydf_NRGT[[paste0(thisresponse,"pred")]] <- predict(thismod,nd_NRGT)
      plot(mydf_NRGT[[paste0(thisresponse,"pred")]],mydf_NRGT[[thisresponse]],
           xlab="predicted (from Reno model)",ylab="observed (GUMO)",main=thisresponse)
    }
    
  }
  
}

#  Make table of parameter values and standard errors:

predictors
predictors2
predictors3


GLMtable <- data.frame(
  ResponseVar = rep(allresponses,2),
  NO2 = rep(c(T,F),each=length(allresponses)),
  Intercept = NA
)

temp <- lapply(predictors,function(t) GLMtable[[t]] <<- NA )

r=1
for(r in 1:length(allresponses)){
  thisresponse <- allresponses[r]
  n=1
  for(n in 1:2){
    no2 = as.logical(n-1)
    
    if(no2){
      thismodel <- allmodels_withno2[[thisresponse]]
    }else{
      thismodel <- allmodels_nono2[[thisresponse]]
    }
   
    sum <- summary(thismodel)
    thesecoefs <- coefficients(thismodel)
    thesevarnames <- names(thesecoefs)
    coefs <- list()
    
    thisrow <- which(GLMtable$ResponseVar==thisresponse&GLMtable$NO2==no2)
    
    thisint <- sum$coefficients["(Intercept)","Estimate"]
    thisintse <- sum$coefficients["(Intercept)","Std. Error"]
    GLMtable$Intercept[thisrow] <- sprintf("%s (%s)",sigfig(thisint,3),sigfig(thisintse,2))
    
    v=1
    for(v in 1:length(predictors)){
      thisvar = predictors[v]
      if(thisvar %in% thesevarnames){
        thiscoef <- sum$coefficients[thisvar,"Estimate"]
        thiscoefse <- sum$coefficients[thisvar,"Std. Error"]
        GLMtable[[thisvar]][thisrow] <- sprintf("%s (%s)",sigfig(thiscoef,3),sigfig(thiscoefse,2))
      }else{
        GLMtable[[thisvar]][thisrow] <- ""
      }
    }
  }
  
}

write.csv(GLMtable,"GLMtable.csv",row.names=F)




# FINAL MODELS --------------


# Pick up here later...

# FOR TOTAL MERCURY -----------------

glm2 <- glm(form2,data=mydf2,family=Gamma(link="log"))
glm1 <- glm(form3,data=mydf2,family=Gamma(link="identity"))
summary(glm1)

theser <- DHARMa::simulateResiduals(glm1)
plot(theser)
DHARMa::testResiduals(theser)   # passes tests

library(effects)

plot(Effect("AirTUNR_C",glm1),type="response",add=T)

plot(exp(TRM_CEM)~AirTUNR_C,data=mydf2)

nd <- data.frame(
  AirTUNR_C = seq(range(mydf2$AirTUNR_C)[1],range(mydf2$AirTUNR_C)[2],length=100)
)

pred <- predict(glm1,newdata = nd,se=T,type="response") 

lines(nd$AirTUNR_C,pred$fit,lwd=2)
lines(nd$AirTUNR_C,pred$fit-2*pred$se.fit,lty=2)
lines(nd$AirTUNR_C,pred$fit+2*pred$se.fit,lty=2)


# Next steps? Try other responses?



