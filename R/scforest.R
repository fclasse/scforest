library(dplyr)
library(tidyr)
library(stringr)
library(partykit)
library(lavaan)
library(caret)
library(foreach)
library(parallel)
library(doParallel)
library(doRNG)
library(varhandle)
library(dHSIC)
library(strucchange)

################################################################################
### Hilfsfunktionen

scdata <- function(data,model,input){
  fitsi <- try(lavaan::cfa(model = model, data = data, estimator = "ML",do.fit=FALSE))
  manifs <- fitsi@Model@dimNames[[1]][[1]] #manifest variables
  latvars <- fitsi@Model@dimNames[[1]][[2]] #latent variables
  data = as.data.frame(data)
  data = data[,c(input,manifs)] %>% replace_na(list("-99"))
  return(list(manifs,latvars,data))
}

scform <- function(scores,input){
  formy = "";  for (i in 1:length(scores)){formy=paste(formy, paste(scores[i],"+"))  }
  formy = substr(formy,1,nchar(formy)-2); form = paste(formy,"~")
  for (i in 1:length(input)){ form=paste(form,paste(input[i],"+")) }
  form = substr(form,1,nchar(form)-2)
  return(form)
}

scsampling <- function(data,dbsamp,bagging){
  if (dbsamp){
    folds <- createFolds(y = rownames(data), k= 2)
    folds1 <- folds$Fold1
    folds2 <- folds$Fold2
  } else if (!is.null(bagging) && bagging < 1) {
    folds1=folds2 <- sample(x=rownames(data),size=bagging*nrow(data))
  } else{folds1=folds2=rownames(data)} 
  return(list(folds1,folds2))
}

scscores <- function(data,model,std.lv){
  csresult <- tryCatch({
    fit_num = lavaan::cfa(model = model, data = data,  estimator = "ML",std.lv = std.lv)
    fit_num_scores <- lavScores(fit_num)
    colnames(fit_num_scores) = str_replace_all(colnames(fit_num_scores), "[^[:alnum:]]", "")
    scores = colnames(fit_num_scores)
    data = cbind(data,fit_num_scores)
    list(data,scores)},error=function(e){warning("scforest warning: Your model does not converge with ML estimator. Try 'direct=FALSE'.");stop(e)})
  return(csresult)
}

sctrees <- function(data,model,input,ntrees,cutoff_rmsea,cutoff_loading,direct,dbsamp,bagging,ordered,std.lv,ctree_control){
  #Preprocess & gather Info
  dt=scdata(data,model,input);manifs=dt[[1]];latvars=dt[[2]];data=dt[[3]]
  
  #Scores ausrechnen wenn "direct"
  if(dbsamp) bagging=NULL
  if(!dbsamp & is.null(bagging)){direct=T}
  if(direct){sc=scscores(data,model,std.lv);data=sc[[1]];scores=sc[[2]]}
  
  #All Iterations Multicore
  ncores <- detectCores()-1
  cl <- makeCluster(spec=ncores) 
  registerDoParallel(cl)
  trees <- foreach(j=1:ntrees, .packages=c("stringr","lavaan","caret","partykit","strucchange"),.export=c("scdata","scform","scsampling","scscores")) %dorng% { 
    
    #### Data Partitioning
    sp = scsampling(data,dbsamp,bagging)
    treedata = data[which(as.numeric(rownames(data)) %in%  sp[[1]]),]
    datafit = data[which(as.numeric(rownames(data)) %in%  sp[[2]]),]
    
    #### MOB mit der ersten Haelfte der Daten
    #### Scores ausrechnen wenn nicht "direct"
    try({
      if(!direct){sc=scscores(data=treedata,model,std.lv);treedata=sc[[1]];scores=sc[[2]]}
      tree <- ctree(as.formula(scform(scores,input)),treedata, control = ctree_control  )
      },silent = T) 
    ####
    
    if(exists("tree")){ #erster error handler
      
      #### Re-Fit & Modelfit table
      fit_ord <- list()
      ni <- nodeids(tree, terminal = TRUE)
      modelfit <- data.frame()
      rls <- partykit:::.list.rules.party(tree)
      types <- sapply(input, function(y) ifelse(class(data[,y]) == "factor","LMuo","maxLM")  )#welcher parameter-stablilitaetstest
      typenames <- names(types)
      
      if(!(sum(rls=="")>0) & length(ni)!=0){ #zweiter Error hanlder
        
        ###Nodes refitten und Ergebnisse in Tabelle schreiben
        for(i in 1:length(ni)){
          mf<-c()
          data_refit <- subset(datafit,eval(parse(text=rls[i])) ) #richtiger Terminal node
          fit_ord[[i]] <- tryCatch({lavaan::cfa(model = model, data = data_refit, start = start, ordered = ordered, estimator = "WLS", std.lv = std.lv, control=list(iter.max=100))},error=function(e){return("n.c.")},warning=function(e){return("n.c.")})
          names(fit_ord)[[i]] <- paste0("node",ni[i])
          if(is.character(fit_ord[[i]])){next}
          pars <- tryCatch({cbind(fit_ord[[i]]@ParTable$lhs,fit_ord[[i]]@ParTable$op,fit_ord[[i]]@ParTable$rhs,fit_ord[[i]]@ParTable$est,fit_ord[[i]]@ParTable$se)},error=function(e){return("n.c.")})
          if(any(pars!="n.c.")) {
            parvars <- as.numeric(pars[pars[,2]=="~~",4]) 
            parload <- pars[pars[,2]=="=~",c(4,5)]
            zv <- as.numeric(parload[,1])/as.numeric(parload[,2]) #z-values 
            zvp <- sapply(zv[zv!=Inf],function(y) pnorm(q=y, lower.tail=FALSE)) #pvalues
            } else {next} #nur Varianzen & factor loadings
          if(any(parvars<0 | parvars>30) | any(zvp>0.05 & zvp!=1) | any(abs(as.numeric(parload[,1]))<cutoff_loading) ){next} #negative Varianzen oder zu kleine factor loadings --> next iteration
          
          mf[1] <- j; mf[2] <- ni[i]; if(dbsamp){mf[3] <- nrow(subset(treedata,eval(parse(text=rls[i])) ));mf[4] <- fit_ord[[i]]@Data@nobs[[1]]} else { mf[3] <- nrow(data_refit); mf[4] <- NA}
          mf[5] <- fitMeasures(fit_ord[[i]],"rmsea");mf[6] <- fitMeasures(fit_ord[[i]],"rmsea.ci.lower");mf[7] <- fitMeasures(fit_ord[[i]],"rmsea.ci.upper");mf[8] <- fitMeasures(fit_ord[[i]],"pvalue"); mf[9] <- rls[which(names(rls)==ni[i])]
          
          ### Parameter stability 
          if(fitMeasures(fit_ord[[i]],"rmsea")<cutoff_rmsea){  #stability test nur für models mit gutem fit
            datastable <- subset(data,eval(parse(text=rls[i])) ) #richtiger Terminal node aber mit gesamten Daten!
            fit_num <- tryCatch({lavaan::cfa(model = model, data = datastable, estimator = "ML", std.lv = std.lv, control=list(iter.max=100))},error=function(e){return("n.c.")}) 
            if(!is.character(fit_num)) {
              fluc_tests <- tryCatch({mapply(function(x,y) { sctest(fit_num, order.by = datastable[,y],   functional = x)$p.value }, types, typenames  )},error=function(e){return("n.c.")}) #was wenn "solution has not been found"?
              if(is.character(fluc_tests)){mf[10] <- "n.c."} else {mf[10] <- tryCatch({ifelse(any(fluc_tests < (0.05/length(fluc_tests)) & fluc_tests != 0),"unstable","stable")},error=function(e){return("n.c.")})  }
            } else {mf[10] <- "n.c."} #"not converged"
          } else {mf[10] <- "n.e."} #"not execuded"
        
          modelfit <- rbind(modelfit,mf)
        }
        
        if(nrow(modelfit)==0){modelfit<-"none converged"} else {
          modelfit <- as.data.frame(modelfit)
          colnames(modelfit) <- c("tree","node","n_fit","n_refit","RMSEA","RMSEA C.I. lower","RMSEA C.I. upper","p-value_chisq","decision_rule","parameter stability")
          modelfit <- modelfit[order(modelfit[,5],modelfit[,6]),] 
          }
        list(tree,fit_ord,modelfit,sp[[1]],sp[[2]]) 
      } else {tree}
    } else {NA}
  }
  stopCluster(cl)
  return(trees)
}

###################################
### Hauptfunktion

scforest.train <- function(data,model,input,ntrees=100,split=2,minsize=300,cutoff_rmsea=.05,cutoff_loading=.2,direct=F,dbsamp=T,bagging=NULL,ordered=NULL,std.lv=FALSE,ctree_control=ctree_control(minbucket=minsize, mtry= split)){
  
  #### Trees berechnen
  trees=sctrees(data,model,input,ntrees,cutoff_rmsea,cutoff_loading,direct,dbsamp,bagging,ordered,std.lv,ctree_control)
  
  #### Liste benennen
  for(i in 1:ntrees){
    names(trees)[[i]] <- paste0("iteration",i)
    if(length(trees[[i]])==1){names(trees[[i]]) <- paste0("tree",i)}
    if(length(trees[[i]])>1){names(trees[[i]]) <- c(paste0("tree",i),paste0("fit",i),paste0("modelfit",i),paste0("DataGrow",i),paste0("DataRefit",i))}
  }
  
  #### Bestnodes erstellen 
  mofi <- data.frame()
  for(i in 1:length(trees)){
    if ( length(trees[[i]])>1 && !is.character(trees[[i]][[3]]) ) mofi <- rbind(mofi,trees[[i]][[3]])
  }
  mofi = na.exclude(mofi[colSums(!is.na(mofi)) > 0]) 
  mofi_stable = mofi[mofi$RMSEA<cutoff_rmsea & mofi$`parameter stability`=="stable",]
  
  if(!(ncol(mofi_stable)==0)){
    #### trees "bereinigen"
    pred_trees <- trees[ as.numeric( unique(mofi_stable$tree))  ]
    names(pred_trees) <- names( trees[ as.numeric( unique(mofi_stable$tree)) ] )
    mofi_stable = mofi_stable[order(mofi_stable$RMSEA),]
  }
  suc_trees <- trees[sapply(trees, function(y) length(y)==5)]
  names(suc_trees) <- names(trees[sapply(trees, function(y) length(y)==5)])
  
  ##### Information
  info=list()
  info[[1]] <- ntrees; names(info)[[1]] <- "ntrees"
  info[[2]] <- length(suc_trees); names(info)[[2]] <- "successful_iterations"
  if(!(ncol(mofi)==0)){info[[3]] <- length(pred_trees)} else {info[[3]] <- 0}; names(info)[[3]] <- "pred_iterations" 
  info[[4]] <- input; names(info)[[4]] <- "input"
  info[[5]] <- direct; names(info)[[5]] <- "direct"
  info[[6]] <- dbsamp; names(info)[[6]] <- "dbsamp"
  if(is.null(bagging)){info[[7]] <- "NULL"} else {info[[7]] <- bagging}; names(info)[[7]] <- "bagging"
  info[[8]] <- cutoff_rmsea; names(info)[[8]] <- "cutoff_rmsea"
  info[[9]] <- cutoff_rmsea; names(info)[[9]] <- "cutoff_loading"
  info[[10]] <- ordered; names(info)[[10]] <- "ordered"
  info[[11]] <- model; names(info)[[11]] <- "model"

  #### Results erstellen
  scfor <- list()
  scfor[[1]] <- info; names(scfor)[[1]] <- "Info"
  scfor[[2]] <- mofi; names(scfor)[[2]] <- "Nodes"
  if(!(ncol(mofi_stable)==0)){scfor[[3]] <- pred_trees; names(scfor)[[3]] <- "Pred_trees"}
  if(!(ncol(mofi_stable)==0)){scfor[[4]] <- mofi_stable; names(scfor)[[4]] <- "Bestnodes"}
  scfor[[5]] <- suc_trees; names(scfor)[[5]] <- "Suc_trees"
  
  return(scfor)
}





#######################################################################
### scforest predict!!
scpredconf <- function(i,j,preds,latvar,conf,bestnodes,tree,indmethod){ #irgendwelche Probleme mit dopar... alles exportiert?
  try({
    types <- sapply(conf, function(y) ifelse(class(data[,y]) == "factor","discrete","gaussian")  )
    typenames <- names(types)
    
    if(nrow(preds)>0){
      varlist <- lapply(names(types),function(y){as.double(preds[,y])})
      ind_tests <- mapply(function(x,y) { dhsic.test(Y=as.matrix(preds[,latvar]),X=as.matrix(varlist[[which(names(types)==y)]]),matrix.input=F,method=indmethod,pairwise = F,kernel=c("gaussian",x))$p.value }, types, typenames  )
    } else {ind_tests <- NA}
    unconf_info<- c()
    unconf_info[1] <- tree[j]
    unconf_info[2] <- bestnodes[bestnodes$tree==tree[j],"node"][i]
    if(nrow(preds)>0){unconf_info[3] <- if(any( ind_tests < (0.05/length(ind_tests)) )){"confounded"} else {"not confounded"} } else {unconf_info[3] <- "no predictions"}
    unconf_info[4] <- bestnodes[bestnodes$tree==tree[j],"decision_rule"][i]
  },silent=T )#unabhängigkeit von allen input variablen (multivariate distribution)
  return(list(ind_tests,unconf_info))
}

scpredtree <- function(data,bestnodes,bestfits,tree,idvar,manifs,latvar,conf,indmethod,exclude_unconf,stdscores){
  ncores <- detectCores()-1
  cl <- makeCluster(spec=ncores) 
  registerDoParallel(cl)
  predis <- foreach(j=1:length(tree), .packages=c("lavaan","dHSIC","dplyr"), .export=c("scpredconf")) %dopar% {
    ind_tree=list()
    totpreds=data.frame()
    for(i in 1:length(bestnodes[bestnodes$tree==tree[j],2]) ){
      data_node <- data %>% filter(eval(parse(text=  bestnodes[bestnodes$tree==tree[j],"decision_rule"][i]  )) ) %>% select(all_of(c(idvar,manifs,conf)))
      preds <- tryCatch({ 
        scores <- lavPredict(bestfits[[j]][[i]],newdata = data_node)
        if(stdscores){scoressd <- apply(scores, 2, sd);scores <- t(  apply(scores, 1, function(x) x/scoressd)  )} 
        cbind(data_node,scores)   
        }, error=function(cond){return(data.frame())})
      
      ###Independence-test für jeden Node
      ind_tree[[i]] <- scpredconf(i,j,preds,latvar,conf,bestnodes,tree,indmethod)
      
      ###Preds nur verwenden wenn parameter unconfounded
      if(exclude_unconf){
        if( !(ind_tree[[i]][[2]][3] %in%c("no predictions","confounded"))){ totpreds <- rbind(totpreds,preds) }   
      } else { totpreds <- rbind(totpreds,preds) }
    }
    list(totpreds,ind_tree)
  }
  stopCluster(cl)
  return(predis)
}

scpredcompile <- function(predis,data,idvar,latvar,tree){
  goodpreds <- as.data.frame(data[,idvar]);colnames(goodpreds)=idvar;v=1
  for (i in 1:length(predis) ){ #every tree
    if(nrow(predis[[i]][[1]])>0) {preds <- predis[[i]][[1]][,c(idvar,latvar)]} else {next};v=v+1
    diffid <- setdiff(data[,idvar],preds[,idvar])
    diffid <- as.data.frame(cbind(diffid,rep(NA, length(diffid)) ))
    colnames(diffid) <- colnames(preds)
    preds <- as.data.frame(rbind(preds,diffid))
    if(nrow(preds) == nrow(data)  ){ #error handler
      preds <- preds[order(match(preds[,idvar], data[,idvar] )), ]
      goodpreds <- cbind(goodpreds,preds[,latvar])
      colnames(goodpreds)[v]<-paste0("tree",tree[i])
    }
  }
  goodpreds <- as.data.frame(goodpreds)
  if(ncol(goodpreds)>2){goodpreds$mean <- rowMeans(goodpreds[,2:ncol(goodpreds)], na.rm = TRUE)}
  
  unconf_table <-  data.frame()  
  for(i in 1:length(tree)){
    for(j in 1:length(predis[[i]][[2]])){
      unconf_var <- names(predis[[i]][[2]][[j]][[1]])[which(as.numeric(predis[[i]][[2]][[j]][[1]]) < (0.05/length(predis[[i]][[2]][[j]][[1]]))  )]
      if( !(identical(unconf_var, character(0)) |  is.null(unconf_var))  ) {unconf_row <- c(predis[[i]][[2]][[j]][[2]], unconf_var )} else {unconf_row <- c(predis[[i]][[2]][[j]][[2]], NA )}
      unconf_table <- rbind(unconf_table, unconf_row) 
    }
  }
  colnames(unconf_table) <- c("tree","node","confounded params","decision_rule","confounder")

  gathered <- list()
  gathered[[1]] <- goodpreds
  names(gathered)[[1]] <- "Goodpreds"
  gathered[[2]] <- unconf_table
  names(gathered)[[2]] <- "Unconfoundedness_table"
  return(gathered)
}

###Hauptfunktion
scforest.predict <- function(trained,data,idvar,conf=trained$Info$input,latvar=latvars[1],indmethod="gamma",exclude_unconf=TRUE,stdscores=FALSE){
  ### Warnung
  if(nrow(trained$Bestnodes) == 0){stop("scforest error: no suitable nodes produced in training.")}
  
  #### Gather Information
  data = as.data.frame(data)
  model = trained$Info$model
  ordered = trained$Info$ordered
  manifs = trained[[3]][[1]][[2]][[  which(sapply(trained[[3]][[1]][[2]], function(y) class(y)) == "lavaan")[1]  ]]@Model@dimNames[[1]][[1]]
  latvars = trained[[3]][[1]][[2]][[ which(sapply(trained[[3]][[1]][[2]], function(y) class(y)) == "lavaan")[1]  ]]@Model@dimNames[[1]][[2]]
  bestnodes = trained$Bestnodes[order(as.numeric(trained$Bestnodes$tree),as.numeric(trained$Bestnodes$node)),]
  tree = unique(bestnodes$tree)
  
  #### Pick models used for prediction out of trained model
  bstfts = sapply(trained[[3]], function(y) y[[2]]);bestfits = list()
  for(i in 1:length(unique(bestnodes$tree))){  bestfits[[i]] = bstfts[[i]][  which( nodeids(trained[[3]][[i]][[1]], terminal = TRUE)  %in%  bestnodes[bestnodes$tree==tree[i],2])  ] }
  names(bestfits) = names(trained[[3]])
  
  ### Warnungen
  if(length(bestfits) != length(tree)  ){stop("scforest fatal error: length of best_fits in trained model does not correspond with information in bestnodes.")} #weiter zum nächsten tree wenn liste mit best fits und bestnodes nicht übereinstimmen!
  if( is.null(data) ){stop("scforest error: please define data set for prediction including same manifest variables and partitioning variables as data set for training.")}
  if( is.null(latvar) ){stop("scforest error: please define latent variable of interest for LV-scores.")}
  
  ### Prediction & Independence tests...
  predis <- scpredtree(data,bestnodes,bestfits,tree,idvar,manifs,latvar,conf,indmethod,exclude_unconf,stdscores)
  
  ### Compile results
  scfor <- scpredcompile(predis,data,idvar,latvar,tree)
  return(scfor)
}








