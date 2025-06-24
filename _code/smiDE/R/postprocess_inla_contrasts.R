

postprocess_inla_contrasts <- function(lclists, mod, inlanames, realnames, modtime=NULL, error_msg){
  ## inla model summary 
  if(!is.null(error_msg)){
    emm_list <- lapply(lclists$terms, function(xx){
      data.table(msg = error_msg, term = xx)
    })
    pw_list <- lapply(lclists$terms, function(xx){
      data.table(msg = error_msg, term = xx)
    })
    onevrest_list <- lapply(lclists$terms, function(xx){
      data.table(msg = error_msg, term = xx)
    })
    onevall_list <- lapply(lclists$terms, function(xx){
      data.table(msg = error_msg, term = xx)
    })
    msumm <- data.table(msg = error_msg)
    
    return(
      list(model_summary = msumm
           ,pairwise_list = pw_list 
           ,onevrest_list = onevrest_list
           ,onevall_list =  onevall_list
           ,emm_list = emm_list
           ,terms = lclists[["terms"]] 
      )
    )
      
  } else {
    modsumm <- summary(mod)
    msumm <- rbindlist(list(
      data.table(modsumm$fixed, keep.rownames=TRUE)
      ,data.table(modsumm$hyperpar, keep.rownames=TRUE)
    )
    ,use.names=TRUE,fill=TRUE
    )
    setnames(msumm, "rn", "term")
    names(realnames) <- inlanames
    msumm[term %in% inlanames,term:=realnames[term]]
    
    if(!is.null(modtime)){
      msumm <- 
        tryCatch({
          suppressMessages(
            cbind(msumm, rbind(modtime))[,`:=`(user.child=NULL,sys.child=NULL)]
          )
        }
        ,error=function(e){return(msumm)}
        )
    }
    
    
    lincomb_ests <- data.table(mod[["summary.lincomb.derived"]], keep.rownames=TRUE)
    lincomb_ests[,ID:=NULL]
    setnames(lincomb_ests, "rn", "contrast")
    
    if(mod$family == "nbinom2"){
      ### log-link results to exponentiate
      logc <- c("mean", grep("quant$", names(lincomb_ests), value=TRUE), "mode")
      lincomb_ests[,(c(logc)):=lapply(.SD, exp), .SDcols=c(logc)]
    }
      
    outl <- list("emm_list" = list(), "pairwise_list" = list(), "onevrest_list" = list(), "onevall_list" = list())
    keepnames <- c("level", "term", "category", "wts")
    for(term in names(lclists[["emm_list"]])){
      outl[["emm_list"]][[term]] <- 
        merge(lincomb_ests
              ,lclists[["emm_list"]][[term]][,keepnames,with=FALSE]
              ,by.x="contrast"
              ,by.y="level"
        )
      setnames(outl[["emm_list"]][[term]], "contrast", "level")     
    }
    
    
    keepnames <- c("contrast", "counts_1", "propnz_1", "ncells_1", "counts_2", "propnz_2", "ncells_2", "term")
    for(contrasttype in c("pairwise_list", "onevrest_list", "onevall_list")){
      if(contrasttype %in% names(lclists)){
        for(term in names(lclists[[contrasttype]])){
          havnames <- intersect(keepnames, names(lclists[[contrasttype]][[term]]))
          tmp <- 
            merge(lincomb_ests
                  ,lclists[[contrasttype]][[term]][,havnames,with=FALSE]
                  ,by="contrast"
            )
          #if(linkused=="identity" & contrasttype=="pairwise_list"){
          #  fc_pw <- 
          #    merge(emm_list[["jj"]][,c("level", "response", "ky","category"),with=FALSE]
          #          ,emobdt_jj[,c("level", "response", "ky", "category"),with=FALSE]
          #          ,by="ky",allow.cartesian=TRUE, suffixes=c("_1", "_2"))
          #  fc_pw <- fc_pw[category_1 < category_2]
          #  fc_pw[,fold_change:=response_1 / response_2]
          #  fc_pw[,contrast:= paste0(level_1, " / ", level_2)] 
          #  pw[,level_1:=tstrsplit(contrast," - ")[[1]]]
          #  pw[,level_2:=tstrsplit(contrast," - ")[[2]]]
          #  pw <- merge(pw, fc_pw[,.(level_1, level_2, fold_change)]
          #              ,by=c("level_1", "level_2")
          #  )
          #  
          #  rm_names <- c("level_1", "level_2")
          #  pw[,(rm_names):=NULL]
          #} 
          outl[[contrasttype]][[term]] <- tmp
        }
      }
    }
    
    return(
      list(model_summary = msumm
           ,pairwise_list = outl[["pairwise_list"]] 
           ,onevrest_list = outl[["onevrest_list"]]
           ,onevall_list = outl[["onevall_list"]]
           ,emm_list = outl[["emm_list"]]
           ,terms = lclists[["terms"]] 
      )
    ) 
  }
  
}