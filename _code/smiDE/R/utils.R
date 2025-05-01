
#' Function which prints a message using shell echo; useful for printing messages from inside mclapply when running in Rstudio
#' @param ... a list of messages pasted together.
message_parallel <- function(...){
  #system(sprintf('echo "\n%s\n"', paste0(..., collapse="")))
  system(sprintf('echo "%s"', paste0(..., collapse="")))
}


#' Function to suppress output from cat() in getting spaMM summary
#' Borrowed / copied from purrr package, (don't want all the purrr dependencies)
capture_output <- function (code) 
{
  warnings <- character()
  wHandler <- function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
  messages <- character()
  mHandler <- function(m) {
    messages <<- c(messages, conditionMessage(m))
    invokeRestart("muffleMessage")
  }
  temp <- file()
  sink(temp)
  on.exit({
    sink()
    close(temp)
  })
  result <- withCallingHandlers(code, warning = wHandler, 
                                message = mHandler)
  output <- paste0(readLines(temp, warn = FALSE), collapse = "\n")
  list(result = result, output = output, warnings = warnings, 
       messages = messages)
}



spamm_summary <- function(obj){
  #browser()
  ss <- capture_output(summary(obj))$result
  
  beta_summ <- data.table::data.table(ss$beta_table, keep.rownames=TRUE)
  beta_summ[,termtype:="fixed_effect"]
  data.table::setnames(beta_summ, "rn", "term")
  setnames(beta_summ
           ,c("Estimate", "Cond. SE", "t-value")
           ,c("est", "se", "t"))
  lambda_summ <- data.table::data.table(ss$lambda_table, keep.rownames=TRUE)
  setnames(lambda_summ, "Var.", "est")
  lambda_summ[,term:=paste0(Group, " ",Term, " Var.")]
  lambda_summ[,termtype:="random_effect"]
  lambda_summ[,`:=`(rn=NULL, Group=NULL, Term=NULL)]
  corrPars <- spaMM::get_ranPars(obj, which = "corrPars")
  corrPars <- corrPars[sort(names(corrPars))]
  cP <- unlist(corrPars)
  cP_summ <- data.table::data.table(term=names(cP), est = cP, termtype="random_effect_corrPars")
  
  ll_summ <- data.table(term=names(ss$likelihoods), est=ss$likelihoods)
  ll_summ[,termtype:="logLik"]
  msumm <- 
  rbindlist(list(
    beta_summ
    ,lambda_summ
    ,cP_summ
    ,ll_summ
  ),use.names=TRUE,fill=TRUE) 
  return(list(msumm=msumm)) 
}
nebula_summary <- function(obj){
  msumm <- data.table::melt(as.data.table(obj$summary), id.vars=c("gene_id"))
  msumm[,term:=gsub("^logFC_","", gsub("^p_", "", gsub("^se_", "", variable)))]
  msumm[,termidx:=.GRP,by=term]
  msumm[,idx:=1:.N]
  msumm[,attr:=gsub(paste0("_",term),"",variable,fixed = TRUE),by=.(idx)]
  disp_summ <- data.table(term = names(obj$overdispersion)
                          ,est = unlist(obj$overdispersion)
                          ,termtype="dispersion")
  neb_b <- dcast(msumm, gene_id + term + termidx ~ attr, value.var="value")[order(termidx)]
  neb_b[,`:=`(gene_id=NULL, termidx=NULL)]
  setnames(neb_b, c("logFC", "se", "p"), c("est", "se", "pval"))
  neb_b[,z:=est/se]
  neb_b[,termtype:="fixed_effect"]
  setcolorder(neb_b, c("term", "est", "se", "z", "pval"))
  neb_cov <- matrix(NA,nrow(neb_b), nrow(neb_b))
  neb_cov[lower.tri(neb_cov,diag=T)] = as.numeric(obj$covariance[1,])
  neb_cov[upper.tri(neb_cov)] = t(neb_cov)[upper.tri(neb_cov)]
  dimnms <- gsub("logFC_", "", grep("logFC", names(obj$summary), value=TRUE))
  dimnames(neb_cov) <- list(c(dimnms), c(dimnms))
  neb_b <- rbindlist(list(neb_b, disp_summ), use.names=TRUE, fill=TRUE)
  neb_b[,msg:=paste0(obj$convergence)] 
  
  return(list(msumm=neb_b,neb_cov=neb_cov)) 
}

.coef_summary <- function(mod, original_formula, fitmethod){
  summ_dtl <- list()
  #browser()
  #for(component in c("cond", "zi")){
  
    summ <- coef(summary(mod))
    if(!is.null(summ)){
      summ <- data.table(term=rownames(summ),summ)
      #c("est", "se", "z", "pval")
      newnames <- switch(fitmethod, 
                         "MASS::glmmPQL" = c("est", "se", "df", "t", "pval")
                         ,"MASS::glm.nb" = c("est", "se", "z", "pval")
                         ,"stats::glm" = c("est", "se", "t", "pval")
                         ,"lme4::lmer" = c("est", "se", "t")
                         )
      setnames(summ, old=names(summ)[2:ncol(summ)], new=newnames)
      
      re_terms <- lme4::findbars(original_formula) 
      has_re <- length(re_terms) > 0
       
      if(has_re & fitmethod %in% c("lme4::lmer", "MASS::glmmPQL")){
        rand_summ <- .random_summary(mod,original_formula, fitmethod) 
        summ <- rbindlist(list(summ[,termtype:="fixed"]
                               ,rand_summ
        ),use.names=TRUE,fill=TRUE)
      }
      if(fitmethod == "MASS::glm.nb"){
        summ <- rbindlist(list(summ, data.table(term="theta", est = mod$theta, se = mod$SE.theta)), use.names=TRUE,fill=TRUE)
      }
      if(fitmethod == "MASS::glmmPQL"){
         if(!is.null(mod$theta)) summ <- rbindlist(list(summ, data.table(term="theta", est = mod$theta)), use.names=TRUE,fill=TRUE)
      }
      summ_dtl <- data.table::copy(summ)
    }
  #}
  
  summ_out <- summ_dtl#[["cond"]]
  #mod_stats <- data.table(term=c("logLik", "dispersion")
  #                      ,est = c(as.numeric(logLik(mod)), summary(mod)$sigma)
  #                      ,termtype="model_stats")  
  #
  
  mod_stats <- data.table(term=c("logLik")
                        ,est = c(as.numeric(logLik(mod)))
                        ,termtype="model_stats")  
  
  #if(!is.null(summ_dtl[["zi"]])){
  #  summ_out <- rbindlist(list(summ_dtl[["cond"]][,component:="cond"]
  #                             ,summ_dtl[["zi"]][,component:="zi"]
  #                             ),use.names=TRUE,fill=TRUE)
  #  mod_stats[,component:="overall"]
  #}

  summ_out <- rbindlist(list(summ_out, mod_stats), use.names = TRUE,fill=TRUE)
  summ_out[is.na(termtype),termtype:="fixed"]
   
  return(summ_out)
}


.random_summary <- function(mod,original_formula, fitmethod){
  #browser()
  mtchdt <- 
    rbindlist(lapply(lme4::findbars(original_formula)
                     ,function(x) data.table(term=deparse(x)
                                             ,Groups=gsub("^.*\\ | ", ""
                                                          ,gsub("\\)",""
                                                                , deparse(x))
                                             )
                     )
    ))
 
  vc <- lme4::VarCorr(mod)
  reterm <- switch(fitmethod,
                  "MASS::glmmPQL" = lapply(attr(vc, "title"), function(x) gsub("\\ =.*$", "", x))[1]
                  ,"lme4::lmer" = attributes(vc)$names[1]
                  )
  if(fitmethod=="MASS::glmmPQL"){
    rand_dt <- data.table(Std.Dev.=as.matrix(vc)[,"StdDev"])
    rand_dt[,Name:=rownames(as.matrix(vc))]
    rand_dt[1:(.N-1),Groups:=reterm]
    rand_dt[is.na(Groups),Groups:="Residual"]
    
    setnames(rand_dt, old=c("Std.Dev.")
             , new=c("std.dev")
             ,skip_absent = TRUE
    )
    #if(any(grepl("corr", names(rand_dt)))){
    #  rand_dt[,lead1_corr:=shift(corr, 1, type="lead")]
    #  rand_dt[corr=="",corr:=lead1_corr]
    #  rand_dt[,lead1_corr:=NULL]
    #}
    #rand_dt[,lag1_grp:=shift(Groups, 1, type="lag")]
    #rand_dt[Groups=="",Groups:=lag1_grp]
    #rand_dt[,lag1_grp:=NULL]
    rand_dt[Groups=="Residual", Name:="Residual"]
  } else if (fitmethod=="lme4::lmer") {
    rand_dt <- lme4:::formatVC(lme4::VarCorr(mod))
    rand_dt <- data.table(rand_dt)
    setnames(rand_dt, old=c("Std.Dev.", "Corr")
             , new=c("std.dev", "corr")
             ,skip_absent = TRUE
    )
    
  
  }
  #rand_dt <- data.table(Std.Dev.=as.matrix(vc)[,"StdDev"])
  rand_dt[,Name:=rownames(as.matrix(vc))]
  rand_dt[1:(.N-1),Groups:=reterm]
  rand_dt[is.na(Groups),Groups:="Residual"]
  
  setnames(rand_dt, old=c("Std.Dev.")
           , new=c("std.dev")
           ,skip_absent = TRUE
  )
  if(any(grepl("corr", names(rand_dt)))){
    rand_dt[,lead1_corr:=shift(corr, 1, type="lead")]
    rand_dt[corr=="",corr:=lead1_corr]
    rand_dt[,lead1_corr:=NULL]
  }
  rand_dt[,lag1_grp:=shift(Groups, 1, type="lag")]
  rand_dt[Groups=="",Groups:=lag1_grp]
  rand_dt[,lag1_grp:=NULL]
  rand_dt[Groups=="Residual", Name:="Residual"]
  #keepnms <- intersect(c("Groups", "Name", "std.dev", "corr"), names(rand_dt))
  keepnms <- intersect(c("Groups", "Name", "std.dev"), names(rand_dt))
  mrand <- melt(rand_dt[,keepnms, with=FALSE]
                , id.vars=c("Groups", "Name")
                ,value.name="est")[est!=""]
  mrand[,est:=as.numeric(est)]
  
  mrand_dt <- merge(mrand,mtchdt, by="Groups",all.x=TRUE)
  mrand_dt[Groups=="Residual", term:="Residual"]
  mrand_dt[,termtype:=paste0("random,",Groups, ",", Name, ",", variable)]
  mrand_dt <- mrand_dt[,.(est,term,termtype)]
  return(mrand_dt) 
}

.replace_empty_rows_one_zero <- function(sm){
  dgt <- as(sm, "TsparseMatrix") 
  missing_rows <- setdiff(0:(nrow(dgt)-1L), dgt@i)
  dgt@i <- c(dgt@i, missing_rows)
  dgt@j <- c(dgt@j, rep(0L, length(missing_rows)))
  dgt@x <- c(dgt@x, rep(0L, length(missing_rows)))
  ret <- as(dgt, "CsparseMatrix")
  return(ret)
}

.replace_empty_cols_one_zero <- function(sm){
  dgt <- as(sm, "TsparseMatrix") 
  missing_cols <- setdiff(0:(ncol(dgt)-1L), dgt@j)
  dgt@j <- c(dgt@j, missing_cols)
  dgt@i <- c(dgt@i, rep(0L, length(missing_cols)))
  dgt@x <- c(dgt@x, rep(0L, length(missing_cols)))
  ret <- as(dgt, "CsparseMatrix")
  return(ret)
}

.replace_empty_sm <- function(sm){
  dims <- dim(sm)
  outsm <- new("dgCMatrix")
  outsm@Dim <- dims
  outsm@p <- integer(dims[2] + 1L)
  dimnames(outsm) <- dimnames(sm)
  return(outsm)
}


