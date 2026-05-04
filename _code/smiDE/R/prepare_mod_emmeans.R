

#' Function for creating a input list for emmeans cov.reduce argument
make_cov_reduce_list <- function(fixedterms,dat){
  cov_reduce_list <- list()
  for(termm in fixedterms){
    if(!is.factor(dat[[termm]])){
      cov_reduce_list[[termm]] <- function(x){ c(mean(x) + sd(x), mean(x))}   
    }
  }
  return(cov_reduce_list)
}


#' Function for converting nebula model to emmeans object
#' 
prepare_mod_emmeans <- function(mod, term, fitmethod, fixed_formula, dat, cov_reduce_list){
  ## nebula (neg binomial or poisson mixed models) 
  ## requires special output coercion to emmeans format 
  if(fitmethod == "nebula::nebula"){
    neb_summ <- nebula_summary(mod)
    bhat <- neb_summ$msumm$est
    names(bhat) <- neb_summ$msumm$term
    tt <- terms(fixed_formula)
    allvar <- as.character(attr(tt, "variables"))[-1]
    offsetv <- allvar[attr(tt,"offset")] 
    offsetuse <- 1
    if(length(offsetv) > 0) {
      offsetuse <- mean(eval(parse(text = offsetv), envir = as.data.frame(dat)))
    }

    mod <- 
      emmeans::qdrg(formula=fixed_formula
                    ,data=dat
                    ,vcov = neb_summ$neb_cov
                    ,coef=bhat
                    ,link='log'
                    ,offset=offsetuse
                    ,cov.reduce = cov_reduce_list
      )
    return(mod) 
  } else if (fitmethod=="spaMM::fitme"){
    #browser()
    ss <- capture_output(summary(mod))$result
    spamm_summ <- data.table::data.table(ss$beta_table, keep.rownames=TRUE)
    data.table::setnames(spamm_summ, "rn", "term")
    bhat <- spamm_summ[["Estimate"]]
    names(bhat) <- spamm_summ$term
    tt <- terms(fixed_formula)
    allvar <- as.character(attr(tt, "variables"))[-1]
    offsetv <- allvar[attr(tt,"offset")] 
    offsetuse <- 1
    if(length(offsetv) > 0){
      offsetuse <- mean(eval(parse(text=offsetv), envir=as.data.frame(dat)))
    }
    
    vcov_mat <- vcov(mod)
    attr(vcov_mat,"class") <- c("matrix", "array")
    
    mod <- 
      emmeans::qdrg(formula=fixed_formula
                    ,data=dat
                    ,vcov = vcov_mat
                    ,coef=bhat
                    ,link=family(mod)$link
                    ,offset=offsetuse
                    ,cov.reduce = cov_reduce_list
      )
    
  } else {
    return(mod) 
  }
}
