




spatial_durbin_nnls_wrap <- function(y, X,signs,Wrs=NULL, wts=NULL, unconstrained = FALSE, lag_response = TRUE, lag_predictors = TRUE){
  
  ## y should be centered, so no intercept needed
  if(!abs(mean(y)) < 1e-8){
    warning("response variable should be centered for nnls in order to fit without intercept")
  } 
  Xsigned <- X %*% Matrix::Diagonal(x=signs, names = colnames(X))

  if(!is.null(Wrs)){
    if(lag_predictors){
      WXsigned <- Wrs %*% Xsigned
      colnames(WXsigned) <- paste0("lag.", colnames(WXsigned))
      Xsigned <- cbind(Xsigned, WXsigned)
      signs <- c(signs, signs)
      rm(WXsigned)
    }
    
    if(lag_response){
      Wy <- Wrs %*% y
      colnames(Wy) <- "lag.y"
      Xsigned <- cbind(Wy, Xsigned)
      signs <- c(1, signs)
      rm(Wy)
    } 
  }  
 
  gc() ;
  #signs <- c(1, signs, signs) 
  #Xsigned <- cbind(Wy, Xsigned, WXsigned)
  #rm(WXsigned,Wy); gc()
  
  if(!is.null(wts)){
    a <- Matrix::crossprod(Matrix::Diagonal(x=sqrt(wts)) %*% Xsigned)
    b <- Matrix::crossprod(Matrix::Diagonal(x=wts) %*% Xsigned , y)
  } else {
    a <- Matrix::crossprod(Xsigned)
    b <- Matrix::crossprod(Xsigned, y)
  }
  if(unconstrained){
    beta <- solve(as.matrix(a), as.matrix(b))
  } else {
    beta <- RcppML::nnls(as.matrix(a), as.matrix(b))
  }
  #chk <- lm(y ~ as.matrix(Xsigned)-1, w = wts)
  # sqrt(diag(vcov(chk)))[1:10]
  # pred <- predict(chk, weights = wts, interval = "confidence")
  # pred2 <- predict(chk, weights = wts, interval = "prediction")
  # 
  xbeta <- Xsigned %*% beta ## get prediction
  rownames(beta) <- colnames(Xsigned)
  #res <- c(y - xbeta[,1])
  #if(!is.null(wts)){ 
  #  ### xbeta_var matches 'pred' SE above from wls (unconstrained)
  #  sigma2 <- (sum(wts * res^2)/(length(wts[wts > 0]) - ncol(Xsigned)))
  #  xbeta_var <- Xsigned %*% solve(a)
  #  xbeta_var <- Matrix::rowSums(xbeta_var * Xsigned)*sigma2
  #  sigma2 %*% solve(a)
  #  sig2_i <- sigma2 / wts
  #  acheeze <- 
  #    Matrix::crossprod(Matrix::Diagonal(x=sqrt(pmax(wts, 1e-16)*sig2_i)) %*% Xsigned)
  #  varbetainv  <- a %*% acheeze %*% a
  #  solve(varbetainv)     
  #} else {
  #  
  #}
  beta <- beta * signs ## convert back
  
  return(list(beta = beta
              ,xbeta = xbeta
              ,res = y - xbeta)) 
}