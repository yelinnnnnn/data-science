#multiple_linear_regression , using only the "data.frame" type of dataset
#x=independent, column name or number / if multiple -> binding "c()"  
#y=dependent, column name or number 
my_lm<-function(x,y,data) { 
  Coefficient_s <-c("intercept",x)        
  p<-length(x)  
  xx<-rep(1,times=nrow(data[,x]))   
  x1<-cbind(xx,data[,x]) ; y1<-data[,y]    
  x<-as.matrix(x1) ; y<-as.matrix(y1)      
  B<-c(solve(t(x)%*%x) %*% t(x) %*% y)  
  y_hat<-x%*%B 
  e<-y-y_hat
  sse<-t(e) %*% e 
  r.v<- sse/(nrow(data)-p-1) 
  r.v<-as.double(r.v)
  cov.mat<-r.v * solve(t(x)%*%x)   
  S.E<-sqrt(diag(cov.mat)) 
  t_value<-B/S.E  ;  t_value<-as.vector(t_value) ; S.E<-as.vector(S.E)
  df_<-nrow(data)-2
  p_value_t<-c() 
  p_value_t<-2 * pt(-abs(t_value), df_)
  Estimate<-B      
  lmlm<-data.frame(Coefficient_s,Estimate,S.E,t_value,p_value_t)
  
  #R^2
  I<-diag(nrow(data))
  J<-matrix(data=1,nrow=nrow(data),ncol=nrow(data))
  sst<- (t(y) %*% y) - ((1/nrow(data)) * t(y)%*%J%*%y)
  R_sq<- 1-(sse/sst)
  R_sq<-as.vector(R_sq)
  #adjusted R^2
  ad_R_sq<- 1-((sse/sst)*(nrow(data)-1)/(nrow(data)-p-1))
  ad_R_sq<-as.vector(ad_R_sq)
  
  R_squared<-c("R^2"=R_sq, "adjusted_R^2"=ad_R_sq)
  
  #residuals standard error
  rse<-sqrt(r.v)
  
  #ANOVA table
  Sources<-c("Regression","Error","Total")
  SS<-c(sst-sse,sse,sst)
  df__<-c(p,nrow(data)-p-1,nrow(data)-1)
  MS<-c((sst-sse)/p ,r.v,"-")
  F_value<-c(((sst-sse)/p)/r.v,"-","-")
  p_f<- 1 - pf(as.double(F_value[1]), p, nrow(data)-p-1)
  p_value_f<- c(p_f,"-","-")
  
  anova_table<-data.frame(Sources,SS,df__ ,MS,F_value,p_value_f)
  
  ##residuals analysis(?)
  par(mfrow=c(2,2))
  plot(as.vector(y_hat),as.vector(e), main = "Residuals vs Fitted",
       xlab = "Fitted values", ylab = "Residuals")
  abline(h=0,col="blue")
  H<-x%*%solve(t(x)%*%x)%*%t(x)
  st_residual<- e/as.vector(sqrt(r.v*diag(I-H)))
  st_residual<-as.vector(st_residual)
  
  ranked_residuals <- sort(st_residual)
  theoretical_quantiles <- qnorm(ppoints(length(st_residual))) # theoretical quantiles 계산
  # Q-Q plot
  plot(theoretical_quantiles, ranked_residuals, main = "Q-Q Plot of Residuals",
       xlab = "Theoretical Quantiles", ylab = "Standardized residuals")
  abline(0, 1, col = "blue") 
  
  plot(as.vector(y_hat),abs(st_residual), main = "Scale - Location",
       xlab = "Fitted values", ylab = "sqrt(Standardized residuals)")
  leverage<-as.vector(diag(H))
  plot(leverage,as.vector(st_residual), main = "Residuals vs Leverage",
       xlab = "Leverage", ylab = "Standardized residuals")
  
  My_list <- list("Coefficients"=lmlm,"R-squared"=R_squared,"RSE"=rse, "ANOVA"=anova_table)
  return(My_list) }



#example
x_x<-rnorm(100,47.3,19)
x_x1<-rnorm(100,22.6,13.5)
y_y<-rnorm(100,67.2,40.1)  
data_f<-data.frame(x_x,x_x1,y_y)


my_lm(c(1,2),3,data=data_f)

lmlmlm<-lm(y_y~x_x+x_x1,data=data_f)
summary(lmlmlm)
par(mfrow=c(2,2))
plot(lmlmlm)