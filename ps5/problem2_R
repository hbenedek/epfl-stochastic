library(Rlab)
library(readr)
## I used the last 5 years of S&P 500

sp500 <- read_table2("sp500.csv")

data = sp500[,'Close*'] #We look at the close price 
data = rev(data$`Close*`) #Sort from oldest to newest

log_return = log(data[2:1258] / data[1:1257])
R = abs(mean(log_return))

positives = subset(log_return, log_return > 0)
p = length(positives) / length(log_return) #p>0.5


Can_I_make_money = function(a, b, R, p){
  P = 0; T = 0;
  
  while ((- b/100 < P) & (P < a/100)) {
    if( rbern(1, p)){
      P = P + R
    } else {P = P-R}
    T=T+1
  }
  
  return(c(P,T))
}

## Try 100 times with a=b=1%
P_1 = data.frame(P = double(), Time= double())
for (i in 1:100) {
  P_1[i,] = Can_I_make_money(1,1, R, p)
}
sum(with(P_1, P>0)) ## numbers of tries with a gain of a% log return
mean(P_1$Time) ##Mean of stopping times



## Try 100 times with a=b=10%
P_10 = data.frame(P = double(), Time= double())
for (i in 1:100) {
  P_10[i,] = Can_I_make_money(10,10, R, p)
}
sum(with(P_10, P>0))
mean(P_10$Time)


## Try 100 times with a=10%, b=1%
P_a = data.frame(P = double(), Time= double())
for (i in 1:100) {
  P_a[i,] = Can_I_make_money(10,1, R, p)
}
sum(with(P_a, P>0))
mean(P_a$Time)



## Try 100 times with a=10%, b=0.1%
P_b = data.frame(P = double(), Time= double())
for (i in 1:100) {
  P_b[i,] = Can_I_make_money(10,0.1, R, p)
}
sum(with(P_b, P>0))
mean(P_b$Time)
