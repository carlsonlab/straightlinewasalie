
coef <- function(H,h,k,khat) {choose(h,khat)*choose(H-h,k-khat)/choose(H,k)}

coef(6000,600,10,0)/sum(sapply(c(0:10),function(i){coef(6000,600,10,i)}))
coef(6000,600,10,1)/sum(sapply(c(1:10),function(i){coef(6000,600,10,i)}))

coef(6000,60,10,0)/sum(sapply(c(0:10),function(i){coef(6000,60,10,i)}))
coef(6000,60,10,1)/sum(sapply(c(1:10),function(i){coef(6000,60,10,i)}))

ratio <- function(H,h,k,khat) {(H-h-k+1)/(h*k)}

ratio(6000,60,10,1)
