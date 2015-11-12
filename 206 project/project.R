setwd('e:\\STA206\\assignment\\project')
abalone = read.table("abalone.txt" , sep =',')
colnames(abalone) = c("sex","length","diameter","height","whole","sucked",
                     "viscera","shell","rings")

a = abalone
a$sex = factor(a$sex)
#
#split data 
#
set.seed(12)
n = nrow(a)
index.s = sample(1:n , size = n*2/3 , replace = FALSE)

a.train = a[index.s,]
a.test = a[-index.s,]

sapply(2:9, function(i) boxplot(a.train[,i],a.test[,i]))
###############################################################

boxplot(rings~sex,data = a.train,
        main='side-by-side boxplots',xlab='factor levels',
        ylab='observation',col=rainbow(6))


#############################################################################
#############################################################################
#transform response variables
a.train$rings = log(a.train$rings)
a.test$rings = log(a.test$rings)
#
#model selection
#
#
#selection of first-order effects
#
#
model.1st = lm(rings~. , data = a.train)
mse = summary(model.1st)$sigma^2

#
#best subsets selection
#
library(leaps)
sub.set = regsubsets(rings~. , data = a.train , nbest = 1, 
                     nvmax = 16 , method = "exhaustive")
sum.sub = summary(sub.set)

#number of parameters in each model
num.p = as.numeric(rownames(sum.sub$which)) + 1L

#parameters in model 
n.train = nrow(a.train)
sse = sum.sub$rss

#aic ,  pic
aic = n.train*log(sse/n) + 2*num.p
bic = n.train*log(sse/n) + log(n)*num.p

sub.table = cbind(sum.sub$which, sse, sum.sub$rsq, sum.sub$adjr2, 
                  sum.sub$cp, aic ,bic)

#null model
fit0 = lm(rings~1, data = a.train)
sse0 = sum(fit0$residuals^2)
p0 = 1
c0 = sse0/mse - (n.train-2*p0)
aic0 = n.train*log(sse0/n.train) + 2*p0
bic0 = n.train*log(sse0/n.train) +log(n.train)*p0
none = c(p0, rep(0,9), sse0, 0, 0, c0, aic0, bic0)

sub.table = rbind(none, sub.table)
colnames(sub.table) = c(colnames(sum.sub$which), "sse", "R^2", "R^2_a", "cp",
                        "aic", "bic")

#
#forward stepwise procedure
#
library(MASS)

step.forward = stepAIC(fit0, scope = list(upper = model.1st, lower = ~1),
                       direction = "both", k=2)
#
#
#selection of first-order and second-order effects
#
#
model.2nd = lm(rings~.^2, data = a.train)
mse2 = summary(model.2nd)$sigma^2
#
#forward stepwise procedure
#
step.forward2 = stepAIC(fit0, scope = list(upper = model.2nd, lower = ~1),
                        direction = "both", k=2)
###############################################
###############################################
#
#model validation
#
#
#
#internal validation
#
model1 = lm(step.forward , data = a.train)
plot(model1, which = 1)
plot(model1, which = 2)

model2 = lm(step.forward2 , data = a.train)
plot(model, which = 1)
plot(model, which = 2)

sse.1st = anova(step.forward)["Residual" , 2]
p.1st = length(step.forward$coefficients)
cp.1st = sse.1st/mse2 - (n.train-2*p.2nd)
press.1st = sum(step.forward$residuals^2/(1-influence(step.forward)$hat)^2)
mse.1st = anova(step.forward)["Residuals",3]
#cp太大

sse.2nd = anova(step.forward2)["Residual" , 2]
p.2nd = length(step.forward2$coefficients)
cp.2nd = sse.2nd/mse2 - (n.train-2*p.2nd)
press.2nd = sum(step.forward2$residuals^2/(1-influence(step.forward2)$hat)^2)
mse.2nd = anova(step.forward2)["Residuals",3]
#(cp是51，和p 24差了些， 可能是在一开始统计数据的时候就少了些重要变量造成了model bias)
#press.2nd/n = 0.00733 , mse.2nd = 0.00706. Little difference between these two variables
#supports the validity of the model. And the mse is small which shows a good ablity of 
#the model

#
#external validation
#
#caculation
model2.v = lm(step.forward2 , data = a.test)

mspr2 =round (mean((predict.lm(model2, a.test)-a.test$rings)^2),3)

press.2nd/n.train

sse_model2 = round(anova(model2)["Residuals",2],3)

sse_model2.v = round(anova(model2.v)["Residuals",2],3)

mse_model2 = round(anova(model2)["Residuals",3],3)

mse_model2.v = round(anova(model2.v)["Residuals",3],3)

model2_R2_a = round(summary(model2)$adj.r.squared,3)

model2_R2_a.v = round(summary(model2.v)$adj.r.squared,3)
#model2
mod_sum_2 = cbind(coef(summary(model2.v))[,1], coef(summary(model2.v))[,2], 
                  coef(summary(model2))[,1],coef(summary(model2))[,2])
colnames(mod_sum_2) = c('coef validation','coef std.err validation',
                        'coef ','coef std.err')

Training_2 = cbind(sse_model2,mse_model2,model2_R2_a,round(press.2nd,3),
                   round(press.2nd/n.train,3),"--")
Validation_2 = cbind(sse_model2.v,mse_model2.v,model2_R2_a.v,"--","--",
                     mspr2)
con_2 = rbind(Training_2,Validation_2)
rownames(con_2) = c('Training','Validation')
colnames(con_2) = c('sse','mse','R2_2','press','press/n','mspr')

#下面这两个是table
mod_sum_2
con_2

#################################################################################
#################################################################################
#
#outlying
#
#outlying y
model.final = lm(step.forward2, data = a)
hii = influence(model.final)$hat
mse = anova(model.final)["Residuals",3]
res = model.final$residuals
stu.res = res/sqrt(mse*(1-hii))   #studentized residuals

res.del = res / (1-hii)   # deleted residuals
library(MASS)
stu.res.del = studres(model.final)  #studentized deleted residuals
bon.thre = qt(1-0.1/(2*n),n-model.final$rank-1)

#residuals vs. fitted values plots
plot(model.final$fitted, stu.res.del , xlab="fitted value", ylab="residual",
     cex.lab=1.5, cex.axis = 1.5, pch = 19, cex = 1.5)
abline(h=0, col = grey(0.8), lwd = 2, lty = 2)
abline(h = bon.thre, lwd = 2, lty = 3)
abline(h = -bon.thre, lwd = 2, lty = 3)

#test for outlying Y
sse = sum((summary(model.final)$residuals)^2)
ti = res*sqrt((nrow(a)-fit$rank-1)/(sse*(1-hii)-res^2))
tt = qt(1-0.1/(2*nrow(a)) , nrow(a)-fit$rank-1 )
any(abs(ti)>tt)
index_outy = which(abs(ti)>tt)


#test for outlying X
any(hii>2*model.final$rank/nrow(a))
index_outx = which(hii>2*model.final$rank/nrow(a))

#cook's distance (outlying influence)
Di = stu.res^2*hii/(model.final$rank*(1-hii))
plot(Di,type="h",ylab = "Cook's distance")

Di = c(Di)
dd = pf(Di , model.final$rank, nrow(a)-model.final$rank)
any(dd>0.5)

#DFFITS DFBETAS
sta = influence.measures(model.final)

#DFFITS
2*sqrt(model.final$rank/n)

#DFBETAS
2/sqrt(n)


