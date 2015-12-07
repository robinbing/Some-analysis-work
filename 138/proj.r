# load data
setwd("D:/fall 2015/sta138/final project")
df = read.table('final.dat', header = TRUE)

length(which(df$dav==1))

# VIF
diag(solve(cor(df[,-1])))

# goodness of fit test
library('ResourceSelection')
sat = glm(dav~ .^2, family = binomial, data = df)
hoslem.test(df$dav, sat$fitted.values)

# model selection
backward = step(sat, direction = 'both', trace = FALSE)

# goodness of fit
model2 = glm(backward, family = binomial, data = df)
anova(sat, model2, test='Chisq')

# goodness-of-link
z = -(1 + 1/model2$fitted.values * log(1 - model2$fitted.values))
model3 = glm(formula = dav ~ mcs + beck + pgend + age + educat + z, 
             family = binomial, data = df)
summary(model3)

# selec model2
summary(model2)

# sd pearson resid
peadRed = residuals(model2, 'pearson')/(1-influence(model2)$hat)
plot(peadRed, ylab='Standardized Pearson Resuduals')
abline(h=2, col='red')

# outlier
outlying.y = which(peadRed>2)
outlying.y

outlying.x = which(influence(model2)$hat > 2*6/400)
outlying.x

# influential
outlier = union(outlying.x,outlying.y)
reg = function(x){
  glm(formula = dav ~ mcs + beck + pgend + age + educat, family = binomial, data=df[-x,])
}
Diff = sapply(outlier, function(x) mean(abs(reg(x)$fitted - model2$fitted[-x])))
summary(Diff)
