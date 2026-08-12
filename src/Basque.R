
rm(list = ls())
source("~/Documents/GitHub/DRoSC/src/helpers.R")

library(ggcorrplot)
library(Synth)
### Basque ###
### T0 = 15, N = 16, T1 = 28
data("basque") 

T0 <- 15
T1 <- 28
Y <- basque[basque$regionname=="Basque Country (Pais Vasco)",]$gdpcap
X <- matrix(basque[basque$regionname!="Basque Country (Pais Vasco)" & basque$regionname!="Spain (Espana)",]$gdpcap,
            nrow = T0+T1,byrow = F)
colnames(X) <- unique(basque[basque$regionname!="Basque Country (Pais Vasco)"& basque$regionname!="Spain (Espana)",]$regionname)
colnames(X) <-c("Andalucia","Aragon","Principado De Asturias","Baleares","Canarias","Cantabria",
                "Castilla Y Leon","Castilla-La Mancha","Cataluna","Comunidad Valenciana","Extremadura","Galicia",
                "Madrid","Murcia","Navarra","Rioja"
                
)


Y0 <- Y[1:T0]
Y1 <- Y[-(1:T0)]
X0 <- X[1:T0,]
X1 <- X[-(1:T0),]

N <- ncol(X)
# standard SC #
SC.fit <- sc(Y0,X0)
SC.beta <- SC.fit$w.hat
tau.t <- Y1-X1%*%SC.beta
tau.SC<- mean(Y1-X1%*%SC.beta)
# Correlation plot (Figure 1) #
SC.names <- names(SC.beta[which(SC.beta>0)][order(SC.beta[which(SC.beta>0)],decreasing = F)])
cor.plot <- ggcorrplot(cor(X0)[,SC.names])+ theme( legend.text = element_text(size = 14),text = element_text(size = 15))
cor.plot
mu <- apply(X1,2,mean)

l.cand <- seq(0,0.06,0.001)
dt.pe<- data.frame()
for (i in 1:length(l.cand)) {
  set.seed(1)
  l <- l.cand[i]
  DRSC.fit <- DRoSC(Y0,Y1,X0,X1,lambda = l,Inference = T)
  
  dt.pe <- data.frame(lambda = l,tauHat = DRSC.fit$tauHat,
                      CI.tau.1 = DRSC.fit$CI.tau[,1],
                      CI.tau.2 = DRSC.fit$CI.tau[,2])
}
dt.all <- dt.pe


offset <- 0.00 # Shift amount for side-by-side plotting
lambda.plot <- ggplot(dt.all[dt.all$lambda %in% c(0,0.015,0.03,0.045,0.06),], aes(lambda))+
  geom_errorbar(aes(x = lambda-offset, ymin = CI.tau.1, ymax = CI.tau.2), 
                width = 0.001, linewidth = 0.4, color = "orange")+
  geom_line(aes(y = tauHat),colour = 'black')+
  geom_hline(yintercept = 0,colour = 'black',linetype = 'dashed')+
  geom_hline(yintercept = tau.SC,colour = 'blue',linetype = 'dashed')+
  geom_vline(xintercept = dt.all$lambda[min(which(round(dt.all$tauHat,4)==0))],colour = 'red',linetype = 'dashed')+
  xlab(latex2exp::TeX("$\\lambda$"))+
  ylab(latex2exp::TeX("$\\hat{\\tau}$"))+
  theme_minimal()+
  theme(text = element_text(size = 17)) + guides(fill="none")
lambda.plot

apply(X1,2,mean)-apply(X0,2,mean)
mean(Y1)-mean(Y0)
log(3.309/0.01,base = 1.25)

Y0
X0
X1

source("~/Documents/GitHub/DRoSC/src/helpers.R")
set.seed(1)
DRSC.fit <- DRoSC(Y0,Y1,X0,X1,lambda = 0.0,Inference = T,M=500)
DRSC.fit$lambda.pert
df <- data.frame(
  id = c(0:(nrow(DRSC.fit$Int.mat)-1),nrow(DRSC.fit$Int.mat)+2),
  lower = c(DRSC.fit$Int.mat[, 1],DRSC.fit$CI.tau[,1]),
  upper = c(DRSC.fit$Int.mat[, 2],DRSC.fit$CI.tau[,2])
)

df$agg <- df$id == max(df$id)
ggplot(df, aes(x = id, ymin = lower, ymax = upper, color = agg)) +
  geom_linerange(size = 0.7) +
  geom_hline(yintercept = 0,colour = 'darkgreen',linetype = 'dashed')+
  geom_hline(yintercept = DRSC.fit$tauHat,colour = 'mediumblue',linetype = 'longdash')+
  scale_color_manual(
    values = c("FALSE" = "black", "TRUE" = "red"),
    labels = c("Individual intervals", "Aggregated interval"),
    name = NULL
  ) +
  labs(x = "Perturbation index", y = "Perturbed intervals") +
  theme_minimal()+ theme(legend.position = "none")

idx <- which(apply(DRSC.fit$Int.mat,1,mean)<0.9&apply(DRSC.fit$Int.mat,1,mean)>-0.9)

DRSC.fit$Int.mat[-idx,]
apply(DRSC.fit$mu.valid[-idx,],1,function(x){x-mu})
DRSC.fit$beta.valid[-idx,]
DRSC.fit$betaHat
rowSums(DRSC.fit$mu.valid[-idx,]*DRSC.fit$beta.valid[-idx,])-sum(mu*DRSC.fit$betaHat)

