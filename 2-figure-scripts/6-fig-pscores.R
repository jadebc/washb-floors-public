# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

# load main data
bd <- readRDS(clean_bdata_path)
ke <- readRDS(clean_kdata_path)

# load propensity scores
bpred = readRDS(bpred_path)
bpred$dataid = as.character(bpred$dataid)

kpred = readRDS(kpred_path)

# merge propensity scores into main data
b = full_join(bd, bpred, by = "dataid")
k = full_join(ke, kpred, by = "hhid")

# bright color blind palette:  https://personal.sron.nl/~pault/ 
cblack <- "#000004FF"
cblue <- "#3366AA"
cteal <- "#11AA99"
cgreen <- "#66AA55"
cchartr <- "#CCCC55"
cmagent <- "#992288"
cred <- "#EE3333"
corange <- "#EEA722"
cyellow <- "#FFEE33"
cgrey <- "#777777"

#----------------------------------
# Plot distributions of P(floor|W)
# Bangladesh
#----------------------------------
btrim_lower = min(b$pred, na.rm=TRUE) + 0.01
btrim_upper = max(b$pred, na.rm=TRUE) - 0.01

pdf(paste0(fig_path, "/wbb-Pfloor-dists.pdf"),height=5,width=10)
lo <- graphics::layout(mat=matrix(1:2,nrow=1,ncol=2))
hist(b$pred[b$floor==0],
     breaks=seq(0,1,by=0.01),
     col=alpha(cblack,alpha=0.3),
     border=NA,
     ylim=c(0,2800),
     xlim=c(0,1),
     bty='n',
     las=1,xlab='P(improved floor|W)',main='Bangladesh')

hist(b$pred[b$floor==1],
     xlim=c(0,1),
     breaks=seq(0,1,by=0.01),
     col=alpha(cblue,alpha=0.5),
     border=NA,
     add=TRUE)

abline(v=btrim_lower, lty = "dashed")
abline(v=btrim_upper, lty = "dashed")

legend("top",legend=c("Soil floor","Improved floor"),fill=c(alpha(cblack,alpha=0.3),alpha(cblue,alpha=0.5)),bty='n')

hist(log(b$pred[b$floor==0]),
     breaks=100,
     col=alpha(cblack,alpha=0.3),
     border=NA,
     # ylim=c(0,2800),
     bty='n',
     las=1,xlab='log P(improved floor|W)',main='')

hist(log(b$pred[b$floor==1]),
     breaks=100,
     col=alpha(cblue,alpha=0.5),
     border=NA,
     add=TRUE)

abline(v=log(btrim_lower), lty = "dashed")
abline(v=log(btrim_upper), lty = "dashed")

graphics::layout(mat=1)
dev.off()



#----------------------------------
# Plot distributions of P(floor|W)
# Kenya
#----------------------------------
ktrim_lower = min(k$pred, na.rm=TRUE) + 0.01
ktrim_upper = max(k$pred, na.rm=TRUE) - 0.01

pdf(paste0(fig_path, "/wbk-Pfloor-dists.pdf"),height=5,width=10)
lo <- graphics::layout(mat=matrix(1:2,nrow=1,ncol=2))
hist(k$pred[k$floor==0],
     breaks=seq(0,1,by=0.01),
     col=alpha(cblack,alpha=0.3),
     border=NA,
     # ylim=c(0,200),
     xlim=c(0,1),
     bty='n',
     las=1,xlab='P(improved floor|W)',main='Kenya')

hist(k$pred[k$floor==1],
     xlim=c(0,1),
     breaks=seq(0,1,by=0.01),
     col=alpha(cblue,alpha=0.5),
     border=NA,
     add=TRUE)

abline(v=ktrim_lower, lty = "dashed")
abline(v=ktrim_upper, lty = "dashed")

legend("topright",legend=c("Soil floor","Improved floor"),fill=c(alpha(cblack,alpha=0.3),alpha(cblue,alpha=0.5)),bty='n')

hist(log(k$pred[k$floor==0]),
     breaks=100,
     col=alpha(cblack,alpha=0.3),
     border=NA,
     # ylim=c(0,200),
     bty='n',
     las=1,xlab='log P(improved floor|W)',main='')

hist(log(k$pred[k$floor==1]),
     breaks=100,
     col=alpha(cblue,alpha=0.5),
     border=NA,
     add=TRUE)

abline(v=log(ktrim_lower), lty = "dashed")
abline(v=log(ktrim_upper), lty = "dashed")

graphics::layout(mat=1)
dev.off()

