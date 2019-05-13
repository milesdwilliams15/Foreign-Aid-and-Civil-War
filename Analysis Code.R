library(haven)
mig_data <- read_dta("finaldata.dta")

# Estimate ML Tobit models
library(dplyr)
library(survival)
library(lmtest)
library(ggplot2)
library(dotwhisker)
library(gridExtra)
library(stargazer)
mig_data <- mig_data %>% group_by(dyad) %>%
  mutate(istock_lag=lag(istock,order.by='year')) %>%
  ungroup() %>%
  mutate(istock_change=istock-istock_lag) %>%
  group_by(dyad) %>%
  mutate(istock_change_lag=lag(istock_change,order.by='year'),
         civilwar_lag=lag(civilwar,order.by='year'),
         disaster_lag=lag(disaster,order.by='year'),
         usmil_lag=lag(usmil,order.by='year'),
         gdpcap_lag=lag(gdpcap,order.by='year'),
         pop_lag=lag(population,order.by='year'),
         fh_lag=lag(fh,order.by='year'),
         trade_lag=lag(trade,order.by='year'),
         dualcitizenship_lag=lag(dualcitizenship,order.by='year'),
         any_votingrc_lag=lag(any_votingrc,order.by='year'))
mig_data$istock_change[mig_data$istock_change<0] <- 0
mig_data$istock_change_lag[mig_data$istock_change_lag<0] <- 0
mltob1 <- survreg(Surv(log(commit3a+1),log(commit3a+1)>0,type='left') ~ 
                    civilwar_lag + log(istock_lag+1) + 
                    log(trade_lag+1) +
                    log(usmil_lag+1) +
                    log(distance) +
                    colony + 
                    log(disaster_lag+1) + 
                    log(gdpcap_lag) + log(pop_lag) + fh_lag +
                    donor + as.factor(year) +
                    frailty.gaussian(dyad),
                  mig_data %>% 
                    filter(year>= 1994 & donor != 'US'),
                  dist='gaussian')

broom::tidy(coeftest(mltob1))[2:11,] %>% 
  mutate(term=c('Civil War','Migration (ln)',
                'Trade (ln)','US Military Aid (ln)',
                'Distance (ln)','Colony','Disaster (ln)',
                'Income (ln)','Population (ln)','Democracy')) %>%
  dwplot(dot_args=list(color='black')) + geom_vline(xintercept=0,linetype=2) +
  xlab('Estimated Coefficient\n(95% Confidence Intervals Shown)') +
  theme_bw() +
  theme(legend.position = 'none',
        text=element_text(family='serif'),
        axis.text = element_text(color='black')) 

cw_b <- c()
cw_se <- c()
mig_vals <- with(mig_data %>% filter(year>=1994),
                 seq(round(min(log(istock_lag+1),na.rm=T)),
                     round(max(log(istock_lag+1),na.rm=T)),len=7))
for(i in 1+mig_vals){
  mltob <- survreg(Surv(log(commit3a+1),log(commit3a+1)>0,type='left') ~ 
                     civilwar_lag*lmig+log(distance) + log(usmil_lag+1) +
                     colony + log(disaster_lag+1) + fh_lag +
                     log(gdpcap_lag) + log(pop_lag) + 
                     log(trade_lag+1) + 
                     donor+as.factor(year) + 
                     frailty.gaussian(dyad),
                   mig_data %>% filter(year>= 1994) %>%
                     mutate(lmig=log(istock_lag+1)-(i-1)),
                   dist='gaussian')
  smry <- coeftest(mltob)
  cw_b[i] <- smry[2,1]
  cw_se[i] <- smry[2,2]
}
cw_b2 <- c()
cw_se2 <- c()
trade_vals <- with(mig_data %>% filter(year>=1994),
                   seq(round(min(log(trade_lag+1),na.rm=T)),
                       round(max(log(trade_lag+1),na.rm=T)),len=7))
for(i in 1+trade_vals){
  mltob <- survreg(Surv(log(commit3a+1),log(commit3a+1)>0,type='left') ~ 
                     civilwar_lag*ltrade+log(istock_lag+1)+log(distance) + 
                     log(usmil_lag+1) +
                     colony + log(disaster_lag+1) + fh_lag +
                     log(gdpcap_lag) + log(pop_lag) + 
                     donor+as.factor(year) + 
                     frailty.gaussian(dyad),
                   mig_data %>% filter(year>= 1994) %>%
                     mutate(ltrade=log(trade_lag+1)-(i-1)),
                   dist='gaussian')
  smry <- coeftest(mltob)
  cw_b2[i] <- smry[2,1]
  cw_se2[i] <- smry[2,2]
}
cw_b3 <- c()
cw_se3 <- c()
usmil_vals <- with(mig_data %>% filter(year>=1994),
                   seq(round(min(log(usmil_lag+1),na.rm=T)),
                       round(max(log(usmil_lag+1),na.rm=T)),len=7))
for(i in 1+usmil_vals){
  mltob <- survreg(Surv(log(commit3a+1),log(commit3a+1)>0,type='left') ~ 
                     civilwar_lag*lusmil+log(istock_lag+1)+log(distance) + 
                     log(trade_lag+1) +
                     colony + log(disaster_lag+1) + fh_lag +
                     log(gdpcap_lag) + log(pop_lag) + 
                     donor+as.factor(year) + 
                     frailty.gaussian(dyad),
                   mig_data %>% filter(year>= 1994) %>%
                     mutate(lusmil=log(usmil_lag+1)-(i-1)),
                   dist='gaussian')
  smry <- coeftest(mltob)
  cw_b3[i] <- smry[2,1]
  cw_se3[i] <- smry[2,2]
}
cw_b4a <- c()
cw_se4a <- c()
cw_b4b <- c()
cw_se4b <- c()
mig_vals2 <- with(mig_data %>% filter(year>=1994),
                  seq(round(min(log(istock_lag+1),na.rm=T)),
                      round(max(log(istock_lag+1),na.rm=T)),len=3))
for(i in 1+mig_vals2){
  mltob1 <- survreg(Surv(log(commit3a+1),log(commit3a+1)>0,type='left') ~ 
                      civilwar_lag*lmig*I(dualcitizenship_lag*any_votingrc_lag)+
                      log(distance) + log(usmil_lag+1) +
                      colony + log(disaster_lag+1) + fh_lag +
                      log(gdpcap_lag) + log(pop_lag) + 
                      log(trade_lag+1) + 
                      donor+as.factor(year) + 
                      frailty.gaussian(dyad),
                    mig_data %>% filter(year>= 1994) %>%
                      mutate(lmig=log(istock_lag+1)-(i-1)),
                    dist='gaussian')
  smry1 <- coeftest(mltob1)
  cw_b4a[i] <- smry1[2,1]
  cw_se4a[i] <- smry1[2,2]
  mltob2 <- survreg(Surv(log(commit3a+1),log(commit3a+1)>0,type='left') ~ 
                      civilwar_lag*lmig*I(dualcitizenship_lag*any_votingrc_lag-1)+
                      log(distance) + log(usmil+1) +
                      colony + log(disaster_lag+1) + fh_lag +
                      log(gdpcap_lag) + log(pop_lag) + 
                      log(trade_lag+1) + 
                      donor+as.factor(year) + 
                      frailty.gaussian(dyad),
                    mig_data %>% filter(year>= 1994) %>%
                      mutate(lmig=log(istock_lag+1)-(i-1)),
                    dist='gaussian')
  smry2 <- coeftest(mltob2)
  cw_b4b[i] <- smry2[2,1]
  cw_se4b[i] <- smry2[2,2]
}
cw_b5a <- c()
cw_se5a <- c()
cw_b5b <- c()
cw_se5b <- c()
for(i in 1+mig_vals2){
  mltob1 <- survreg(Surv(log(commit3a+1),log(commit3a+1)>0,type='left') ~ 
                      civilwar_lag*lmig*dualcitizenship_lag+
                      any_votingrc_lag +
                      log(distance) + log(usmil_lag+1) +
                      colony + log(disaster_lag+1) + fh_lag +
                      log(gdpcap_lag) + log(pop_lag) + 
                      log(trade_lag+1) + 
                      donor+as.factor(year) + 
                      frailty.gaussian(dyad),
                    mig_data %>% filter(year>= 1994) %>%
                      mutate(lmig=log(istock_lag+1)-(i-1)),
                    dist='gaussian')
  smry1 <- coeftest(mltob1)
  cw_b5a[i] <- smry1[2,1]
  cw_se5a[i] <- smry1[2,2]
  mltob2 <- survreg(Surv(log(commit3a+1),log(commit3a+1)>0,type='left') ~ 
                      civilwar_lag*lmig*I(dualcitizenship_lag-1)+
                      any_votingrc_lag +
                      log(distance) + log(usmil+1) +
                      colony + log(disaster_lag+1) + fh_lag +
                      log(gdpcap_lag) + log(pop_lag) + 
                      log(trade_lag+1) + 
                      donor+as.factor(year) + 
                      frailty.gaussian(dyad),
                    mig_data %>% filter(year>= 1994) %>%
                      mutate(lmig=log(istock_lag+1)-(i-1)),
                    dist='gaussian')
  smry2 <- coeftest(mltob2)
  cw_b5b[i] <- smry2[2,1]
  cw_se5b[i] <- smry2[2,2]
}
cw_b6a <- c()
cw_se6a <- c()
cw_b6b <- c()
cw_se6b <- c()
for(i in 1+mig_vals2){
  mltob1 <- survreg(Surv(log(commit3a+1),log(commit3a+1)>0,type='left') ~ 
                      civilwar_lag*lmig*any_votingrc_lag+
                      dualcitizenship_lag +
                      log(distance) + log(usmil_lag+1) +
                      colony + log(disaster_lag+1) + fh_lag +
                      log(gdpcap_lag) + log(pop_lag) + 
                      log(trade_lag+1) + 
                      donor+as.factor(year) + 
                      frailty.gaussian(dyad),
                    mig_data %>% filter(year>= 1994) %>%
                      mutate(lmig=log(istock_lag+1)-(i-1)),
                    dist='gaussian')
  smry1 <- coeftest(mltob1)
  cw_b6a[i] <- smry1[2,1]
  cw_se6a[i] <- smry1[2,2]
  mltob2 <- survreg(Surv(log(commit3a+1),log(commit3a+1)>0,type='left') ~ 
                      civilwar_lag*lmig*I(any_votingrc_lag-1)+
                      dualcitizenship_lag +
                      log(distance) + log(usmil+1) +
                      colony + log(disaster_lag+1) + fh_lag +
                      log(gdpcap_lag) + log(pop_lag) + 
                      log(trade_lag+1) + 
                      donor+as.factor(year) + 
                      frailty.gaussian(dyad),
                    mig_data %>% filter(year>= 1994) %>%
                      mutate(lmig=log(istock_lag+1)-(i-1)),
                    dist='gaussian')
  smry2 <- coeftest(mltob2)
  cw_b6b[i] <- smry2[2,1]
  cw_se6b[i] <- smry2[2,2]
}

na.omit(data.frame(cw_b=100*(exp(cw_b)-1),cw_se=100*(exp(cw_se)-1))) %>%
  mutate(mig_level=mig_vals) %>%
  ggplot(aes(mig_level,cw_b)) + geom_point() + geom_line() +
  geom_errorbar(aes(ymin=cw_b-1.96*cw_se,ymax=cw_b+1.96*cw_se),
                width=.3) +
  geom_hline(yintercept=0,linetype=2) +
  labs(x="Migrant Stock (ln)",
       y="% Change in Expected Bilateral Aid \nConditional on Recipient Civil War\n(95% C.I.s shown)",
       title="Migration") +
  theme_bw() +
  theme(plot.title = element_text(hjust=.5),
        panel.grid.minor = element_blank(),
        text=element_text(family='serif')) +
  scale_x_continuous(breaks=mig_vals,labels=round(mig_vals,2))

mig_data %>% 
  filter(year>=1994) %>% ungroup() %>%
  mutate("Civil War"=car::recode(civilwar_lag,"1='Yes';0='No';else=NA")) %>%
  filter(`Civil War`!='NA') %>%
  ggplot(aes(log(istock_lag),log(commit3a),linetype=`Civil War`)) + 
  geom_smooth(method='lm',alpha=.5) +
  theme_bw() +
  labs(x='Bilateral Migration (ln)',
       y='Bilateral Aid (ln)') +
  facet_wrap(~donor,scales='free')

na.omit(data.frame(cw_b2=100*(exp(cw_b2)-1),cw_se2=100*(exp(cw_se2)-1))) %>%
  mutate(trade_level=trade_vals) %>%
  ggplot(aes(trade_level,cw_b2)) + geom_point() + geom_line() +
  geom_errorbar(aes(ymin=cw_b2-1.96*cw_se2,ymax=cw_b2+1.96*cw_se2),
                width=.3) +
  geom_hline(yintercept=0,linetype=2) +
  labs(x="Bilateral Trade (ln)",
       y="% Change in Expected Bilateral Aid \nConditional on Recipient Civil War\n(95% C.I.s shown)",
       title="Trade") +
  theme_bw() +
  theme(plot.title=element_text(hjust=.5),
        panel.grid.minor = element_blank(),
        text=element_text(family='serif')) +
  scale_x_continuous(breaks=trade_vals,
                     labels=round(trade_vals,3))

mig_data %>% 
  filter(year>=1994) %>% ungroup() %>%
  mutate("Civil War"=car::recode(civilwar_lag,"1='Yes';0='No';else=NA")) %>%
  filter(`Civil War`!='NA') %>%
  ggplot(aes(log(trade_lag),log(commit3a),linetype=`Civil War`)) + 
  geom_smooth(method='lm',alpha=.5) +
  theme_bw() +
  labs(x='Bilateral Trade (ln)\n(Yearly Average per Dyad)',
       y='Bilateral Aid (ln)\n(Yearly Average per Dyad)') +
  facet_wrap(~donor,scales='free')

na.omit(data.frame(cw_b3=100*(exp(cw_b3)-1),
                   cw_se3=100*(exp(cw_se3)-1))) %>%
  mutate(usmil_level=usmil_vals) %>%
  ggplot(aes(usmil_level,cw_b3)) + geom_point() + geom_line() +
  geom_errorbar(aes(ymin=cw_b3-1.96*cw_se3,ymax=cw_b3+1.96*cw_se3),
                width=.3) +
  geom_hline(yintercept=0,linetype=2) +
  labs(x="US Military Aid (ln)",
       y="% Change in Expected Bilateral Aid \nConditional on Recipient Civil War\n(95% C.I.s Shown)",
       title="Strategic Military Value") +
  theme_bw() +
  theme(plot.title=element_text(hjust=.5),
        panel.grid.minor = element_blank(),
        text=element_text(family='serif')) +
  scale_x_continuous(breaks=usmil_vals,labels=round(usmil_vals,2))

na.omit(data.frame(cw_b=c(cw_b4a,cw_b4b),cw_se=c(cw_se4a,cw_se4b))) %>%
  mutate(mig_level=c(mig_vals2,mig_vals2),
         Mobilization=c(rep(c('no','yes'),each=length(mig_vals2)))) %>%
  ggplot(aes(mig_level,cw_b,linetype=Mobilization)) + geom_point(position = position_dodge(1)) + 
  geom_line(position=position_dodge(1)) +
  geom_errorbar(aes(ymin=cw_b-1.96*cw_se,ymax=cw_b+1.96*cw_se),
                width=.3,
                position=position_dodge(1)) +
  geom_hline(yintercept=0,linetype=2) +
  labs(x="Migrant Stock (ln)",
       y="Marginal Effect of Civil War\n(95% Confidence Intervals Shown)",
       title="Migrant Mobilization") +
  theme_bw() +
  theme(plot.title = element_text(hjust=.5),
        panel.grid.minor.x = element_blank(),
        legend.position = c(.2,.9)) +
  scale_x_continuous(breaks=mig_vals2,labels=round(mig_vals2,2))
na.omit(data.frame(cw_b=c(cw_b5a,cw_b5b),cw_se=c(cw_se5a,cw_se5b))) %>%
  mutate(mig_level=c(mig_vals2,mig_vals2),
         'Dual Citizenship'=c(rep(c('no',
                                    'yes'),each=length(mig_vals2)))) %>%
  ggplot(aes(mig_level,cw_b,linetype=`Dual Citizenship`)) + geom_point(position = position_dodge(1)) + 
  geom_line(position=position_dodge(1)) +
  geom_errorbar(aes(ymin=cw_b-1.96*cw_se,ymax=cw_b+1.96*cw_se),
                width=.3,
                position=position_dodge(1)) +
  geom_hline(yintercept=0,linetype=2) +
  labs(x="Migrant Stock (ln)",
       y="Marginal Effect of Civil War\n(95% Confidence Intervals Shown)",
       title="Migrant Mobilization") +
  theme_bw() +
  theme(plot.title = element_text(hjust=.5),
        panel.grid.minor.x = element_blank(),
        legend.position = c(.2,.9)) +
  scale_x_continuous(breaks=mig_vals2,labels=round(mig_vals2,2))
na.omit(data.frame(cw_b=c(cw_b6a,cw_b6b),cw_se=c(cw_se6a,cw_se6b))) %>%
  mutate(mig_level=c(mig_vals2,mig_vals2),
         'Voting Rights'=c(rep(c('no',
                                 'yes'),each=length(mig_vals2)))) %>%
  ggplot(aes(mig_level,cw_b,linetype=`Voting Rights`)) + geom_point(position = position_dodge(1)) + 
  geom_line(position=position_dodge(1)) +
  geom_errorbar(aes(ymin=cw_b-1.96*cw_se,ymax=cw_b+1.96*cw_se),
                width=.3,
                position=position_dodge(1)) +
  geom_hline(yintercept=0,linetype=2) +
  labs(x="Migrant Stock (ln)",
       y="Marginal Effect of Civil War\n(95% Confidence Intervals Shown)",
       title="Migrant Mobilization") +
  theme_bw() +
  theme(plot.title = element_text(hjust=.5),
        panel.grid.minor.x = element_blank(),
        legend.position = c(.2,.9)) +
  scale_x_continuous(breaks=mig_vals2,labels=round(mig_vals2,2))

mf <- model.frame(lm(commit3a ~ 
                       civilwar_lag + istock_lag + 
                       trade_lag +
                       usmil_lag +
                       distance +
                       colony + 
                       disaster_lag + 
                       gdpcap_lag + pop_lag + fh_lag,
                     mig_data %>% filter(year>=1994)))
stargazer(mf,header=F,
          title="Summary Statistics (1994-2008)",
          covariate.labels = c("Aid Commitments\n(in thousands)",
                               "Civil War",
                               "Migrants",
                               "Trade",
                               "US Military Aid",
                               "Distance\n(kilometers)",
                               "Colony",
                               "Disaster",
                               "Income",
                               "Population",
                               "Democracy\n(Freedom House)"),
          column.separate = 0)

donors <- matrix(sort(unique(mig_data$donor)),
                 ncol=2)
stargazer(donors,header=F,
          title="OECD Donor Countries")

recipients <- matrix(sort(unique(mig_data$recipient)),
                     ncol=4)
stargazer(recipients,header=F,
          title="Recipient Countries")
