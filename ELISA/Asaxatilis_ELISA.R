############################################################################
#################### ELISA Analysis of Oxidative Stress ####################
####################### Abudefduf Saxatilis - liver ########################
############################################################################

setwd("C:/Users/allys/Box/Auburn/Scientific Writing/Analyses and Data/ELISA")

library(ggplot2)
library(car)
library(nlme)

elisa <- as.data.frame(read.csv("A.saxatilis_ELISA-w-attributes.csv"))
elisa$Tower = as.factor(elisa$Tower)

plot(data = elisa, Protein.Carbonyl~Protein.Raw)
check = lm(data = elisa, Protein.Carbonyl~Protein.Raw)
summary(check)
abline(check)

plot(data = elisa, conc~Protein.Raw)
check = lm(data = elisa, conc~Protein.Raw)
summary(check)
abline(check)

mean(subset(elisa, elisa$Temperature == "C")$conc)
mean(subset(elisa, elisa$Temperature == "HW")$conc)



#There is no effect of raw protein concentrations on detected protein carbonyl concentrations
#All analyses are run using protein carbonyl concentrations from ELISA output
#results = lme(Protein.Carbonyl~Temperature+Complexity+Temperature:Complexity, data = elisa, random=~1|Tower/Placement)
  #Results not impacted by random effects of tower, placement within the tower, individual weight, or individual length
#results = lme(Protein.Carbonyl~Temperature+Complexity+Temperature:Complexity, data = elisa, random=~1|Weight)
#hist(residuals(results))
#summary(results)
#intervals(results)

results = lme(conc~Weight+Temperature+Complexity+Temperature:Complexity, data = elisa, random=~1|Tower)
summary(results)
#Results not impacted by random effects of tower, placement within the tower, individual weight, or individual length
hist(residuals(results))
summary(results)
intervals(results)

results = lm(conc~Complexity, data = elisa)
summary(results)
hist(residuals(results))
summary(results)
intervals(results)

results = lm(conc~Temperature, data = elisa)
summary(results)
hist(residuals(results))
summary(results)
intervals(results)

TukeyHSD(conc~Weight+Temperature+Complexity+Temperature:Complexity, data = elisa, random=~1|Tower)


par(mfrow = c(1,2))
plot(residuals(results)~fitted(results), main = "Residuals vs Fitted")
qqnorm(residuals(results))



#Plot results

bp = ggplot(elisa, aes(y = conc, x = Temperature, fill = Temperature)) +
  geom_point(aes(y = conc, color = Temperature), position = position_jitter(width = .25), size = 3, alpha = 0.6) +
  geom_boxplot(width = .5, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 2) +
  guides(fill = "none") +
  guides(color = "none") +
  scale_fill_manual(values=c("#009E73", "#D55E00")) +
  scale_colour_manual(values=c("#009E73", "#D55E00")) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 15, colour = "black")) +
  ylab("Protein Carbonyl Concentration (nmol/mg)") +
  xlab("Treatment") +
  annotate("text", x=2.1, y=0.1, label= "p = 0.007", size = 8)
bp

pdf("ELISA_liver_Temp_boxplot.pdf")
bp
dev.off()
