setwd("~/renyimeng/GAGNAR")
library(dplyr)
library(ggplot2)
city_int <- read.csv("city_interval.csv")
city_int <- as.matrix(city_int)
city_long <- rbind(city_int[,1:3],city_int[,4:6],
                   city_int[,7:9],city_int[,10:12]) %>% as.data.frame()
colnames(city_long) <- c("est", "lower", "upper")
city_long$group <- c(rep(1, 7), rep(2, 7), rep(3, 7), rep(4, 7)) %>% factor(levels = c(4,3,2,1))
city_long$para <- rep(c("beta0","beta1","beta2","gamma1","gamma2","gamma3","gamma4"), 4) %>%
  factor(levels = rev(c("beta0","beta1","beta2","gamma1","gamma2","gamma3","gamma4")))

mycolor = c("#fac943", "#203653","#7ca729","#c65b66")
pdf(width=5, height=4, file = "./report/city_interval.pdf")
ggplot(city_long, aes(x = para, group = group)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, color = group), 
                position = position_dodge(0.6), size = 0.8) +
  scale_colour_manual(values = mycolor, limits = c(1,2,3,4)) +
  geom_point(aes(y = est, group = group), color = "red", size = 1, position = position_dodge(0.6)) +
  theme_bw() +
  theme(legend.position = c(0.93, 0.18),
        legend.background=element_blank(),
        legend.title = element_text(face="italic"),
        axis.text = element_text(size=13)) +
  labs(x = "", y = "Estimates") + ylim(-1,1) +
  coord_flip() +
  scale_x_discrete(labels = rev(expression(beta[0], beta[net], beta[mom],
                                         gamma[pop], gamma[gdp1], gamma[sav], gamma[fore])))
dev.off()

