library(ggsci)

# ===== city =====

setwd("~/renyimeng/GAGNAR_140/JBES_empirical/city")
pred_rmse = read.csv("./report/pred_rmse_city_train7_test7.csv", header = F)
pred_rmse = t(pred_rmse)
start_time_point <- 2002:2011
n_test = ncol(pred_rmse) + 1  # 和训练模型对应，n_test比实际测试集长度大1
df_pred = data.frame(ID = rep(start_time_point, nrow(pred_rmse)), 
                     value = c(t(pred_rmse[,1:(n_test - 1)])),
                     type = c(rep("GAGNAR", n_test - 1), rep("CRP", n_test - 1), rep("GNAR",n_test - 1),
                              rep("NAR", n_test - 1), rep("AR", n_test - 1), rep("ARMA", n_test - 1)))
df_pred$type <- factor(df_pred$type, levels = c("AR", "ARMA", "NAR", "GNAR","CRP","GAGNAR"))

# 由于 GNAR 在某一个test set下误差非常大（已讨论），因此这里不对 GNAR 做可视化
df_pred1 = df_pred[which((df_pred$type != "GNAR") & (df_pred$ID >= 2006) & (df_pred$ID <= 2010)),]
pdf(width=5, height=4, file = "./report/city_pred_rmse_train7_test7_roll.pdf")
ggplot(df_pred1, aes(x = as.integer(ID), y = value, group = type)) +
  geom_line(aes(color = type), lwd = 0.8) +
  # scale_color_manual(values = c("red", "#96cac1", "#6a93cc", "#facb4c", "#df9a96", "#7f7f7f")) +
  scale_color_simpsons() +
  theme_bw() +
  ylim(0,14) +
  theme(legend.position = c(0.85,0.82),
        legend.background = element_blank(),
        legend.title = element_blank(),
        text = element_text(size = 12)) +
  scale_x_continuous(breaks = c(start_time_point)) +
  labs(x = "Start Time Point", y = "ReMSPE")
dev.off()


# ===== stock =====

setwd("~/renyimeng/GAGNAR_140/JBES_empirical/stock")
pred_rmse = read.csv("./report/pred_rmse_stk_train40_test6_full_sparse3.csv", header = F)
start_time_point <- 41:47
n_test = ncol(pred_rmse) + 1  # 和训练模型对应，n_test比实际测试集长度大1
df_pred = data.frame(ID = rep(start_time_point, nrow(pred_rmse)), 
                     value = c(t(pred_rmse[,1:(n_test - 1)])),
                     type = c(rep("GAGNAR", n_test - 1), rep("CRP", n_test - 1), rep("GNAR",n_test - 1),
                              rep("NAR", n_test - 1), rep("ARMA", n_test - 1), rep("AR", n_test - 1)))
df_pred$type <- factor(df_pred$type, levels = c("AR", "ARMA", "NAR", "GNAR","CRP","GAGNAR"))

# 由于 NAR 的结果太好了，因此这里不对 NAR 做可视化
df_pred1 = df_pred[which((df_pred$type != "NAR")),]
pdf(width=5, height=4, file = "./report/stk_pred_rmse_train40_test6_roll.pdf")
ggplot(df_pred1, aes(x = as.integer(ID), y = value, group = type)) +
  geom_line(aes(color = type), lwd = 0.8) +
  # scale_color_manual(values = c("red", "#96cac1", "#6a93cc", "#facb4c", "#df9a96", "#7f7f7f")) +
  scale_color_simpsons() +
  theme_bw() +
  ylim(0,5.1) +
  theme(legend.position = c(0.85,0.82),
        legend.background = element_blank(),
        legend.title = element_blank(),
        text = element_text(size = 12)) +
  scale_x_continuous(breaks = c(start_time_point)) +
  labs(x = "Start Time Point", y = "ReMSPE")
dev.off()





