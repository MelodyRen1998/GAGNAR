city_interval_hdi_95 <- round(as.data.frame(city_interval_hdi_95), 3)
city_theta = round(city_theta, 3)
city_est_char = t(city_theta)
city_est_char_1 = city_est_char
for (kk in 1:4) {
  if(kk == 1) {colindex = c(2,3)}
  if(kk == 2) {colindex = c(5,6)}
  if(kk == 3) {colindex = c(8,9)}
  if(kk == 4) {colindex = c(11,12)}
  # colindex <- ifelse(kk == 1, c(2,3), ifelse(kk == 2, c(5,6), ifelse(kk == 3, c(8,9), c(11,12))))
  city_est_char_1[,kk] = paste0(as.numeric(city_est_char[,kk]), " (", 
                                as.numeric(city_interval_hdi_95[, colindex[1]]), ", ", 
                                as.numeric(city_interval_hdi_95[, colindex[2]]), ")")
}
write.csv(city_est_char_1, "~/renyimeng/GAGNAR_140/JBES_empirical/city_est_ci.csv", row.names = F)


stk_interval_hdi_95 <- round(as.data.frame(stk_interval_95*10), 3)  # stock
stk_theta = round(Theta*10, 3)
stk_est_char = t(stk_theta)
stk_est_char_1 = stk_est_char
for (kk in 1:5) {
  if(kk == 1) {colindex = c(2,3)}
  if(kk == 2) {colindex = c(5,6)}
  if(kk == 3) {colindex = c(8,9)}
  if(kk == 4) {colindex = c(11,12)}
  if(kk == 5) {colindex = c(14,15)}
  # colindex <- ifelse(kk == 1, c(2,3), ifelse(kk == 2, c(5,6), ifelse(kk == 3, c(8,9), c(11,12))))
  stk_est_char_1[,kk] = paste0(as.numeric(stk_est_char[,kk]), " (", 
                                as.numeric(stk_interval_hdi_95[, colindex[1]]), ", ", 
                                as.numeric(stk_interval_hdi_95[, colindex[2]]), ")")
}
write.csv(stk_est_char_1, "~/renyimeng/GAGNAR_140/JBES_empirical/stk_est_ci_times10.csv", row.names = F)
