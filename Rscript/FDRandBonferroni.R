

# 1. 设置工作目录并读取数据
setwd("circ_plot/")
data <- read.csv("HNC2IF.csv", stringsAsFactors = FALSE)

# 2. 确认你的IVW p值所在列名 
p_col_name <- "p" 

# 3. 设定总体的检验次数 (IC个数)
total_tests <- 91

# 4. 计算 Bonferroni 校正
data$Bonferroni_p <- p.adjust(data[[p_col_name]], method = "bonferroni", n = total_tests)

# 5. 计算 FDR (False Discovery Rate, 通常即为 BH 方法) 校正
data$FDR_p <- p.adjust(data[[p_col_name]], method = "fdr", n = total_tests)

# 6. 检查前几行以确认计算成功
head(data[, c(p_col_name, "Bonferroni_p", "FDR_p")])

# 7. 将结果输出并保存为新的 CSV 文件
write.csv(data, "HNC2IF_adjusted_p.csv", row.names = FALSE)
