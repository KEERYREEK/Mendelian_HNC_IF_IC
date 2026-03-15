
library(stringr)
library(data.table)
library(readr)
library(dplyr)
library(RMediation)


#######IF2IC2HNC

record_1 <- fread("IC_IF/2HNC/IF2IC/recorded_table.csv")
record_IF <- fread("IF_HNC/IF2HNC/recorded_table.csv")
record_IC <- fread("HNC_IC/IC2HNC/recorded_table.csv")
record_data <- as.data.frame(t(as.data.frame(strsplit(record_1$name, "_TO_", fixed = TRUE))))
IF <- record_data$V1
IC <- record_data$V2
results <- data.frame()
results1 <- data.frame()
results2 <- data.frame()
results3 <- data.frame()
# i <- 1
for (i in 1:length(IF)) {
  iif <- record_IF[record_IF$name == paste0(IF[i], "_TO_HNC"), ]
  iic <- record_IC[record_IC$name == paste0(IC[i], "_TO_HNC"), ]
  
  if (any(iif$threshold == as.character(record_1[i, "threshold"]))) {
    if (any(iic$threshold == as.character(record_1[i, "threshold"]))) {
      result <- c(IF[i],IC[i],record_1[i, c("b", "b_se")],
                  iif[iif$threshold == as.character(record_1[i, "threshold"]), c("b", "b_se")],
                  iic[iic$threshold == as.character(record_1[i, "threshold"]), c("b", "b_se")],
                  record_1[i, "threshold"])
      colnames(results) <- colnames(result)
      results <- rbind(results, result)
      results1 <- rbind(results1,record_1[i, ])
      results2 <- rbind(results2,iif[iif$threshold == as.character(record_1[i, "threshold"]), ])
      results3 <- rbind(results3,iic[iic$threshold == as.character(record_1[i, "threshold"]), ])
    }
  }
}
colnames(results) <- c("IF","IC","b","b_se","IF_b", "IF_b_se", "IC_b", "IC_b_se","threshold")
tem <- data.frame()
for (j in 1:nrow(results)){
  te <- medci(mu.x = results$b[j], mu.y = results$IC_b[j], se.x = results$b_se[j], se.y = results$IC_b_se[j], rho = 0, alpha = .05, type = "all", plot = F, plotCI = F)
  tem <- rbind(tem, c(te$`Monte Carlo`$`95% CI`[1], te$`Monte Carlo`$`95% CI`[2],te$`Monte Carlo`$Estimate,te$`Monte Carlo`$SE))
}
colnames(tem) <- c("Med_lower","Med_upper","Med","Med_SE")
results <- cbind(results,tem)
dir.create("IF2IC2HNC_results")
fwrite(results, "IF2IC2HNC_results/IF2IC2HNC_results.csv")
fwrite(results1, "IF2IC2HNC_results/IF2IC_results.csv")
fwrite(results2, "IF2IC2HNC_results/IF2HNC_results.csv")
fwrite(results3, "IF2IC2HNC_results/IC2HNC_results.csv")


########HNC2IF2IC

record_1 <- fread("IC_IF/HNC2/IF2IC/recorded_table.csv")
record_IF <- fread("IF_HNC/HNC2IF/recorded_table.csv")
record_IC <- fread("HNC_IC/HNC2IC/recorded_table.csv")

record_data <- as.data.frame(t(as.data.frame(strsplit(record_1$name, "_TO_", fixed = TRUE))))
IF <- record_data$V1
IC <- record_data$V2
results <- data.frame()
results1 <- data.frame()
results2 <- data.frame()
results3 <- data.frame()
# i <- 1
for (i in 1:length(IF)) {
  iif <- record_IF[record_IF$name == paste0("HNC_TO_", IF[i]), ]
  iic <- record_IC[record_IC$name == paste0("HNC_TO_", IC[i]), ]
  
  if (any(iif$threshold == as.character(record_1[i, "threshold"]))) {
    if (any(iic$threshold == as.character(record_1[i, "threshold"]))) {
      result <- c(IF[i],IC[i],record_1[i, c("b", "b_se")],
                  iif[iif$threshold == as.character(record_1[i, "threshold"]), c("b", "b_se")],
                  iic[iic$threshold == as.character(record_1[i, "threshold"]), c("b", "b_se")],
                  record_1[i, "threshold"])
      colnames(results) <- colnames(result)
      results <- rbind(results, result)
      results1 <- rbind(results1,record_1[i, ])
      results2 <- rbind(results2,iif[iif$threshold == as.character(record_1[i, "threshold"]), ])
      results3 <- rbind(results3,iic[iic$threshold == as.character(record_1[i, "threshold"]), ])
    }
  }
}
colnames(results) <- c("IF","IC","b","b_se","IF_b", "IF_b_se", "IC_b", "IC_b_se","threshold")
tem <- data.frame()
for (j in 1:nrow(results)){
  te <- medci(mu.x = results$IF_b[j], mu.y = results$b[j], se.x = results$IF_b_se[j], se.y = results$b_se[j], rho = 0, alpha = .05, type = "all", plot = F, plotCI = F)
  tem <- rbind(tem, c(te$`Monte Carlo`$`95% CI`[1], te$`Monte Carlo`$`95% CI`[2],te$`Monte Carlo`$Estimate,te$`Monte Carlo`$SE))
}
colnames(tem) <- c("Med_lower","Med_upper","Med","Med_SE")
results <- cbind(results,tem)
dir.create("HNC2IF2IC_results")
fwrite(results, "HNC2IF2IC_results/HNC2IF2IC_results.csv")
fwrite(results1, "HNC2IF2IC_results/IF2IC_results.csv")
fwrite(results2, "HNC2IF2IC_results/HNC2IF_results.csv")
fwrite(results3, "HNC2IF2IC_results/HNC2IC_results.csv")


########IC2IF2HNC

record_1 <- fread("IC_IF/2HNC/IC2IF/recorded_table.csv")
record_IF <- fread("IF_HNC/IF2HNC/recorded_table.csv")
record_IC <- fread("HNC_IC/IC2HNC/recorded_table.csv")
record_data <- as.data.frame(t(as.data.frame(strsplit(record_1$name, "_TO_", fixed = TRUE))))
IC <- record_data$V1
IF <- record_data$V2
results <- data.frame()
results1 <- data.frame()
results2 <- data.frame()
results3 <- data.frame()
# i <- 1
for (i in 1:length(IC)) {
  iif <- record_IF[record_IF$name == paste0(IF[i], "_TO_HNC"), ]
  iic <- record_IC[record_IC$name == paste0(IC[i], "_TO_HNC"), ]
  
  if (any(iif$threshold == as.character(record_1[i, "threshold"]))) {
    if (any(iic$threshold == as.character(record_1[i, "threshold"]))) {
      result <- c(IC[i],IF[i],record_1[i, c("b", "b_se")],
                  iic[iic$threshold == as.character(record_1[i, "threshold"]), c("b", "b_se")],
                  iif[iif$threshold == as.character(record_1[i, "threshold"]), c("b", "b_se")],
                  record_1[i, "threshold"])
      colnames(results) <- colnames(result)
      results <- rbind(results, result)
      results1 <- rbind(results1,record_1[i, ])
      results2 <- rbind(results2,iic[iic$threshold == as.character(record_1[i, "threshold"]), ])
      results3 <- rbind(results3,iif[iif$threshold == as.character(record_1[i, "threshold"]), ])
    }
  }
}
colnames(results) <- c("IC","IF","b","b_se", "IC_b", "IC_b_se", "IF_b", "IF_b_se", "threshold")
tem <- data.frame()
for (j in 1:nrow(results)){
  te <- medci(mu.x = results$b[j], mu.y = results$IF_b[j], se.x = results$b_se[j], se.y = results$IF_b_se[j], rho = 0, alpha = .05, type = "all", plot = F, plotCI = F)
  tem <- rbind(tem, c(te$`Monte Carlo`$`95% CI`[1], te$`Monte Carlo`$`95% CI`[2],te$`Monte Carlo`$Estimate,te$`Monte Carlo`$SE))
}
colnames(tem) <- c("Med_lower","Med_upper","Med","Med_SE")

results <- cbind(results,tem)
dir.create("IC2IF2HNC_results")
fwrite(results, "IC2IF2HNC_results/IC2IF2HNC_results.csv")
fwrite(results1, "IC2IF2HNC_results/IC2IF_results.csv")
fwrite(results2, "IC2IF2HNC_results/IC2HNC_results.csv")
fwrite(results3, "IC2IF2HNC_results/IF2HNC_results.csv")


###HNC2IC2IF

record_1 <- fread("IC_IF/HNC2/IC2IF/recorded_table.csv")
record_IF <- fread("IF_HNC/HNC2IF/recorded_table.csv")
record_IC <- fread("HNC_IC/HNC2IC/recorded_table.csv")

record_data <- as.data.frame(t(as.data.frame(strsplit(record_1$name, "_TO_", fixed = TRUE))))
IC <- record_data$V1
IF <- record_data$V2
results <- data.frame()
results1 <- data.frame()
results2 <- data.frame()
results3 <- data.frame()
# i <- 1
for (i in 1:length(IC)) {
  iif <- record_IF[record_IF$name == paste0("HNC_TO_",IF[i]), ]
  iic <- record_IC[record_IC$name == paste0("HNC_TO_",IC[i]), ]
  
  if (any(iif$threshold == as.character(record_1[i, "threshold"]))) {
    if (any(iic$threshold == as.character(record_1[i, "threshold"]))) {
      result <- c(IC[i],IF[i],record_1[i, c("b", "b_se")],
                  iic[iic$threshold == as.character(record_1[i, "threshold"]), c("b", "b_se")],
                  iif[iif$threshold == as.character(record_1[i, "threshold"]), c("b", "b_se")],
                  record_1[i, "threshold"])
      colnames(results) <- colnames(result)
      results <- rbind(results, result)
      results1 <- rbind(results1,record_1[i, ])
      results2 <- rbind(results2,iic[iic$threshold == as.character(record_1[i, "threshold"]), ])
      results3 <- rbind(results3,iif[iif$threshold == as.character(record_1[i, "threshold"]), ])
    }
  }
}
colnames(results) <- c("IC","IF","b","b_se", "IC_b", "IC_b_se", "IF_b", "IF_b_se", "threshold")
tem <- data.frame()
for (j in 1:nrow(results)){
  te <- medci(mu.x = results$IC_b[j], mu.y = results$b[j], se.x = results$IC_b_se[j], se.y = results$b_se[j], rho = 0, alpha = .05, type = "all", plot = F, plotCI = F)
  tem <- rbind(tem, c(te$`Monte Carlo`$`95% CI`[1], te$`Monte Carlo`$`95% CI`[2],te$`Monte Carlo`$Estimate,te$`Monte Carlo`$SE))
}
colnames(tem) <- c("Med_lower","Med_upper","Med","Med_SE")

results <- cbind(results,tem)
dir.create("HNC2IC2IF_results")
fwrite(results, "HNC2IC2IF_results/HNC2IC2IF_results.csv")
fwrite(results1, "HNC2IC2IF_results/IC2IF_results.csv")
fwrite(results2, "HNC2IC2IF_results/HNC2IC_results.csv")
fwrite(results3, "HNC2IC2IF_results/HNC2IF_results.csv")
