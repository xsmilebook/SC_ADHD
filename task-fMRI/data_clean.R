# 读取RDS文件
filepath <- "D:/code/SC_ADHD/datasets/ABCD/task-fMRI/activation/SCdata_Yeo17_CV75_Activation.sum.msmtcsd.combat_TD_ADHDall.rds"
data <- readRDS(filepath)

# 查看原始数据信息
cat("原始数据维度:", dim(data)[1], "行", dim(data)[2], "列\n")

# 计算每列的空值数量
na_count <- sapply(data, function(x) sum(is.na(x)))
print("各列空值计数:")
print(na_count)

# 仅剔除所有值都为NA的行（整行都是NA的行）
# 使用complete.cases()或手动判断
all_na_rows <- apply(is.na(data), 1, all)
data_clean <- data[!all_na_rows, ]

# 或者使用另一种方法：
# data_clean <- data[!apply(is.na(data), 1, all), ]

# 查看清理后的结果
cat("\n清理结果:\n")
cat("原始行数:", nrow(data), "\n")
cat("清理后行数:", nrow(data_clean), "\n")
cat("删除的全NA行数:", sum(all_na_rows), "\n")

# 保存到CSV文件
output_path <- "D:/code/SC_ADHD/datasets/ABCD/task-fMRI/activation/SCdata_Yeo17_CV75_Activation.sum.msmtcsd.combat_TD_ADHDall_clean.csv"
write.csv(data_clean, output_path, row.names = FALSE)

cat("\n数据已保存到:", output_path, "\n")