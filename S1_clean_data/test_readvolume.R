library(tidyverse)

# 直接指定那个有问题的文件路径
problem_file <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/CCNP/processed/schaefer400_nodevolume/sub-CCNPPEK0245_ses-01_Volume7.txt'

# 尝试读取它
volume_data <- read_table(problem_file, col_names = FALSE)

# 检查返回的对象是什么
print(class(volume_data))
print(dim(volume_data))
str(volume_data) # 查看详细结构
print(volume_data$X2)