
library(readxl)
library(dplyr)
library(openxlsx)
library(stringr)
library(readr)


base_path <- "D:/code/SC_ADHD/datasets"


# --- 1. EFNY data---

behavior_path <- file.path(base_path, "EFNY/behavior")
demo_path_efny <- file.path(base_path, "EFNY/demography")


demo_df <- read_csv(file.path(behavior_path, "MRItable20250418.csv"))
demo_df <- demo_df %>%
  mutate(
    sub_id = str_remove_all(数据ID, "_"),   
    sub_id = str_replace(sub_id, "^", "sub-")  
  )

demo_df_selected <- demo_df %>%
  select(sub_id, 年龄, 性别)


ADHD_df <- read_xlsx(file.path(demo_path_efny, "QC_cibir.xlsx"))
SC_df_efny <- read_csv(file.path(demo_path_efny, "EFNY_SC.csv"))
ICV_BP <- read.csv(paste0(behavior_path, "/ICV_BP.txt"), header = F)

ADHD_df <- ADHD_df %>%
  mutate(ADHD = ifelse(str_detect(Group, "ADHD"), 1, 0))
ADHD_df$ADHD[is.na(ADHD_df$ADHD)] <- 0
ADHD_df$sub_id <- ADHD_df$ID


efny_final_df <- merge(SC_df_efny, demo_df_selected, by = "sub_id")
efny_final_df <- efny_final_df %>% 
  left_join(ADHD_df %>% select(sub_id, ADHD), by = "sub_id") %>%
  select(-`sc_path`) 


efny_final_df <- efny_final_df %>% 
  dplyr::rename(ID = `sub_id`, age = `年龄`, sex = `性别`) %>%
  mutate(age = as.numeric(age)) 



efny_final_df$sex <- ifelse(efny_final_df$sex == "男", "M", "F")

# manually modify
efny_final_df <- efny_final_df %>%
  mutate(ADHD = ifelse(ID == "sub-THU20240526276LZH", 1, ADHD))


# --- 2. PKU6 data ---

demo_path_pku6 <- file.path(base_path, "PKU6/demography")
file_path_pku6 <- file.path(demo_path_pku6, "20250421_ADHD_NC_demo_adhdrs_cantab_forCIBR.xlsx")
ICV_PKU6 <- read.csv(paste0(demo_path_pku6, "/ICV_PKU6.txt"), header = F)

SC_df_pku6 <- read_csv(file.path(demo_path_pku6, "PKU6_SC.csv"))
basic_demo.ADHD <- read_excel(file_path_pku6, sheet = "ADHD_demo")
basic_demo.NC <- read_excel(file_path_pku6, sheet = "NC_demo")


combined_df <- bind_rows(
  basic_demo.ADHD %>% mutate(ADHD = 1),
  basic_demo.NC %>% mutate(ADHD = 0)
)


combined_df$age = as.numeric(combined_df$`age(m)`) / 12
combined_df <- combined_df %>% select(-`age(m)`, -`IA-3`, -`HI-3`, -`TO-3`)
combined_df <- combined_df %>% dplyr::rename(ID = id)

# delete sub without age,sex
combined_df <- combined_df %>% 
  filter(!is.na(sex) & !is.na(age))
combined_df$ID <- paste0("sub-", combined_df$ID)

# merge data
pku6_final_df <- merge(combined_df, SC_df_pku6, by.x = "ID", by.y = "sub_id")
pku6_final_df <- pku6_final_df %>% select(-`sc_path`)

# sex
pku6_final_df$sex <- ifelse(pku6_final_df$sex == 1, "M", "F")


# ---3. CCNP data ---
demopath_CCNP = file.path(base_path, "CCNP/demography")
basic_demo.CCNP = read.csv(file.path(demopath_CCNP, "basic_demo_devCCNPPEK.csv"))


basic_demo.CCNP <- basic_demo.CCNP %>%
  select(c("scanID", "Session", "Sex", "ScanAge", "subID", "mean_fd", "ICV"))


basic_demo.CCNP <- basic_demo.CCNP %>%
  filter(!(scanID %in% c("sub-CCNPPEK0058_ses-01", "sub-CCNPPEK0059_ses-01")))


# add ICV

# ICV_BP$V3 <- NULL
names(ICV_BP) <- c("ID", "ICV")
names(ICV_PKU6) <- c("ID", "ICV")
ICV_BP <- ICV_BP %>%
  select(where(~!all(is.na(.))))
ICV_PKU6 <- ICV_PKU6 %>%
  select(where(~!all(is.na(.))))
ICV_PKU6 <- ICV_PKU6 %>%
  mutate(ID = paste0("sub-", ID))

efny_final_df <- efny_final_df %>%
  left_join(ICV_BP, by = "ID")

pku6_final_df <- pku6_final_df %>%
  left_join(ICV_PKU6, by = "ID")

cat("EFNY ICV NA number:", sum(is.na(efny_final_df$ICV)), "\n")
cat("PKU6 ICV NA number:", sum(is.na(pku6_final_df$ICV)), "\n")

# --- 3. merge datasets ---
efny_final_df <- efny_final_df %>%
  dplyr::rename(
    subID = ID,
    mean_fd = FD
  )
efny_final_df <- efny_final_df %>% 
  mutate(
    site = "EFNY",
    scanID = subID,
    Session = as.numeric(1),
  )

pku6_final_df <- pku6_final_df %>%
  rename(
    subID = ID,
    mean_fd = FD
  )
pku6_final_df <- pku6_final_df %>% 
  mutate(
    site = "PKU6",
    scanID = subID,
    Session = as.numeric(1)
  )


CCNP_final_df <- basic_demo.CCNP %>%
  rename(
    age = ScanAge,
    sex = Sex
  )
CCNP_final_df <- CCNP_final_df %>% 
  mutate(
    site = "CCNP",
    ADHD = 0
  )
CCNP_final_df$sex <- factor(CCNP_final_df$sex, levels = c("Male", "Female"), labels = c("M", "F"))


all_sites_df <- bind_rows(efny_final_df, pku6_final_df, CCNP_final_df)

all_sites_df <- all_sites_df %>%
  filter(!is.na(ICV))
rows_after_icv <- nrow(all_sites_df)


all_sites_df <- all_sites_df %>%
  filter(age <= 16.5 & age >= 6.5)
rows_after_age <- nrow(all_sites_df)
cat("删除", rows_after_icv - rows_after_age, "行 (因age > 17). 剩余行数:", rows_after_age, "\n")

mean_fd <- mean(all_sites_df$mean_fd, na.rm = TRUE)
std_fd <- sd(all_sites_df$mean_fd, na.rm = TRUE)


lower_bound <- mean_fd - 3 * std_fd
upper_bound <- mean_fd + 3 * std_fd

cat("FD 统计值: Mean =", round(mean_fd, 4), ", SD =", round(std_fd, 4), "\n")
cat("FD 离群值过滤范围 (保留区间): [", round(lower_bound, 4), ",", round(upper_bound, 4), "]\n")

# 使用 between() 函数进行过滤
all_sites_df <- all_sites_df %>%
  filter(between(mean_fd, lower_bound, upper_bound))

rows_after_fd_filter <- nrow(all_sites_df)
cat("删除", rows_after_age - rows_after_fd_filter, "行 (因FD为离群值). 最终剩余行数:", rows_after_fd_filter, "\n")


output_path <- file.path(base_path, "demography", "THREE_SITES_basic_demo.csv")
write.csv(all_sites_df, output_path, row.names = FALSE, na = "") # 使用na=""将NA值写为空字符串

