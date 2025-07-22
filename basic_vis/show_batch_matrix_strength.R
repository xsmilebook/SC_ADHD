homepath <- "D:/code/SC_ADHD"
source(file.path(homepath, "code", "functions", "plotmatrix.R"))

interfileFolder <- paste0(homepath, '/datasets/interfileFolder')
FigureFolder <- paste0(homepath, "/Figures/Yeo17/site_avg")
Yeoresolution <- 17
if (Yeoresolution == 7){
  Yeoresolution.delLM = 6
}else if (Yeoresolution == 17){
  Yeoresolution.delLM = 15
}
edgenum <- Yeoresolution.delLM*(Yeoresolution.delLM+1) /2

# load data
SCdata <- readRDS(paste0(interfileFolder, '/SCdata_Yeo', Yeoresolution, '_CV75_sumSCinvnode.sum.msmtcsd.merge.rds'))

# 分站点提取数据
SCdata.PKU6 <- SCdata %>% filter(site == "PKU6")
SCdata.EFNY <- SCdata %>% filter(site == "EFNY")
SCdata.CCNP <- SCdata %>% filter(site == "CCNP")

# 先计算每个站点的 mean_strength 数据框
compute_mean_strength <- function(df_site) {
  df.mean.strength <- df_site %>%
    select(starts_with("SC.")) %>% 
    pivot_longer(everything(), names_to = "SCname", values_to = "value") %>% 
    group_by(SCname) %>%
    summarise(mean_strength = mean(value, na.rm = TRUE), .groups = 'drop') %>%
    # 新增的排序步骤：
    # 1. sub("SC.", "", SCname) 提取出 "SC." 后面的数字部分（结果为字符串）
    # 2. as.numeric() 将字符串数字转换为真正的数值
    # 3. arrange() 根据这个数值进行排序，确保 "SC.1", "SC.2", ... "SC.120" 的正确顺序
    dplyr::arrange(as.numeric(sub("SC.", "", SCname)))
  
  return(df.mean.strength)
}

# 分别计算三个站点的均值
df.PKU6 <- compute_mean_strength(SCdata.PKU6)
df.EFNY <- compute_mean_strength(SCdata.EFNY)
df.CCNP <- compute_mean_strength(SCdata.CCNP)

# 找出所有数据中的最大绝对值（用于统一 color scale）
max_abs <- max(
  abs(df.PKU6$mean_strength),
  abs(df.EFNY$mean_strength),
  abs(df.CCNP$mean_strength),
  na.rm = TRUE
)

# 设置 lmthr 为统一值
lmthr <- max_abs
NAcol <- "#67001F"

# Yeo 标签定义（同前）
yeo_labels <- c(
  "VisPeri",        # 1
  "SomMotA",        # 2
  "VisCent",        # 3
  "SomMotB",        # 4
  "DorsAttnA",      # 5
  "DorsAttnB",      # 6
  "SalVentAttnA",   # 7
  "ContC",          # 8
  "DefaultC",       # 9
  "TempPar",        # 10
  "ContA",          # 11
  "SalVentAttnB",   # 12
  "DefaultA",       # 13
  "ContB",          # 14
  "DefaultB"        # 15
)

# 修改绘图函数，固定使用统一的 lmthr
generate_matrix_plot <- function(df.mean.strength, site_name) {
  variable <- "mean_strength"
  
  Fig <- plotmatrix(
    dataname = deparse(substitute(df.mean.strength)),
    variable = variable,
    ds.resolution = Yeoresolution.delLM,
    Pvar = NA,
    NAcol = NAcol,
    lmthr = lmthr,
    axeslabels = yeo_labels,
    axeslabelsGap = TRUE
  )
  
  ggsave(
    filename = file.path(FigureFolder, paste0("Yeo17_", site_name, "_mean_strength_matrix.tiff")),
    plot = Fig,
    dpi = 600,
    width = 20,
    height = 18,
    units = "cm"
  )
}

# 保存图像
generate_matrix_plot(df.PKU6, "PKU6")
generate_matrix_plot(df.EFNY, "EFNY")
generate_matrix_plot(df.CCNP, "CCNP")