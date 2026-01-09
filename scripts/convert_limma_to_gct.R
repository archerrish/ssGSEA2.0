#!/usr/bin/env Rscript
#
# 将 limma 差异分析结果（TSV）转换为 PTM-SEA 所需的 GCT 格式
# 
# 重要：PTM-SEA/GSEA 需要**全部位点**的排序，不能过滤！
#
# 输入：LUAD_phospho_paired_limma.tsv
#   - 列: logFC, t, P.Value, adj.P.Val 等
#   - 位点ID格式: ENSG00000003056.8|ENSP00000000412.3|S267|DDQLGEESEERDDHL|1
#
# 输出：GCT 文件，行为位点，列为 signed_logP（单样本）
#

# 设置 CRAN 镜像
options(repos = c(CRAN = "https://cloud.r-project.org"))

# 安装/加载依赖
if(!require(pacman)) install.packages("pacman")
pacman::p_load(data.table, magrittr)

# 使用项目自带的 io.R（包含 GCT 写入功能）
source("src/io.R")

# ============================================
# 参数设置
# ============================================
input_file <- "example/LUAD_phospho_paired_limma.tsv"
output_prefix <- "example/LUAD_phospho_ptmsea"

# ============================================
# 读取数据
# ============================================
cat("读取 limma 结果...\n")
df <- fread(input_file, sep="\t", header=TRUE)

cat("原始位点数:", nrow(df), "\n")
cat("列名:", colnames(df), "\n\n")

# ============================================
# 提取 7AA 序列窗（flanking sequence）
# ============================================
# idx 格式: ENSG|ENSP|S267|DDQLGEESEERDDHL|1
#           ^--- 第4列是 7AA flanking

cat("提取 flanking sequence...\n")
df[, flanking := sapply(strsplit(idx, "\\|"), function(x) x[4])]

# 检查是否有重复的 flanking
if(any(duplicated(df$flanking))) {
  n_dup <- sum(duplicated(df$flanking))
  cat("警告: 检测到", n_dup, "个重复的 flanking sequence\n")
  cat("这些可能是同一位点在不同蛋白异构体/基因上的情况\n")
  cat("将使用 make.unique() 保留所有位点\n\n")
  df$flanking <- make.unique(df$flanking, sep="_")
}

# ============================================
# 提取修饰类型（S/T/Y 等）
# ============================================
df[, mod_type := sapply(strsplit(idx, "\\|"), function(x) {
  site_info <- x[3]
  substr(site_info, 1, 1)  # 第一个字符（S/T/Y）
})]

# 添加 "-p" 后缀（磷酸化）
df[, site_id := paste0(flanking, "-p")]

# ============================================
# 计算 signed -log10(p)
# ============================================
cat("计算 signed -log10(p)...\n")
df[, signed_logP := -log10(P.Value) * sign(logFC)]

# 处理 p=0 的情况（设为 -log10 的最小非零值）
if(any(df$P.Value == 0)) {
  min_nonzero_p <- min(df$P.Value[df$P.Value > 0], na.rm=TRUE)
  n_zero <- sum(df$P.Value == 0)
  cat("警告:", n_zero, "个位点 p-value = 0，将设为 p =", min_nonzero_p, "\n")
  df[P.Value == 0, signed_logP := -log10(min_nonzero_p) * sign(logFC)]
}

cat("signed_logP 范围:", round(range(df$signed_logP, na.rm=TRUE), 2), "\n\n")

# ============================================
# 准备 GCT 矩阵
# ============================================
cat("构建 GCT 对象...\n")

# 数据矩阵（单列，全部位点！）
gct_matrix <- matrix(df$signed_logP, ncol=1)
rownames(gct_matrix) <- df$site_id
colnames(gct_matrix) <- "LUAD_Tumor_vs_NAT"

# 行注释（row descriptions）
rdesc <- data.frame(
  id.original = df$idx,
  flanking = df$flanking,
  mod_type = df$mod_type,
  logFC = df$logFC,
  t = df$t,
  P.Value = df$P.Value,
  adj.P.Val = df$adj.P.Val,
  signed_logP = df$signed_logP,
  stringsAsFactors = FALSE
)
rownames(rdesc) <- df$site_id

# 列注释（column descriptions）
cdesc <- data.frame(
  sample = "LUAD_Tumor_vs_NAT",
  comparison = "Tumor vs NAT",
  method = "limma paired",
  n_pairs = df$n_pairs_nonNA[1],
  stringsAsFactors = FALSE
)
rownames(cdesc) <- "LUAD_Tumor_vs_NAT"

# ============================================
# 创建 GCT 对象
# ============================================
gct <- new("GCT")
gct@mat <- gct_matrix
gct@rid <- df$site_id
gct@cid <- "LUAD_Tumor_vs_NAT"
gct@rdesc <- rdesc
gct@cdesc <- cdesc

# ============================================
# 导出 GCT
# ============================================
output_file <- paste0(output_prefix, ".gct")
cat("写入 GCT 文件:", output_file, "\n")
write.gct(gct, ofile=output_file, appenddim=FALSE)

# ============================================
# 生成统计摘要
# ============================================
cat("\n========== 转换完成 ==========\n")
cat("输出文件:", output_file, "\n")
cat("总位点数:", nrow(gct@mat), "（全部位点，未过滤）\n")
cat("样本数:", ncol(gct@mat), "\n")
cat("\nsigned_logP 统计:\n")
print(summary(gct@mat[,1]))

cat("\n位点分布:\n")
cat("  显著上调 (FDR<0.05 & logFC>0):", sum(df$adj.P.Val < 0.05 & df$logFC > 0, na.rm=TRUE), "\n")
cat("  显著下调 (FDR<0.05 & logFC<0):", sum(df$adj.P.Val < 0.05 & df$logFC < 0, na.rm=TRUE), "\n")
cat("  不显著:", sum(df$adj.P.Val >= 0.05, na.rm=TRUE), "\n")

cat("\nsigned_logP 阈值分布:\n")
cat("  强上调信号 (>10):", sum(gct@mat[,1] > 10, na.rm=TRUE), "\n")
cat("  中等上调 (2~10):", sum(gct@mat[,1] > 2 & gct@mat[,1] <= 10, na.rm=TRUE), "\n")
cat("  弱信号 (-2~2):", sum(abs(gct@mat[,1]) <= 2, na.rm=TRUE), "\n")
cat("  中等下调 (-10~-2):", sum(gct@mat[,1] < -2 & gct@mat[,1] >= -10, na.rm=TRUE), "\n")
cat("  强下调信号 (<-10):", sum(gct@mat[,1] < -10, na.rm=TRUE), "\n")

cat("\n可用于 PTM-SEA 的命令:\n")
cat("Rscript ssgsea-cli.R \\\n")
cat("  -i", output_file, "\\\n")
cat("  -d db/ptmsigdb/v2.0.0/all/ptm.sig.db.all.flanking.human.v2.0.0.gmt \\\n")
cat("  -o outputs/LUAD_ptmsea \\\n")
cat("  -w 0 \\\n")
cat("  -c z.score \\\n")
cat("  -m 5 \\\n")
cat("  -p 1000 \\\n")
cat("  -z .\n")
