#!/usr/bin/env Rscript
#
# 将 PTM-SEA 结果（GCT 格式）转换为易读的汇总表格
#
# 输入: outputs/LUAD_ptmsea-*.gct
# 输出: outputs/LUAD_ptmsea_summary.tsv (包含 Signature, NES, p-value, FDR 等)
#

if(!require(pacman)) install.packages("pacman")
pacman::p_load(data.table, magrittr)

source("src/io.R")

# ============================================
# 参数设置
# ============================================
output_prefix <- "outputs/LUAD_ptmsea"
sample_name <- "LUAD_Tumor_vs_NAT"

# ============================================
# 读取 GCT 文件
# ============================================
cat("读取 PTM-SEA 结果文件...\n")

# NES scores
scores <- parse.gctx(paste0(output_prefix, "-scores.gct"))
nes_values <- scores@mat[, sample_name]
signatures <- scores@rid

# p-values
pval_gct <- parse.gctx(paste0(output_prefix, "-pvalues.gct"))
pval_values <- pval_gct@mat[, sample_name]

# FDR adjusted p-values
fdr_gct <- parse.gctx(paste0(output_prefix, "-fdr-pvalues.gct"))
fdr_values <- fdr_gct@mat[, sample_name]

# 行注释（包含 signature 大小、重叠等信息）
rdesc <- scores@rdesc

# ============================================
# 创建汇总表格
# ============================================
cat("构建汇总表格...\n")

summary_df <- data.frame(
  Signature = signatures,
  NES = round(nes_values, 4),
  p_value = format(pval_values, scientific=TRUE, digits=3),
  FDR = format(fdr_values, scientific=TRUE, digits=3),
  Signature.size = rdesc$Signature.set.size,
  Overlap.n = rdesc[[paste0("Signature.set.overlap.", sample_name)]],
  Overlap.percent = rdesc[[paste0("Signature.set.overlap.percent.", sample_name)]],
  stringsAsFactors = FALSE
)

# 按 |NES| 排序（最显著的在前面）
summary_df <- summary_df[order(abs(summary_df$NES), decreasing=TRUE), ]
rownames(summary_df) <- NULL

# ============================================
# 分类汇总
# ============================================
cat("\n===== 汇总统计 =====\n")
cat("总 Signatures 数:", nrow(summary_df), "\n")

# FDR<0.05 的显著结果
sig_idx <- as.numeric(fdr_values) < 0.05
cat("显著 (FDR<0.05):", sum(sig_idx), "\n")
cat("  - 富集 (NES>0):", sum(sig_idx & nes_values > 0), "\n")
cat("  - 耗尽 (NES<0):", sum(sig_idx & nes_values < 0), "\n")

# NES 范围
cat("\nNES 范围:", round(min(nes_values), 2), "~", round(max(nes_values), 2), "\n")

# ============================================
# 导出结果表格
# ============================================
output_file <- paste0(output_prefix, "_summary.tsv")
cat("\n写入汇总文件:", output_file, "\n")
write.table(summary_df, file=output_file, sep="\t", quote=FALSE, row.names=FALSE)

# ============================================
# 显示前 10 个结果（按 |NES| 排序）
# ============================================
cat("\n===== 前 10 个最显著的 Signatures =====\n")
print(head(summary_df, 10))

cat("\n===== 前 10 个负方向最显著的 Signatures (NES<0) =====\n")
neg_df <- summary_df[as.numeric(summary_df$NES) < 0, ]
print(head(neg_df, 10))

# ============================================
# 显著性分布直方图
# ============================================
cat("\n===== NES 分布 =====\n")
cat("NES < -2:", sum(as.numeric(summary_df$NES) < -2), "\n")
cat("-2 <= NES < -1:", sum(as.numeric(summary_df$NES) >= -2 & as.numeric(summary_df$NES) < -1), "\n")
cat("-1 <= NES < 0:", sum(as.numeric(summary_df$NES) >= -1 & as.numeric(summary_df$NES) < 0), "\n")
cat("0 <= NES < 1:", sum(as.numeric(summary_df$NES) >= 0 & as.numeric(summary_df$NES) < 1), "\n")
cat("1 <= NES < 2:", sum(as.numeric(summary_df$NES) >= 1 & as.numeric(summary_df$NES) < 2), "\n")
cat("NES >= 2:", sum(as.numeric(summary_df$NES) >= 2), "\n")

cat("\n===== Signature 大小分布 =====\n")
cat("大小范围:", min(rdesc$Signature.set.size), "-", max(rdesc$Signature.set.size), "\n")
cat("中位数:", median(as.numeric(rdesc$Signature.set.size)), "\n")
cat("平均数:", round(mean(as.numeric(rdesc$Signature.set.size)), 1), "\n")

# ============================================
# 按类别显示
# ============================================
cat("\n===== 按 Signature 类别统计 =====\n")
sig_category <- sub("-.*", "", signatures)
category_counts <- table(sig_category)
print(sort(category_counts, decreasing=TRUE))

cat("\n✓ 完成！查看结果:\n")
cat("  详细表格 (TSV):", output_file, "\n")
cat("  在 Excel 中打开或用 R 加载:\n")
cat("  > summary <- read.delim('", output_file, "')\n", sep="")
cat("  > head(summary)\n")
