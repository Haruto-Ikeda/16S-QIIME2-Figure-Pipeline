# ancombc2.R — #11・#12 の ANCOM-BC2 (R ANCOMBC・Dunnett 型) を計算する (261018 本人決定。BC2 実装手順書)
# 使い方: Rscript ancombc2.R <作業フォルダ> <対照群のコード> <群のコード (カンマ区切り・先頭が対照群)> <乱数の種>
# 入力: <作業フォルダ>/counts.csv (分類群 × 検体・疑似カウント前)・meta.csv (検体, grp)
#   群名は空白などを含みうるので、Python 側で g1, g2, … のコードに置き換えて渡す
# 出力: res_dunn.csv (2 群のときは主解析 res)・zero_ind.csv・warnings.txt (R の警告。Python 側が [注意] でログに出す)・session.txt
# 設定 (261016〜18 本人決定): pseudo 0 (U1)・struc_zero TRUE / neg_lb FALSE (U2)・p_adj BH (U3)・
#   mdFDR holm / B 1000 (U4)・fix "grp" / rand なし (U5)・Dunnett 型だけ (U6)・prv_cut 0 (O1)・lib_cut 0・s0_perc 0.05・alpha 0.05・n_cl 1
args <- commandArgs(TRUE)
wd <- args[1]; ctrl <- args[2]; lv <- strsplit(args[3], ",")[[1]]; seed <- as.integer(args[4])
suppressMessages({library(ANCOMBC); library(phyloseq)})
cnt <- as.matrix(read.csv(file.path(wd, "counts.csv"), row.names = 1, check.names = FALSE))
meta <- read.csv(file.path(wd, "meta.csv"), row.names = 1, check.names = FALSE)
meta$grp <- factor(meta$grp, levels = c(ctrl, setdiff(lv, ctrl)))       # 先頭 = 対照群 (Dunnett の基準)
ps <- phyloseq(otu_table(cnt, taxa_are_rows = TRUE), sample_data(meta))
warns <- character(0)
set.seed(seed)
out <- withCallingHandlers(
  ancombc2(data = ps, fix_formula = "grp", rand_formula = NULL,
           p_adj_method = "BH", pseudo = 0, pseudo_sens = TRUE,
           prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
           group = "grp", struc_zero = TRUE, neg_lb = FALSE,
           alpha = 0.05, n_cl = 1, verbose = FALSE,
           global = FALSE, pairwise = FALSE, dunnet = length(lv) >= 3, trend = FALSE,
           mdfdr_control = list(fwer_ctrl_method = "holm", B = 1000)),
  warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
# 2 群のときは Dunnett 型が使えない (ANCOMBC が切る)。比較は 1 つなので主解析の結果 (BH だけ・mdFDR なし)を同じ形で書く (261020 本人決定)
res <- if (length(lv) >= 3) out$res_dunn else out$res
write.csv(res, file.path(wd, "res_dunn.csv"), row.names = FALSE)
write.csv(out$zero_ind, file.path(wd, "zero_ind.csv"), row.names = FALSE)
writeLines(unique(warns), file.path(wd, "warnings.txt"))
writeLines(c(paste("ANCOMBC", packageVersion("ANCOMBC")), paste("seed", seed)), file.path(wd, "session.txt"))
