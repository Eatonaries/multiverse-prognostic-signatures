#!/usr/bin/env Rscript
# bio4 STAGE B: reuse one stage-A unit -> score every
# (depool x fs x ngenes x cutoff) spec on the held-out test part (internal
# C-index/HR), and, for "full" units, on each external queue.
#
# Stage-A layout (IMPORTANT):
#   A$train  : genes x n_train samples  (columns correspond to A$train_idx)
#   A$test   : genes x n_test  samples  (columns correspond to A$test_idx; NULL for full)
#   A$rankings[[pool]][[fs]] : gene-id ranking, fitted on TRAIN only
#
# Usage: Rscript --vanilla stageB_grid.R <payload.json>
suppressPackageStartupMessages({
  library(survival); library(data.table); library(jsonlite)
})

payload <- jsonlite::fromJSON(commandArgs(trailingOnly = TRUE)[1], simplifyVector = TRUE)
A <- readRDS(payload$stageA_rds)
out_csv <- payload$out_csv
ROOT <- Sys.getenv("BIO4_ROOT", "D:/deepseek/jieru/_bio4")
slim <- readRDS(file.path(ROOT, "data/processed",
                          paste0(tolower(payload$project), "_", payload$tag), "slim.rds"))
y <- slim$y
gene_names <- A$gene_names
tr <- A$train
tse <- A$test
train_idx <- A$train_idx
test_idx <- A$test_idx
ytr <- y[train_idx, , drop = FALSE]
yte <- if (length(test_idx)) y[test_idx, , drop = FALSE] else NULL

NGENES <- c(5, 10, 20, 50)
POOL_LEVELS <- c("limma", "edger", "deseq2")
FS_LEVELS <- c("lasso", "ridge", "enet", "rsf", "unicox")
CUTOFFS <- c("median", "bestlogrank", "tertile")

# ---- external queues (symbol-keyed) ----
geo_mats <- list()
for (gname in names(payload$geo_dirs)) {
  d <- payload$geo_dirs[[gname]]
  expr <- fread(file.path(d, "expr_genes_by_symbol.csv"))
  sym <- expr[[1]]
  m <- as.matrix(expr[, -1, with = FALSE]); rownames(m) <- sym
  colnames(m) <- gsub('"', '', colnames(m), fixed = TRUE)
  sv <- fread(file.path(d, "survival.csv"))
  geo_mats[[gname]] <- list(mat = m, surv = sv[sv$sample_id %in% colnames(m), ])
}

# ---- scoring helpers ----
zfit <- function(mat) {                       # genes x samples -> z, mu, sd
  mu <- rowMeans(mat); sd <- apply(mat, 1, stats::sd)
  ok <- is.finite(sd) & sd > 0
  list(Z = (mat[ok, , drop = FALSE] - mu[ok]) / sd[ok], mu = mu[ok], sd = sd[ok])
}
fit_cox <- function(Xtr, ytr, Xte = NULL) {
  d <- data.frame(Xtr); d$.t <- ytr$time; d$.e <- ytr$event
  f <- tryCatch(suppressWarnings(coxph(Surv(.t, .e) ~ ., data = d, x = TRUE)),
                error = function(e) NULL)
  if (is.null(f) || length(coef(f)) == 0 || any(!is.finite(coef(f)))) return(NULL)
  lp <- if (!is.null(Xte))
    tryCatch(as.numeric(predict(f, newdata = data.frame(Xte), type = "lp")),
             error = function(e) NULL) else NULL
  list(model = f, lp = lp)
}
cidx <- function(lp, time, event) {
  if (is.null(lp) || length(lp) < 5 || sd(lp) == 0 || !all(is.finite(lp))) return(NA_real_)
  concordance(Surv(time, event) ~ I(-lp))$concordance
}
hr_p <- function(lp, time, event, cutoff) {
  g <- tryCatch({
    if (cutoff == "median") {
      factor(lp > median(lp))
    } else if (cutoff == "tertile") {
      b <- stats::quantile(lp, c(1/3, 2/3), names = FALSE)
      factor(cut(lp, breaks = c(-Inf, b, Inf)))
    } else {  # bestlogrank
      u <- sort(unique(lp))
      cand <- u[seq_along(u) > 5 & seq_along(u) < length(u) - 5]
      if (length(cand) == 0) return(NULL)
      cand <- cand[seq(1, length(cand), length.out = min(80, length(cand)))]
      bb <- Inf; bg <- NULL
      for (thr in cand) {
        gg <- factor(lp > thr)
        if (length(unique(gg)) < 2 || any(is.na(gg))) next
        f <- tryCatch(coxph(Surv(time, event) ~ gg), error = function(e) NULL)
        if (is.null(f)) next
        s <- summary(f)$coefficients
        p <- if (nrow(s) >= 1) as.numeric(s[1, "Pr(>|z|)"]) else 1
        if (is.finite(p) && p < bb) { bb <- p; bg <- gg }
      }
      bg
    }
  }, error = function(e) NULL)
  if (is.null(g) || length(g) == 0 || length(unique(g)) < 2) return(c(NA_real_, NA_real_))
  f <- tryCatch(coxph(Surv(time, event) ~ g), error = function(e) NULL)
  if (is.null(f)) return(c(NA_real_, NA_real_))
  s <- summary(f)$coefficients
  if (nrow(s) == 0) return(c(NA_real_, NA_real_))
  c(exp(s[1, "coef"]), as.numeric(s[1, "Pr(>|z|)"]))
}

# ---- main loop ----
rows <- list(); k <- 0L
add <- function(rec) { k <<- k + 1L; rows[[k]] <<- rec }
for (pool in POOL_LEVELS) {
  for (fs in FS_LEVELS) {
    rk <- A$rankings[[pool]][[fs]]
    if (is.null(rk) || length(rk) == 0) next
    genes_ok <- rk[rk %in% rownames(tr)]
    if (length(genes_ok) < 2) next
    for (ng in NGENES) {
      genes <- head(genes_ok, ng)
      if (length(genes) < min(ng, 2)) next
      zf <- zfit(tr[genes, , drop = FALSE])
      Xtr <- t(zf$Z)
      c_in <- NA_real_; lp_te <- NULL
      if (!is.null(tse)) {
        common <- intersect(rownames(zf$Z), rownames(tse))
        if (length(common) >= 2) {
          zt <- sweep(sweep(tse[common, , drop = FALSE], 1, zf$mu[common], "-"),
                      1, zf$sd[common], "/")
          Xte <- t(zt)
          fit <- fit_cox(Xtr[, , drop = FALSE], ytr, Xte)
          if (!is.null(fit)) {
            lp_te <- fit$lp
            c_in <- cidx(lp_te, yte$time, yte$event)
          }
        }
      }
      for (cutoff in CUTOFFS) {
        hp <- c(NA_real_, NA_real_)
        if (!is.null(lp_te) && !is.na(c_in)) hp <- hr_p(lp_te, yte$time, yte$event, cutoff)
        add(list(project = payload$project, norm = A$norm, depool = pool, fs = fs,
                 ngenes = ng, cutoff = cutoff, split = A$split_kind, fold = A$fold,
                 seed = A$seed, extval = "internal", n_train = length(train_idx),
                 n_test = length(test_idx), cindex = c_in, hr = hp[1],
                 p_logrank = hp[2]))
      }
      # external: only from full units, once per spec (median cutoff row)
      if (A$split_kind == "full" && length(geo_mats) > 0) {
        for (gname in names(geo_mats)) {
          gm <- geo_mats[[gname]]
          g2 <- genes[genes %in% names(gene_names)]
          s2 <- as.character(gene_names[g2])
          names(s2) <- g2
          syms <- unique(intersect(as.character(s2), rownames(gm$mat)))
          if (length(syms) < 2) next
          # fold TCGA signature genes to symbol level using ALL samples (full unit)
          tr_sym <- t(sapply(syms, function(sm) {
            id <- g2[s2 == sm]
            if (length(id) == 1) tr[id, , drop = FALSE] else
              colMeans(tr[id, , drop = FALSE])
          }))
          zz <- zfit(tr_sym)
          # predict() matches newdata columns BY NAME: keep gene names explicit
          gz <- rownames(zz$Z)
          Xtr_s <- t(zz$Z)
          colnames(Xtr_s) <- gz
          fit <- fit_cox(Xtr_s, ytr)
          if (is.null(fit)) next
          Xg <- sweep(sweep(gm$mat[gz, , drop = FALSE], 1, zz$mu, "-"),
                      1, zz$sd, "/")
          newd <- as.data.frame(t(Xg))
          colnames(newd) <- gz
          lp_g <- tryCatch(as.numeric(predict(fit$model, newdata = newd,
                                              type = "lp")),
                           error = function(e) NULL)
          if (is.null(lp_g)) next
          # align risk score to the survival table BY sample id
          lp_df <- data.frame(sample_id = colnames(Xg), lp = lp_g,
                              stringsAsFactors = FALSE)
          mg <- merge(gm$surv, lp_df, by = "sample_id")
          if (nrow(mg) < 10) next
          cg <- cidx(mg$lp, mg$time, mg$event)
          cg_css <- if ("event_css" %in% names(mg))
            cidx(mg$lp, mg$time, mg$event_css) else NA_real_
          if (!is.na(cg))
            add(list(project = payload$project, norm = A$norm, depool = pool, fs = fs,
                     ngenes = ng, cutoff = "median", split = "full", fold = 0L,
                     seed = A$seed, extval = gname, n_train = nrow(y),
                     n_test = nrow(gm$surv), cindex = cg, hr = NA_real_,
                     p_logrank = NA_real_, cindex_css = cg_css))
        }
      }
    }
  }
}
df <- rbindlist(lapply(rows, as.data.frame), fill = TRUE)
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
data.table::fwrite(df, out_csv)
cat(sprintf("STAGE B %s__%s__f%02d: %d rows (internal=%d external=%d)\n",
            A$split_kind, A$norm, A$fold, nrow(df),
            sum(df$extval == "internal"), sum(df$extval != "internal")))
