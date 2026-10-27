#!/usr/bin/env Rscript
# bio4 STAGE A: (normalization x split-fold) unit -> train/test matrices,
# gene pools (per depool level), feature rankings (per fs level), all fitted on
# the TRAINING part only. Stage B reuses these for the cheap spec dimensions
# (ngenes x cutoff x extval), turning 11,520 expensive fits into ~84 units +
# cheap aggregation.
#
# Usage: Rscript --vanilla stageA_grid.R <payload.json>
# payload: {project, tag, norm, kind, fold, seed, outdir}
#   kind: "cv" (fold 1..10 = 5-fold x 2 repeats) | "s50"/"s70" (fold 1..5 repeats)
#         | "full" (train on everything)
suppressPackageStartupMessages({
  library(survival); library(glmnet); library(randomForestSRC)
  library(limma); library(edgeR); library(DESeq2); library(data.table)
  library(jsonlite)
})

payload <- jsonlite::fromJSON(commandArgs(trailingOnly = TRUE)[1], simplifyVector = TRUE)
project <- payload$project; tag <- payload$tag
kindn <- payload$norm; split_kind <- payload$kind
fold <- as.integer(payload$fold); seed <- as.integer(payload$seed)
outdir <- payload$outdir

ROOT <- Sys.getenv("BIO4_ROOT", "D:/deepseek/jieru/_bio4")
cdir <- file.path(ROOT, "data/processed",
                  paste0(tolower(project), "_", tag))
slim <- readRDS(file.path(cdir, "slim.rds"))
C <- slim$counts; T <- slim$tpm; y <- slim$y
stopifnot(identical(rownames(C), rownames(T)), identical(rownames(C), rownames(y)))
n <- nrow(C)
gene_names <- setNames(slim$gene_ann$gene_name, slim$gene_ann$gene_id)

POOL_LEVELS <- c("limma", "edger", "deseq2")
MIN_FRAC <- 0.20; GENE_POOL <- 2000L

# ---------------- folds ----------------
make_split <- function() {
  if (split_kind == "full") {
    list(train = seq_len(n), test = integer(0))
  } else if (split_kind == "cv") {
    set.seed(seed)
    f <- sample(rep(seq_len(5), length.out = n))
    set.seed(seed + 777)
    f2 <- sample(rep(seq_len(5), length.out = n))
    te <- if (fold <= 5) which(f == fold) else which(f2 == fold - 5)
    list(train = setdiff(seq_len(n), te), test = te)
  } else if (grepl("^s", split_kind)) {   # s50 / s60 / s70
    frac <- as.numeric(substr(split_kind, 2, 3)) / 100
    set.seed(seed + 101 * fold)
    te <- order(runif(n))[seq_len(round(n * (1 - frac)))]
    list(train = setdiff(seq_len(n), te), test = te)
  } else stop("bad split kind")
}
sp <- make_split()
train_idx <- sp$train; test_idx <- sp$test

# ---------------- normalization ----------------
rank_quantile <- function(v, ref) {
  p <- rank(v, ties.method = "average") / (length(v) + 1)
  approx(seq_along(ref), ref, xout = p, rule = 2)$y
}
normalize_one <- function(M) {  # samples x genes TPM -> genes x samples
  t(log2(M + 1))
}
quant_ref <- NULL
normalize <- function() {
  ct <- C[train_idx, , drop = FALSE]; tt <- T[train_idx, , drop = FALSE]
  ce <- if (length(test_idx)) C[test_idx, , drop = FALSE] else NULL
  te_ <- if (length(test_idx)) T[test_idx, , drop = FALSE] else NULL
  if (kindn == "log2tpm1") {
    list(train = normalize_one(tt), test = if (length(test_idx)) normalize_one(te_) else NULL)
  } else if (kindn == "voom") {
    nf_fun <- if (exists("normLibSizes", where = asNamespace("edgeR"), mode = "function"))
      get("normLibSizes", asNamespace("edgeR")) else get("calcNormFactors", asNamespace("edgeR"))
    d <- DGEList(counts = t(ct)); d <- nf_fun(d)
    keep <- colMeans(ct > 1) >= MIN_FRAC
    d <- d[keep, , keep.lib.sizes = FALSE]
    v <- suppressMessages(voom(d, design = model.matrix(~ y$event[train_idx]), plot = FALSE))
    te2 <- NULL
    if (length(test_idx)) {
      lib <- rowSums(ce)
      cpm <- sweep(ce, 1, lib, "/") * 1e6
      nf <- d$samples$norm.factors[match(rownames(ce), rownames(d))]
      if (anyNA(nf)) nf <- rep(1, nrow(ce))
      cpm <- sweep(cpm, 1, nf, "/")
      te2 <- t(log2(cpm[, keep, drop = FALSE] + 1))
    }
    list(train = v$E, test = te2)
  } else if (kindn == "deseq2vst") {
    suppressWarnings({
      cd <- data.frame(cond = factor(y$event[train_idx]))
      dds <- suppressMessages(DESeqDataSetFromMatrix(round(t(ct)), cd, ~ cond))
      dds <- estimateSizeFactors(dds)
      tr <- assay(vst(dds, blind = TRUE))
      te2 <- NULL
      if (length(test_idx)) {
        # test side: own size factors, blind VST (no outcome information used)
        cd2 <- data.frame(cond = factor(y$event[test_idx]))
        dds2 <- suppressMessages(DESeqDataSetFromMatrix(round(t(ce)), cd2, ~ cond))
        dds2 <- estimateSizeFactors(dds2)
        te2 <- assay(vst(dds2, blind = TRUE))
      }
      list(train = tr, test = te2)
    })
  } else if (kindn == "quantile") {
    # Classic quantile normalisation across SAMPLES (Bolstad-style):
    #   sort each gene's values across samples; average the k-th smallest
    #   value across genes; assign that averaged rank profile back to every
    #   gene (each gene keeps its own sample ORDER, only the values change).
    #   The averaged rank profile is estimated on TRAIN only; test samples
    #   are mapped onto the same profile.
    tr0 <- normalize_one(tt)                  # genes x train samples
    avg_profile <- local({
      S <- t(apply(tr0, 1, sort))             # genes x n, each row sorted
      colMeans(S, na.rm = TRUE)               # n points: mean rank profile
    })
    map_to_profile <- function(M) {
      if (is.null(M) || !nrow(M)) return(NULL)
      out <- matrix(NA_real_, nrow(M), ncol(M), dimnames = dimnames(M))
      for (g in seq_len(nrow(M))) {
        v <- M[g, ]
        if (all(is.na(v))) next
        ord <- order(v)
        vals <- v[ord]
        # interpolate the train reference profile onto this gene's own
        # quantiles, then restore the original sample order
        p <- (rank(v, ties.method = "average")) / (length(v) + 1)
        refp <- seq_along(avg_profile) / (length(avg_profile) + 1)
        newv <- approx(refp, avg_profile, xout = p, rule = 2)$y
        out[g, ] <- newv
      }
      out
    }
    tr <- map_to_profile(tr0)
    te2 <- if (length(test_idx)) map_to_profile(normalize_one(te_)) else NULL
    list(train = tr, test = te2)
  } else stop("bad norm")
}
nm <- normalize()

# ---------------- gene pools (TRAIN only) ----------------
zfit <- function(mat) {
  mu <- rowMeans(mat); sdv <- apply(mat, 1, stats::sd)
  ok <- is.finite(sdv) & sdv > 0
  list(Z = (mat[ok, , drop = FALSE] - mu[ok]) / sdv[ok], mu = mu, sd = sdv, ok = ok)
}
build_pools <- function(norm_tr) {
  # Each DE framework is given its CANONICAL input (field practice):
  #   limma  -> the normalized log-expression matrix
  #   edgeR  -> raw counts (negative-binomial GLM on counts)
  #   DESeq2 -> raw counts (negative-binomial Wald on counts)
  # This also avoids the incompatibility between quantile normalisation
  # (which equalises all gene distributions) and count-based dispersion
  # estimation.
  ct <- C[train_idx, , drop = FALSE]              # samples x genes, raw counts
  frac <- colMeans(ct > 1)
  keep_expr <- intersect(names(frac)[frac >= MIN_FRAC], rownames(norm_tr))
  stopifnot(length(keep_expr) > 500)
  ytr <- y[train_idx, , drop = FALSE]
  PC_ids <- slim$gene_ann$gene_id
  out <- list()
  out$limma <- local({
    sub <- norm_tr[keep_expr, , drop = FALSE]
    fit <- suppressWarnings(suppressMessages(
      eBayes(lmFit(sub, model.matrix(~ ytr$event)), trend = TRUE)))
    tt <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
    head(intersect(rownames(tt), PC_ids), GENE_POOL)
  })
  out$edger <- tryCatch(local({
    # edgeR/DESeq2 take genes x samples counts
    sub <- t(as.matrix(ct[, keep_expr, drop = FALSE]))
    grp <- factor(ifelse(ytr$event == 1, "e", "c"))
    design <- model.matrix(~ grp)
    d <- DGEList(counts = round(sub))
    d <- suppressMessages(estimateDisp(d, design))
    fit <- suppressMessages(glmFit(d, design))
    lrt <- glmLRT(fit, coef = 2)
    tt <- topTags(lrt, n = Inf)$table
    head(intersect(rownames(tt), PC_ids), GENE_POOL)
  }), error = function(e) character(0))
  out$deseq2 <- tryCatch(local({
    sub <- t(as.matrix(ct[, keep_expr, drop = FALSE]))
    cd <- data.frame(cond = factor(ytr$event))
    dds <- suppressMessages(DESeqDataSetFromMatrix(round(sub), cd, ~ cond))
    dds <- suppressMessages(DESeq(dds, quiet = TRUE))
    rn <- resultsNames(dds)
    rname <- rn[grep("^cond_", rn)][1]
    res <- as.data.frame(results(dds, name = rname))
    res <- res[order(res$pvalue), ]
    head(intersect(rownames(res), PC_ids), GENE_POOL)
  }), error = function(e) character(0))
  out
}
pools <- build_pools(nm$train)

# ---------------- feature rankings (TRAIN only, per pool) ----------------
fs_rankings <- function(norm_tr, pool) {
  if (length(pool) < 30) return(list())
  zf <- zfit(norm_tr[pool, , drop = FALSE])
  X <- t(zf$Z)
  ytr <- y[train_idx, , drop = FALSE]
  S <- Surv(ytr$time, ytr$event)
  r <- list()
  r$unicox <- local({
    pv <- apply(X, 2, function(v) {
      f <- try(suppressWarnings(coxph(Surv(time, event) ~ v,
              data = data.frame(v = v, time = ytr$time, event = ytr$event))), silent = TRUE)
      if (inherits(f, "try-error")) 1 else {
        cf <- summary(f)$coefficients[1, ]
        p <- suppressWarnings(as.numeric(cf["Pr(>|z|)"]))
        if (!is.finite(p)) 1 else p
      }
    })
    names(sort(pv))
  })
  fitnet <- function(alpha) tryCatch(local({
    cv <- suppressWarnings(cv.glmnet(X, S, family = "cox", alpha = alpha,
                                     nfolds = 3, standardize = TRUE))
    b1 <- as.matrix(coef(cv, s = "lambda.1se"))
    b0 <- as.matrix(coef(cv, s = "lambda.min"))
    score <- abs(b1[, 1]); z <- score == 0
    score[z] <- abs(b0[z, 1]) - 1e-12
    names(sort(score, decreasing = TRUE))
  }), error = function(e) character(0))
  r$lasso <- fitnet(1); r$ridge <- fitnet(0); r$enet <- fitnet(0.5)
  r$rsf <- tryCatch(local({
    # rfsrc needs the Surv() call inside the formula (a pre-built Surv object
    # is rejected with "formula is incorrectly specified")
    df <- data.frame(X, time = ytr$time, event = ytr$event, check.names = FALSE)
    rf <- suppressWarnings(rfsrc(Surv(time, event) ~ .,
                                 data = df,
                                 ntree = 500, mtry = ceiling(sqrt(ncol(X))),
                                 nodesize = 15, seed = seed,
                                 importance = "permute", nsplit = 1))
    if (is.null(rf$importance)) return(character(0))
    imp <- rf$importance
    # order() returns a bare integer vector; index the names by it
    names(imp)[order(imp, decreasing = TRUE)]
  }), error = function(e) { cat("  [rsf] error:", conditionMessage(e), "\n"); character(0) })
  lapply(r, function(v) if (length(v)) intersect(v, pool) else character(0))
}
rankings <- lapply(POOL_LEVELS, function(p) fs_rankings(nm$train, pools[[p]]))
names(rankings) <- POOL_LEVELS

# ---------------- save ----------------
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
obj <- list(project = project, tag = tag, norm = kindn, split_kind = split_kind,
            fold = fold, seed = seed, train_idx = train_idx, test_idx = test_idx,
            train = nm$train, test = nm$test, pools = pools, rankings = rankings,
            gene_names = gene_names)
f <- file.path(outdir, sprintf("%s__%s__f%02d.rds", split_kind, kindn, fold))
saveRDS(obj, f, compress = "gzip")
cat(sprintf("SAVED %s train=%d test=%d pools=%s\n", basename(f),
            length(train_idx), length(test_idx),
            paste(vapply(pools, length, integer(1)), collapse = "/")))
