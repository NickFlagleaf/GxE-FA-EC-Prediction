################################ Summary function for fitted Factor Analytic mixed models in AS-REML#################################################
################### This work © 2025 by Nick Fradgley, CSIRO is licensed under CC BY 4.0################################
## Some code has been adapted from Supplementary Material in Smith et al. 2021 (https://doi.org/10.3389/fpls.2021.737462)

fa.sum <- function(fa.obj, addGfac = "GID", nonaddGfac = "GID", Efac = "Env", add = T, non.add = T) {
  fa.sum <- summary(fa.obj)
  vc <- fa.sum$varcomp
  # get FA loadings
  gss <- vc[grepl(addGfac, rownames(vc)), ]
  # for additive G effects
  if (add == T) {
    addpsi <- gss[grepl(pattern = paste(Efac, ":vm\\(", addGfac, sep = ""), x = rownames(gss)), "component"] # specific additive G vars
    addLL <- gss[grepl(pattern = "!fa", x = rownames(gss)) & grepl(pattern = paste("vm\\(", addGfac, sep = ""), x = rownames(gss)), "component"]
    addL <- matrix(addLL, nrow = length(addpsi), ncol = length(addLL) / length(addpsi), byrow = F)
    names(addpsi) <- unique(fa.obj$mf$Env)
    rownames(addL) <- unique(fa.obj$mf$Env)
  } else {
    addL <- NA
  }
  # for nonadditive G effects
  if (non.add == T) {
    nonaddpsi <- gss[grepl(pattern = paste(Efac, ":", nonaddGfac, "!", sep = ""), x = rownames(gss)), "component"] # specific nonadditive G vars
    nonaddLL <- gss[grepl(pattern = paste(nonaddGfac, "!.*!fa", sep = ""), x = rownames(gss)), "component"]
    nonaddL <- matrix(nonaddLL, nrow = length(nonaddpsi), ncol = length(nonaddLL) / length(nonaddpsi), byrow = F)
    names(nonaddpsi) <- unique(fa.obj$mf$Env)
    rownames(nonaddL) <- unique(fa.obj$mf$Env)
  } else {
    nonaddL <- NA
  }

  process.FA.loads <- function(psi, L) {
    svdL <- svd(L) # rotate loadings
    u1neg <- 100 * length(svdL$u[, 1][svdL$u[, 1] < 0]) / dim(L)[[1]]
    if (u1neg > 50) {
      svdL$u <- -1 * svdL$u
      svdL$v <- -1 * svdL$v
    } # make svdL$d a diagonal matrix here - especially important
    # for nfac=1 case - also save squared singular values
    svdL$d2 <- diag(svdL$d^2, nrow = length(svdL$d))
    svdL$d <- diag(svdL$d, nrow = length(svdL$d))
    Lam <- svdL$u
    Dmat <- svdL$d2

    iClasses <- apply(Lam, 1, function(x) paste(c("n", "p")[as.numeric(x > 0) + 1], collapse = ""))
    names(iClasses) <- levels(as.data.frame(fa.obj$mf)[, Efac])

    dimnames(Lam) <- list(unique(fa.obj$mf$Env), paste("FA", 1:ncol(Lam), sep = ""))
    dimnames(Dmat) <- list(paste("FA", 1:ncol(Lam), sep = ""), paste("FA", 1:ncol(Lam), sep = ""))
    gcov <- Lam %*% Dmat %*% t(Lam) + diag(psi)
    gcor <- cov2cor(gcov)

    paf.site <- matrix(0, nrow = nrow(Lam), ncol = ncol(Lam), dimnames = dimnames(Lam))
    for (i in 1:ncol(Lam)) {
      paf.site[, i] <- 100 * (Dmat[i, i] * diag(Lam[, i] %*% t(Lam[, i]))) / diag(gcov)
    }
    if (ncol(Lam) > 1) {
      all <- 100 * diag(Lam %*% Dmat %*% t(Lam)) / diag(gcov)
      paf.site <- cbind(paf.site, all)
    }
    paf.mod <- 100 * sum(diag(Lam %*% Dmat %*% t(Lam))) / sum(diag(gcov))
    paf.fac <- 100 * diag(Dmat) / sum(diag(gcov))
    lamout <- list(
      "Lam" = Lam, "iClasses" = iClasses, "Dmat" = Dmat, "svdL" = svdL, "Gcor" = gcor, "Psi" = psi,
      "Gcov" = gcov, "varexpl_components" = paf.site,
      "factor %vaf" = paf.fac, "total %vaf" = paf.mod
    )
    return(lamout)
  }
  if (add == T) {
    add.G.effs <- process.FA.loads(psi = addpsi, L = addL)
  } else {
    add.G.effs <- NA
  }
  if (non.add == T) {
    nonadd.G.effs <- process.FA.loads(psi = nonaddpsi, L = nonaddL)
  } else {
    nonadd.G.effs <- NA
  }


  # get all additive BLUPS per environment
  if (add == T) {
    Envs <- as.character(levels(fa.obj$mf$Env))
    Gids <- as.character(levels(fa.obj$mf$GID))
    nenvs <- length(Envs)
    ngeno <- length(Gids)
    add.n.fac <- ncol(add.G.effs$Lam)
    GIDsinGmat <- fa.obj$G.param[[grep(pattern = paste("rr\\(", Efac, "\\, ", add.n.fac, "\\)\\:vm", sep = ""), x = names(fa.obj$G.param))]][[3]]$levels

    CVEs <- sapply(Envs, function(x) {
      fa.obj$coef$random[grep(pattern = paste("rr\\(", Efac, "\\, ", add.n.fac, "\\)_", x, "\\:vm", sep = ""), rownames(fa.obj$coef$random)), ]
    })
    SVEs <- sapply(Envs, function(x) {
      fa.obj$coef$random[grep(pattern = paste(Efac, "_", x, "\\:vm", sep = ""), rownames(fa.obj$coef$random)), ]
    })
    rownames(CVEs) <- GIDsinGmat
    rownames(SVEs) <- GIDsinGmat
    add.EBLUPS <- CVEs + SVEs
    rownames(add.EBLUPS) <- GIDsinGmat
  } else {
    add.EBLUPS <- NA
  }


  # Get FA loading additive genotype slope scores
  coef.all <- fa.obj$coef$random
  scores <- coef.all[grepl(paste("Comp.*:vm\\(", addGfac, sep = ""), x = rownames(coef.all))] # Get the additive FA genotype scores
  score.mat <- matrix(scores, ncol = add.n.fac)
  score.mat <- score.mat %*% add.G.effs$svdL$v %*% add.G.effs$svdL$d
  dimnames(score.mat) <- list(GIDsinGmat, paste("Comp", rep(1:add.n.fac)))

  iclass.list <- sapply(levels(factor(add.G.effs$iClasses)), function(x) names(add.G.effs$iClasses[add.G.effs$iClasses == x]))
  Iclass.mean.Lams <- lapply(iclass.list, function(x) sapply(1:ncol(add.G.effs$Lam), function(f) mean(add.G.effs$Lam[x, f]))) # get mean Lams per iclass
  iClass.OP.scores <- sapply(Iclass.mean.Lams, function(x) rowSums(score.mat %*% x)) # add up genotype scores at each mean iclass mean lams

  # Calculate Overall Performance and RMSD stability
  OP <- mean(add.G.effs$Lam[, 1]) * score.mat[, 1]
  Eij <- t(apply(CVEs, 1, function(x) resid(lm(x ~ add.G.effs$Lam[, 1]))))
  RMSD <- apply(Eij, 1, function(x) mean(abs(x)))

  iClass.info <- list(
    "iClass.list" = iclass.list,
    "iClassOP" = iClass.OP.scores
  )
  GID.prfmc <- list("OP" = OP, "RMSD" = RMSD)

  # Get environment fixed main effects
  E.main.effs <- fa.obj$coef$fixed
  int <- E.main.effs["(Intercept)", ]
  Ms <- E.main.effs[!rownames(E.main.effs) == "(Intercept)", ]
  E.main.effs <- Ms + int
  names(E.main.effs) <- gsub("Env_", "", names(E.main.effs))

  out <- list(
    "Summary" = fa.sum, "Varcoms" = vc,
    "Fixed.effects" = E.main.effs,
    "Additive_G_effects" = add.G.effs,
    "Non-additive_G_effects" = nonadd.G.effs,
    "Additive_EGBLUPs" = list(
      "CVEs" = CVEs,
      "SVEs" = SVEs,
      "BLUPs" = add.EBLUPS
    ),
    "Additive_genotype_scores" = score.mat,
    "iClass.info" = iClass.info,
    "GID.performance" = GID.prfmc
  )

  return(out)
}
