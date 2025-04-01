## This work © 2024 by Nick Fradgley, CSIRO is licensed under CC BY 4.0
# Expand pedigree function----
# Converts Purdy format pedigree strings to a three column table
Expand.ped <- function(ped, # A data frame of pedigree data with columns:
                            # 'GID' - vector of genotype ID names character strings 
                            # 'Cross' - vector of pedigree charcter strings in Purdy format
                       max.depth = NULL # The maximum depth to expand pedigree to
) {
  library(stringr)
  library(data.table)

  cat("Expanding pedigree...\n")
  all.gens <- list()
  if (is.null(max.depth)) {
    max.depth <- 100
  }
  cat("|", sep = "")
  splitlevels <- c(paste("/", 100:3, "/", sep = ""), "//", "/") # Create hierarchy of characters to split purdy strings by
  max.splits <- sapply(ped$Cross, function(x) names(which.min(which(sapply(splitlevels, function(s) grepl(pattern = s, x = x)))))) # Work out the highest level split character for each individual
  parents <- strsplit(x = ped$Cross, split = max.splits) # Split pedigree into parents
  pedstocheck <- ped$Cross[unlist(lapply(parents, length)) > 2] # Check pedigree has been split into no more than two parents
  if (length(pedstocheck) > 0) { # Print out any dodgy pedigrees to manually fix
    cat("Check pedigrees:")
    cat(pedstocheck)
  }
  ped <- ped[unlist(lapply(parents, length)) == 2, ]
  max.splits <- sapply(ped$Cross, function(x) names(which.min(which(sapply(splitlevels, function(s) grepl(pattern = s, x = x))))))
  parents <- strsplit(x = ped$Cross, split = max.splits)
  parents.mat <- matrix(NA, nrow = length(parents), ncol = 2, dimnames = list(1:length(parents), c("P1", "P2"))) # make empty matrix of two parents
  for (i in 1:length(parents)) {
    parents.mat[i, ] <- parents[[i]] # Add each set of parents to matrix
  }

  # Work out if either parent is a back cross
  isP1BC <- grepl(pattern = c("\\*"), x = str_sub(parents.mat[, "P1"], 2, 2)) | grepl(pattern = c("\\*"), x = str_sub(parents.mat[, "P1"], -2, -2))
  isP2BC <- grepl(pattern = c("\\*"), x = str_sub(parents.mat[, "P2"], 2, 2)) | grepl(pattern = c("\\*"), x = str_sub(parents.mat[, "P2"], -2, -2))

  # Add in extra generations for BC parent F1s
  if (sum(isP1BC) > 0) {
    P1.split.pos <- sapply(parents.mat[isP1BC, "P1"], function(x) unlist(gregexpr(pattern = "\\*", x))[unlist(gregexpr(pattern = "\\*", x)) %in% c(2, str_length(x) - 1)])
    recrnt.spltP1 <- lapply(1:length(P1.split.pos), function(x) {
      c(
        substr(parents.mat[isP1BC, "P1"][x], 1, P1.split.pos[x] - 1),
        substr(parents.mat[isP1BC, "P1"][x], P1.split.pos[x] + 1, str_length(parents.mat[isP1BC, "P1"][x]))
      )
    })
    is.no.BCP1 <- lapply(recrnt.spltP1, function(x) sapply(x, function(s) s %in% 1:100))
    newBClvlsP1 <- paste(as.numeric(apply(as.matrix(1:length(recrnt.spltP1)), 1, function(x) recrnt.spltP1[[x]][is.no.BCP1[[x]]])) - 1, "*", sep = "")
    newBClvlsP1[newBClvlsP1 == "1*"] <- ""
    parents.mat[isP1BC, "P1"] <- apply(as.matrix(1:length(recrnt.spltP1)), 1, function(x) recrnt.spltP1[[x]][!is.no.BCP1[[x]]])
    newsep <- splitlevels[unlist(sapply(ped$Cross[isP1BC], function(x) which(sapply(splitlevels, function(s) grepl(pattern = s, x = x)))[1]))]
    newsep[is.na(newsep)] <- "/"
    parents.mat[isP1BC, "P2"] <- paste(parents.mat[isP1BC, "P2"], newsep,
      paste(newBClvlsP1, apply(as.matrix(1:length(recrnt.spltP1)), 1, function(x) recrnt.spltP1[[x]][!is.no.BCP1[[x]]]), sep = ""),
      sep = ""
    )
  }

  # Add in extra generations for BC parent 2s
  if (sum(isP2BC) > 0) {
    P2.split.pos <- sapply(parents.mat[isP2BC, "P2"], function(x) unlist(gregexpr(pattern = "\\*", x))[unlist(gregexpr(pattern = "\\*", x)) %in% c(2, str_length(x) - 1)])
    recrnt.spltP2 <- lapply(1:length(P2.split.pos), function(x) {
      c(
        substr(parents.mat[isP2BC, "P2"][x], 1, P2.split.pos[x] - 1),
        substr(parents.mat[isP2BC, "P2"][x], P2.split.pos[x] + 1, str_length(parents.mat[isP2BC, "P2"][x]))
      )
    })
    is.no.BCP2 <- lapply(recrnt.spltP2, function(x) sapply(x, function(s) s %in% 1:100))
    newBClvlsP2 <- paste(as.numeric(apply(as.matrix(1:length(recrnt.spltP2)), 1, function(x) recrnt.spltP2[[x]][is.no.BCP2[[x]]])) - 1, "*", sep = "")
    newBClvlsP2[newBClvlsP2 == "1*"] <- ""
    parents.mat[isP2BC, "P2"] <- apply(as.matrix(1:length(recrnt.spltP2)), 1, function(x) recrnt.spltP2[[x]][!is.no.BCP2[[x]]])
    newsep <- splitlevels[unlist(sapply(ped$Cross[isP2BC], function(x) which(sapply(splitlevels, function(s) grepl(pattern = s, x = x)))[1]))]
    newsep[is.na(newsep)] <- "/"
    parents.mat[isP2BC, "P1"] <- paste(parents.mat[isP2BC, "P1"], newsep,
      paste(newBClvlsP2, apply(as.matrix(1:length(recrnt.spltP2)), 1, function(x) recrnt.spltP2[[x]][!is.no.BCP2[[x]]]), sep = ""),
      sep = ""
    )
  }

  # Do the same for all the non BC lines
  all.gens[[1]] <- cbind("ID" = as.character(ped$GID), parents.mat)
  prnts.to.split <- sapply(c(all.gens[[1]][, 2], all.gens[[1]][, 3]), function(x) names(which.min(which(sapply(splitlevels, function(s) grepl(pattern = s, x = x))))))
  donextlevel <- length(unlist(prnts.to.split)) > 0

  # Repeat steps until there are none left to do
  while (donextlevel == T) {
    cat("|", sep = "")
    newIDs <- c(all.gens[[length(all.gens)]][, 2], all.gens[[length(all.gens)]][, 3])
    newIDs <- newIDs[!duplicated(newIDs)]
    newIDs <- newIDs[!newIDs %in% unlist(lapply(all.gens, function(x) x[, 1]))]
    newIDs <- newIDs[sapply(newIDs, function(x) sum(sapply(splitlevels, function(s) grepl(pattern = s, x = x)))) > 0]
    if (length(newIDs) > 0) {
      max.splits <- sapply(newIDs, function(x) names(which.min(which(sapply(splitlevels, function(s) grepl(pattern = s, x = x))))))
      parents <- strsplit(x = newIDs, split = max.splits)
      if (sum(unlist(lapply(parents, length)) > 2) > 0) {
        print("Check pedigrees:")
        print(newIDs[which(unlist(lapply(parents, length)) > 2)])
        stop("Check pedigrees")
      }
      parents.mat <- matrix(NA, nrow = length(parents), ncol = 2, dimnames = list(1:length(parents), c("P1", "P2")))
      for (i in 1:length(parents)) {
        parents.mat[i, ] <- parents[[i]]
      }
    } else {
      parents.mat <- matrix(NA, nrow = 0, ncol = 2, dimnames = list(NULL, c("P1", "P2")))
    }

    if (!length(newIDs) > 0) {
      break
    }

    isP1BC <- (grepl(
      pattern = c("\\*"),
      x = str_sub(parents.mat[, "P1"], 2, 2)
    ) & str_sub(parents.mat[, "P1"], 1, 1) %in% as.character(1:9)) | (grepl(
      pattern = c("\\*"),
      x = str_sub(parents.mat[, "P1"], -2, -2)
    ) & str_sub(parents.mat[, "P1"], -1, -1) %in% as.character(1:9))
    if (sum(isP1BC) > 0) {
      P1.split.pos <- sapply(parents.mat[isP1BC, "P1"], function(x) unlist(gregexpr(pattern = "\\*", x))[unlist(gregexpr(pattern = "\\*", x)) %in% c(2, str_length(x) - 1)])
      recrnt.spltP1 <- lapply(1:length(P1.split.pos), function(x) {
        c(
          substr(parents.mat[isP1BC, "P1"][x], 1, P1.split.pos[x] - 1),
          substr(parents.mat[isP1BC, "P1"][x], P1.split.pos[x] + 1, str_length(parents.mat[isP1BC, "P1"][x]))
        )
      })
      is.no.BCP1 <- lapply(recrnt.spltP1, function(x) sapply(x, function(s) s %in% 1:100))
      newBClvlsP1 <- paste(as.numeric(apply(as.matrix(1:length(recrnt.spltP1)), 1, function(x) recrnt.spltP1[[x]][is.no.BCP1[[x]]])) - 1, "*", sep = "")
      newBClvlsP1[newBClvlsP1 == "1*"] <- ""
      parents.mat[isP1BC, "P1"] <- apply(as.matrix(1:length(recrnt.spltP1)), 1, function(x) recrnt.spltP1[[x]][!is.no.BCP1[[x]]])
      both.prnts <- paste(parents.mat[isP1BC, "P1"], parents.mat[isP1BC, "P2"])
      newsep <- sapply(both.prnts, function(x) splitlevels[which(sapply(splitlevels, function(s) grepl(pattern = s, x = x))) - 1][1])
      newsep[is.na(newsep)] <- "/"
      parents.mat[isP1BC, "P2"] <- paste(parents.mat[isP1BC, "P2"], newsep,
        paste(newBClvlsP1, apply(as.matrix(1:length(recrnt.spltP1)), 1, function(x) recrnt.spltP1[[x]][!is.no.BCP1[[x]]]), sep = ""),
        sep = ""
      )
    }

    isP2BC <- (grepl(
      pattern = c("\\*"),
      x = str_sub(parents.mat[, "P2"], 2, 2)
    ) & str_sub(parents.mat[, "P2"], 1, 1) %in% as.character(1:9)) | (grepl(
      pattern = c("\\*"),
      x = str_sub(parents.mat[, "P2"], -2, -2)
    ) & str_sub(parents.mat[, "P2"], -1, -1) %in% as.character(1:9))

    if (sum(isP2BC) > 0) {
      P2.split.pos <- sapply(parents.mat[isP2BC, "P2"], function(x) unlist(gregexpr(pattern = "\\*", x))[unlist(gregexpr(pattern = "\\*", x)) %in% c(2, str_length(x) - 1)])
      recrnt.spltP2 <- lapply(1:length(P2.split.pos), function(x) {
        c(
          substr(parents.mat[isP2BC, "P2"][x], 1, P2.split.pos[x] - 1),
          substr(parents.mat[isP2BC, "P2"][x], P2.split.pos[x] + 1, str_length(parents.mat[isP2BC, "P2"][x]))
        )
      })
      is.no.BCP2 <- lapply(recrnt.spltP2, function(x) sapply(x, function(s) s %in% 1:100))
      newBClvlsP2 <- paste(as.numeric(apply(as.matrix(1:length(recrnt.spltP2)), 1, function(x) recrnt.spltP2[[x]][is.no.BCP2[[x]]])) - 1, "*", sep = "")
      newBClvlsP2[newBClvlsP2 == "1*"] <- ""
      parents.mat[isP2BC, "P2"] <- apply(as.matrix(1:length(recrnt.spltP2)), 1, function(x) recrnt.spltP2[[x]][!is.no.BCP2[[x]]])
      both.prnts <- paste(parents.mat[isP2BC, "P1"], parents.mat[isP2BC, "P2"])
      newsep <- sapply(both.prnts, function(x) splitlevels[which(sapply(splitlevels, function(s) grepl(pattern = s, x = x))) - 1][1])
      newsep[is.na(newsep)] <- "/"
      parents.mat[isP2BC, "P1"] <- paste(parents.mat[isP2BC, "P1"], newsep,
        paste(newBClvlsP2, apply(as.matrix(1:length(recrnt.spltP2)), 1, function(x) recrnt.spltP2[[x]][!is.no.BCP2[[x]]]), sep = ""),
        sep = ""
      )
    }

    all.gens[[length(all.gens) + 1]] <- cbind("ID" = as.character(newIDs), parents.mat)
    prnts.to.split <- sapply(c(all.gens[[length(all.gens)]][, 2], all.gens[[length(all.gens)]][, 3]), function(x) names(which.min(which(sapply(splitlevels, function(s) grepl(pattern = s, x = x))))))
    donextlevel <- length(unlist(prnts.to.split)) > 0

    if (length(all.gens) > max.depth) {
      donextlevel <- FALSE
    }
  }
  all.gens <- lapply(all.gens, function(x) as.data.frame(x))
  expd.ped <- as.data.frame(rbindlist(all.gens))
  cat(":)\nFinished!")
  return(expd.ped)
}






# Get kinship function----
Get.kin <- function(expnd.ped, # Expanded pedigree in three column format output from Expand.ped function
                    IDs = NULL, # IDs of individual names
                    selhists # Vector of selection history character strings
) {
  library(kinship2)

  if (is.null(IDs)) {
    IDs <- unique(c(expnd.ped$ID, expnd.ped$P1, expnd.ped$P2))
  }
  orphns <- unique(c(expnd.ped$P1, expnd.ped$P2)) # Get list of orphans
  orphns <- orphns[!orphns %in% expnd.ped$ID]
  orphnsP1 <- paste("?", 1:length(orphns), sep = "") # Add in unknown parents for orphans
  orphnsP2 <- paste("?", (length(orphns) + 1):(length(orphns) * 2), sep = "")
  expnd.ped <- data.frame(
    "ID" = c(expnd.ped$ID, orphns),
    "P1" = c(expnd.ped$P1, orphnsP1),
    "P2" = c(expnd.ped$P2, orphnsP2)
  )
  splitlevels <- c(paste("/", 100:3, "/", sep = ""), "//", "/") # Create hierarchy of characters to split purdy strings by
  max.splits <- sapply(expnd.ped$ID, function(x) names(which.min(which(sapply(splitlevels, function(s) grepl(pattern = s, x = x))))))


  # Add inbreeding generations based on selection histories
  selhists.list <- sapply(selhists, function(x) unlist(sapply(1:length(unlist(str_split(x, "-"))), function(z) paste(unlist(str_split(x, "-"))[1:z], collapse = "-")))) # split selection history generations
  selhists.list <- selhists.list[unlist(lapply(selhists.list, function(x) length(x) > 1))] # filter out dodgy selection histories
  selhists.list <- lapply(selhists.list, function(x) x[!x == ""])
  all.selhist.inbreds <- sapply(1:length(selhists.list), function(i) { # Add in intermediate generation in the selection histories
    selhist.inbreds <- data.frame(
      "ID" = NA,
      "P1" = expnd.ped[expnd.ped$ID == names(selhists.list)[i], "P1"],
      "P2" = expnd.ped[expnd.ped$ID == names(selhists.list)[i], "P2"]
    )
    for (g in 1:(length(selhists.list[[i]]) - 1)) {
      selhist.inbreds[g, "ID"] <- selhists.list[[i]][g + 1]
      selhist.inbreds[g + 1, c("P1", "P2")] <- selhists.list[[i]][g + 1]
    }
    selhist.inbreds[nrow(selhist.inbreds), "ID"] <- names(selhists.list)[i]
    selhist.inbreds <- list(selhist.inbreds)
    return(selhist.inbreds)
  })
  all.selhist.inbreds <- as.data.frame(rbindlist(all.selhist.inbreds))


  # Do the same for lines without selection histories and assume 7 generations of selfing
  lines.no.selhist <- expnd.ped[!expnd.ped$ID %in% all.selhist.inbreds$ID, ]
  fill.inselhists <- unlist(apply(lines.no.selhist[, c("P1", "P2")], 1, function(x) paste(paste(x, collapse = "__x__"), "_F1--F2--F3--F4--F5--F6", sep = "")))
  names(fill.inselhists) <- lines.no.selhist$ID
  fillinselhists.list <- sapply(fill.inselhists, function(x) list(sapply(1:length(unlist(str_split(x, "--"))), function(z) paste(unlist(str_split(x, "--"))[1:z], collapse = "--"))))
  all.fillinselhist.inbreds <- sapply(1:length(fillinselhists.list), function(i) {
    selhist.inbreds <- data.frame(
      "ID" = NA,
      "P1" = lines.no.selhist[i, "P1"],
      "P2" = lines.no.selhist[i, "P2"]
    )
    for (g in 1:(length(fillinselhists.list[[i]]) - 1)) {
      selhist.inbreds[g, "ID"] <- fillinselhists.list[[i]][g + 1]
      selhist.inbreds[g + 1, c("P1", "P2")] <- fillinselhists.list[[i]][g + 1]
    }
    selhist.inbreds[nrow(selhist.inbreds), "ID"] <- names(fillinselhists.list)[i]
    selhist.inbreds <- list(selhist.inbreds)
    return(selhist.inbreds)
  })
  all.fillinselhist.inbreds <- as.data.frame(rbindlist(all.fillinselhist.inbreds))

  inbred.ped <- rbind(all.fillinselhist.inbreds, all.selhist.inbreds) # combine selection history lines and full inbred pedigrees
  inbred.ped <- inbred.ped[!duplicated(inbred.ped$ID), ]
  print("Calculating kinship....")
  A <- kinship(id = inbred.ped$ID, dadid = inbred.ped$P1, momid = inbred.ped$P2) # Calculate kinship with kinship2 package
  A <- A[as.character(IDs), as.character(IDs)]
  return(A)
}
