#!/usr/local/bin/Rscript
### Script producing plots from phenotypic analysis
### example to run: ./PlotPhenotype.R

library('flowCore') ### handling FCS files
library('scales') ### for function alpha
library('scatterplot3d') ### for 3D plot
### incorporate function extracting modal population fluorescence values
if (!exists('foo', mode = 'function')) source('PeaksEx.R')
### incorporate function extracting coefficient of variation values
if (!exists('foo', mode = 'function')) source('mCvEx.R')
### incorporate file reordering function
if (!exists("foo", mode="function")) source("ReFiles.R")
### incorporate function calculating distance of point from a line in 2D (plasticity)
if (!exists("foo", mode = "function")) source("PointLine2D.R")
### incorporate function calculating distance of point from a line (plasticity)
if (!exists("foo", mode = "function")) source("PointLine3D.R")

### save list of all possible conditions
### which is further used as a dictionary
cond.ls <- list('_glucose_' = 'Glucose',
                  '_galactose_' = 'Galactose',
                  '_lactose_' = 'Lactose')
### create a similar list as above with shortcuts
short.ls <- list('_glucose_' = 'Glu',
                  '_galactose_' = 'Gal',
                  '_lactose_' = 'Lac')

##################################################
######## DATA LOADING AND PREPROCESSING ##########
##################################################

### set working directories
root.path <- getwd()
### get promoter names and directories with FC datafiles
paths.noRJ <- c(paste0(root.path, '/PhenAnalysis/201019/GMua_env-K12_glucose_all/'),
                paste0(root.path, '/PhenAnalysis/201201/GMua_env-K12_galactose_all/'),
                paste0(root.path, '/PhenAnalysis/201211/GMua_env-K12_lactose_all/'))
paths.wRJ <- c(paste0(root.path, '/PhenAnalysis/201021/GMmv-K12_glucose_all/'),
                paste0(root.path, '/PhenAnalysis/201201/GMmv-K12_galactose_all/'),
                paste0(root.path, '/PhenAnalysis/201211/GMmv-K12_lactose_all/'))
path.chrom <- paste0(root.path, '/PhenAnalysis/201212/mvs_glucose-galactose-lactose/')
path.plas <- c(paths.noRJ, paths.wRJ)

### loop through all datasets and get max kernel densities
### from GFP channel & coefficients of variation
### (i.e. modal population fluorescence & standard deviation / modal population fluorescence)
ndir <- length(unlist(strsplit(path.chrom, split = '/', fixed = T)))
for (path in path.plas) {
  setwd(path)
  if (!file.exists('Peaks.csv')) {
    cat(paste('Analyzing data files in:', path, '\n'))
    ### extract info about last directory
    end.dir <- unlist(strsplit(path, split = '/', fixed = T))[ndir]
    cond <- unlist(strsplit(end.dir, split = '_', fixed = T))[3]
    ### get all files names ending with .fsc
    Files <- Sys.glob('*.fcs')
    nfiles <- length(Files) / 3
    InFiles <- reorder(Files, nfiles)
    ### get modal population fluorescence (max kernel density from GFP channel)
    Fpeaks <- exPeaks(InFiles, (length(Files) / nfiles))
    Fpeaks <- matrix(Fpeaks, nrow = nfiles, ncol = 1, byrow = T)
    rownames(Fpeaks) <- paste(1:nfiles)
    colnames(Fpeaks) <- cond
    write.csv(Fpeaks, 'Peaks.csv', row.names = T)
    ### get coefficient of variation, i.e., standard deviation
    ### from GFP channel / modal population fluorescence
    Fcvs <- exCvs(InFiles, (length(Files) / nfiles))
    Fcvs <- matrix(Fcvs, nrow = nfiles, ncol = 1, byrow = T)
    rownames(Fcvs) <- paste(1:nfiles)
    colnames(Fcvs) <- cond
    write.csv(Fcvs, 'Cvs.csv', row.names = T)
  }
}

setwd(path.chrom)
if (!file.exists('Peaks.csv')) {
  cat(paste('Analyzing data files in:', path.chrom, '\n'))
  ### extract info about last directory
  end.dir <- unlist(strsplit(path, split = '/', fixed = T))[ndir]
  ### get all files names ending with .fsc
  InFiles <- Sys.glob('*.fcs')
  nfiles <- length(InFiles) / 3
  ### get modal population fluorescence (max kernel density from GFP channel)
  Fpeaks <- exPeaks(InFiles, (length(InFiles) / nfiles))
  Fpeaks <- matrix(Fpeaks, nrow = nfiles, ncol = 1, byrow = T)
  rownames(Fpeaks) <- paste(1:nfiles)
  write.csv(Fpeaks, 'Peaks.csv', row.names = T)
  ### get coefficient of variation, i.e., standard deviation
  ### from GFP channel / modal population fluorescence
  Fcvs <- exCvs(InFiles, (length(InFiles) / nfiles))
  Fcvs <- matrix(Fcvs, nrow = nfiles, ncol = 1, byrow = T)
  rownames(Fcvs) <- paste(1:nfiles)
  write.csv(Fcvs, 'Cvs.csv', row.names = T)
}

setwd(root.path)
### sort content of mixed libraries into separate variables
mutTOT.ls <- list()
segTOT.ls <- list()
natTOT.ls <- list()
segTXN.ls <- list()
natTXN.ls <- list()
for (path in paths.noRJ) {
  ### extract info about last directory
  end.dir <- unlist(strsplit(path, split = '/', fixed = T))[ndir]
  cond <- unlist(strsplit(end.dir, split = '_', fixed = T))[3]
  ### load modal fluorescence and coefficient of variation data
  data.m <- read.csv(paste0(path, 'Peaks.csv'), header = T)
  data.m <- data.m[, -1]
  data.n <- read.csv(paste0(path, 'Cvs.csv'), header = T)
  data.n <- data.n[, -1]
  ### extract fluorescence data from mutagenesis
  ### without RiboJ including controls
  mut <- cbind(c(data.m[1:8], data.m[10:18], data.m[20:21], data.m[23:31], data.m[59], data.m[91], data.m[60], data.m[92], data.m[90]),
              c(data.n[1:8], data.n[10:18], data.n[20:21], data.n[23:31], data.n[59], data.n[91], data.n[60], data.n[92], data.n[90]))
  snp <- read.csv('PhenAnalysis/SimMatrixMutagenesisUA.csv', header = T)
  snp <- snp[, 1:2]
  snp <- snp[-1,]
  names <- snp[, 1][which(snp[, 2] > 0 & snp[, 2] < 4)]
  names <- c(names, 'MG1655.TOT', 'MG1655.TXN', 'pUA66', 'pMV001', 'plpp')
  vals <- snp[, 2][which(snp[, 2] > 0 & snp[, 2] < 4)]
  vals <- c(vals, 0, rep(NA, 4))
  df.mut <- cbind(names, vals, mut)
  mutTOT.ls[[paste0('_', cond, '_')]] <- as.data.frame(df.mut)
  ### extract fluorescence data from segr. variants
  ### in MG1655 without RiboJ including controls
  seg <- cbind(c(data.m[41:48], data.m[33], data.m[49:59], data.m[91], data.m[60], data.m[92], data.m[90]),
              c(data.n[41:48], data.n[33], data.n[49:59], data.n[91], data.n[60], data.n[92], data.n[90]))
  names <- c('SC303', 'SC305', 'SC307', 'SC312', 'SC316',
            'SC332', 'SC336', 'SC357', 'SC358', 'SC370',
            'SC371', 'SC381', 'SC382', 'SC401', 'SC400',
            'SC418', 'SC429', 'SC476', 'SC480', 'MG1655.TOT',
            'MG1655.TXN', 'pUA66', 'pMV001', 'plpp')
  df.seg <- cbind(names, seg)
  segTOT.ls[[paste0('_', cond, '_')]] <- as.data.frame(df.seg)
  ### extract fluorescence data from segr. variants
  ### in native background without RiboJ including controls
  nat <- cbind(c(data.m[61], data.m[34], data.m[63:64], data.m[37:40], data.m[59], data.m[91], data.m[60], data.m[92], data.m[90]),
              c(data.n[61], data.n[34], data.n[63:64], data.n[37:40], data.n[59], data.n[91], data.n[60], data.n[92], data.n[90]))
  names <- c('SC312', 'SC358', 'SC400', 'SC418',
            'SC312.TOT', 'SC358.TOT', 'SC400.TOT', 'SC418.TOT',
            'MG1655.TOT', 'MG1655.TXN', 'pUA66', 'pMV001', 'plpp')
  df.nat <- cbind(names, nat)
  natTOT.ls[[paste0('_', cond, '_')]] <- as.data.frame(df.nat)
  ### extract fluorescence data from segr. variants
  ### in MG1655 with RiboJ including controls
  seg <- cbind(c(data.m[73:89], data.m[59], data.m[91], data.m[60], data.m[92], data.m[90]),
              c(data.n[73:89], data.n[59], data.n[91], data.n[60], data.n[92], data.n[90]))
  names <- c('SC303', 'SC305', 'SC307', 'SC312', 'SC316',
            'SC332', 'SC353', 'SC357', 'SC370', 'SC371',
            'SC381', 'SC393', 'SC400', 'SC418', 'SC429',
            'SC476', 'SC480', 'MG1655.TOT', 'MG1655.TXN',
            'pUA66', 'pMV001', 'plpp')
  df.seg <- cbind(names, seg)
  segTXN.ls[[paste0('_', cond, '_')]] <- as.data.frame(df.seg)
  ### extract fluorescence data from segr. variants
  ### in native background with RiboJ including controls
  nat <- cbind(c(data.m[93:96], data.m[69:72], data.m[59], data.m[91], data.m[60], data.m[92], data.m[90]),
              c(data.n[93:96], data.n[69:72], data.n[59], data.n[91], data.n[60], data.n[92], data.n[90]))
  names <- c('SC312', 'SC358', 'SC400', 'SC418',
            'SC312.TXN', 'SC358.TXN', 'SC400.TXN', 'SC418.TXN',
            'MG1655.TOT', 'MG1655.TXN', 'pUA66', 'pMV001', 'plpp')
  df.nat <- cbind(names, nat)
  natTXN.ls[[paste0('_', cond, '_')]] <- as.data.frame(df.nat)
}

mutTXN.ls <- list()
for (path in paths.wRJ) {
  ### extract info about last directory
  end.dir <- unlist(strsplit(path, split = '/', fixed = T))[ndir]
  cond <- unlist(strsplit(end.dir, split = '_', fixed = T))[2]
  ### load modal fluorescence and coefficient of variation data
  data.m <- read.csv(paste0(path, 'Peaks.csv'), header = T)
  data.m <- data.m[, -1]
  data.n <- read.csv(paste0(path, 'Cvs.csv'), header = T)
  data.n <- data.n[, -1]
  ### extract fluorescence data from mutagenesis
  ### without RiboJ including controls
  mut <- cbind(c(data.m[1:85], data.m[87:91], data.m[93:94], data.m[96], data.m[95], data.m[86]),
              c(data.n[1:85], data.n[87:91], data.n[93:94], data.n[96], data.n[95], data.n[86]))
  snp <- read.csv('PhenAnalysis/SimMatrixMutagenesisMV.csv', header = T)
  snp <- snp[, 1:2]
  snp <- snp[-1,]
  names <- snp[, 1][which(snp[, 2] > 0 & snp[, 2] < 4)]
  names <- c(names, 'MG1655.TXN', 'pMV001', 'plpp')
  vals <- snp[, 2][which(snp[, 2] > 0 & snp[, 2] < 4)]
  vals <- c(vals, 0, rep(NA, 2))
  df.mut <- cbind(names, vals, mut)
  mutTXN.ls[[paste0('_', cond, '_')]] <- as.data.frame(df.mut)
}

### extract fluorescence data from seg.
### variants in MG1655 chromosome
chrom.ls <- list()
data.m <- read.csv(paste0(path.chrom, 'Peaks.csv'), header = T)
data.m <- data.m[, -1]
data.n <- read.csv(paste0(path.chrom, 'Cvs.csv'), header = T)
data.n <- data.n[, -1]
glu <- cbind(c(data.m[1], data.m[4], data.m[7], data.m[10], data.m[13], data.m[16], data.m[22], data.m[19]),
            c(data.n[1], data.n[4], data.n[7], data.n[10], data.n[13], data.n[16], data.n[22], data.n[19]))
gal <- cbind(c(data.m[2], data.m[5], data.m[8], data.m[11], data.m[14], data.m[17], data.m[23], data.m[20]),
            c(data.n[2], data.n[5], data.n[8], data.n[11], data.n[14], data.n[17], data.n[23], data.n[20]))
lac <- cbind(c(data.m[3], data.m[6], data.m[9], data.m[12], data.m[15], data.m[18], data.m[24], data.m[21]),
            c(data.n[3], data.n[6], data.n[9], data.n[12], data.n[15], data.n[18], data.n[24], data.n[21]))
names <- c('MG1655', 'SC312', 'SC400', 'SC418', 'SC358', 'MG1655.NEG', 'pUA66', 'plpp')
chrom.ls[[paste0('_glucose_')]] <- as.data.frame(cbind(names, glu))
chrom.ls[[paste0('_galactose_')]] <- as.data.frame(cbind(names, gal))
chrom.ls[[paste0('_lactose_')]] <- as.data.frame(cbind(names, lac))

### set offsets for chromosomal and mutagenesis only plates
offsetTXN.m <- vector()
offsetTXN.n <- vector()
offsetCH.m <- vector()
offsetCH.n <- vector()
offnames <- vector()
### change nseg variable to number of segregating variants instead
for (cond in names(cond.ls)) {
  ### set first dataset to calculate offsets from
  ref <- segTOT.ls[[cond]]
  mut <- mutTXN.ls[[cond]]
  ### save values for offset calculation
  ### (modal population fluorescence)
  pos.ref.m <- as.numeric(ref[, 2][which(grepl('plpp', ref[, 1]))])
  pos.mut.m <- as.numeric(mut[, 3][which(grepl('plpp', mut[, 1]))])
  neg.ref.m <- as.numeric(ref[, 2][which(grepl('pMV001', ref[, 1]))])
  neg.mut.m <- as.numeric(mut[, 3][which(grepl('pMV001', mut[, 1]))])
  mg1655.ref.m <- as.numeric(ref[, 2][which(grepl('MG1655.TXN', ref[, 1]))])
  mg1655.mut.m <- as.numeric(mut[, 3][which(grepl('MG1655.TXN', mut[, 1]))])
  ### save values for offset calculation
  ### (coefficient of variantion)
  pos.ref.n <- as.numeric(ref[, 3][which(grepl('plpp', ref[, 1]))])
  pos.mut.n <- as.numeric(mut[, 4][which(grepl('plpp', mut[, 1]))])
  neg.ref.n <- as.numeric(ref[, 3][which(grepl('pMV001', ref[, 1]))])
  neg.mut.n <- as.numeric(mut[, 4][which(grepl('pMV001', mut[, 1]))])
  mg1655.ref.n <- as.numeric(ref[, 3][which(grepl('MG1655.TXN', ref[, 1]))])
  mg1655.mut.n <- as.numeric(mut[, 4][which(grepl('MG1655.TXN', mut[, 1]))])
  ### calculate offset between the two libraries
  ### (both fluorescence and coefficient of variation)
  off.m <- mean(c((pos.mut.m - pos.ref.m), (neg.mut.m - neg.ref.m), (mg1655.mut.m - mg1655.ref.m)), na.rm = T)
  offsetTXN.m <- c(offsetTXN.m, off.m)
  off.n <- mean(c((pos.mut.n - pos.ref.n), (neg.mut.n - neg.ref.n), (mg1655.mut.n - mg1655.ref.n)), na.rm = T)
  offsetTXN.n <- c(offsetTXN.n, off.n)
  ### set second dataset to calculate offsets from
  chr <- chrom.ls[[cond]]
  ### refresh values for offset calculation
  ### (modal population fluorescence)
  pos.mut.m <- as.numeric(chr[, 2][which(grepl('plpp', chr[, 1]))])
  neg.ref.m <- as.numeric(ref[, 2][which(grepl('pUA66', ref[, 1]))])
  neg.mut.m <- as.numeric(chr[, 2][which(grepl('pUA66', chr[, 1]))])
  ### refresh values for offset calculation
  ### (coefficient of variantion)
  pos.mut.n <- as.numeric(chr[, 3][which(grepl('plpp', chr[, 1]))])
  neg.ref.n <- as.numeric(ref[, 3][which(grepl('pUA66', ref[, 1]))])
  neg.mut.n <- as.numeric(chr[, 3][which(grepl('pUA66', chr[, 1]))])
  ### calculate offset between the two libraries
  ### (both fluorescence and coefficient of variation)
  off.m <- mean(c((pos.mut.m - pos.ref.m), (neg.mut.m - neg.ref.m)), na.rm = T)
  offsetCH.m <- c(offsetCH.m, off.m)
  off.n <- mean(c((pos.mut.n - pos.ref.n), (neg.mut.n - neg.ref.n)), na.rm = T)
  offsetCH.n <- c(offsetCH.n, off.n)
  ### save conditions to assign appropriate names to offsets
  offnames <- c(offnames, cond)
}
names(offsetTXN.m) <- offnames
names(offsetTXN.n) <- offnames
names(offsetCH.m) <- offnames
names(offsetCH.n) <- offnames

### save names shared by UA & MV random variants
### that differ in sequence
mut.bad <- c('038', '044', '048', '070', '136')

##################################################
################ PLOTTING BEGINS #################
##################################################

##################################################
################# MAPPING PLOTS ##################
################# SUPP FIGURE 3 ##################
##################################################

### plot Supp Figure 3a (mapping the effect of random SNPs
### relative to MG1655 variant to position within promoter sequence)
cat(paste('Producing Supplementary Figure 3a.1\n'))
pdf(file = 'SupplementaryFigure_3a.1.pdf', width = 10, height = 5)
  par(mfcol = c(2, 3),
      mar = c(2.1, 4, 1.5, 1),
      las = 1)

  ### obtain data about SNP positions for each random variant
  snp_mapTXN <- read.csv('PhenAnalysis/1SNPmapMV.csv', header = T)
  snp_mapTOT <- read.csv('PhenAnalysis/1SNPmapUA.csv', header = T)
  lenTXN <- length(snp_mapTXN[, 2])
  lenTOT <- length(snp_mapTOT[, 2])
  ### obtain info about promoter annotations (TF and so on)
  anns <- read.csv('PhenAnalysis/AnnotationsBasic.csv', header = T)
  ### loop through all three envrionments
  for (cond in names(cond.ls)) {
    ### save modal fluorescence values for random variants and the MG1655 varaint
    plTXN <- as.numeric(mutTXN.ls[[cond]][, 3][which(as.numeric(mutTXN.ls[[cond]][, 2]) == 1)])
    names(plTXN) <- mutTXN.ls[[cond]][, 1][which(as.numeric(mutTXN.ls[[cond]][, 2]) == 1)]
    mgTXN <- as.numeric(mutTXN.ls[[cond]][, 3][which(as.numeric(mutTXN.ls[[cond]][, 2]) == 0)])
    plTOT <- as.numeric(mutTOT.ls[[cond]][, 3][which(as.numeric(mutTOT.ls[[cond]][, 2]) == 1)])
    names(plTOT) <- mutTOT.ls[[cond]][, 1][which(as.numeric(mutTOT.ls[[cond]][, 2]) == 1)]
    mgTOT <- as.numeric(mutTOT.ls[[cond]][, 3][which(as.numeric(mutTOT.ls[[cond]][, 2]) == 0)])
    ### plotting
    start <- anns[length(anns[, 1]), 3]
    end <- anns[length(anns[, 1]), 4]
    positions <- seq(0 - floor(start / 50) * 50, floor((end - start) / 50) * 50, 50)
    if (cond == names(cond.ls)[1]) {
      par(mar = c(1, 4, 2, 0))
    } else if (cond == names(cond.ls)[2]) {
      par(mar = c(1, 3, 2, 1))
    } else {
      par(mar = c(1, 2, 2, 2))
    }
    plot(3, 3, type = 'n', xlim = c(1, snp_mapTXN[lenTXN, 2]),
          ylim = c(2, 5), xaxt = 'n',
          main = cond.ls[[cond]],
          cex.main = 1.5, xlab = '', ylab = '')
    for (an in 1:length(anns[, 1])) {
      rect(anns[an, 3], 2, anns[an, 4], 5, col = alpha(anns[an, 5], 0.25), border = NA)
    }
    arrows(1, mgTXN, snp_mapTXN[lenTXN, 2], mgTXN, length = 0)
    for (p in 1:length(plTXN)) {
      if (!is.na(plTXN[p])) {
        if (names(plTXN[p]) %in% names(plTOT)) {
          arrows(snp_mapTXN[p, 2], mgTXN, snp_mapTXN[p, 2], plTXN[p],
                length = 0, col = alpha('black', 0.8), lwd = 1.5)
        } else {
          arrows(snp_mapTXN[p, 2], mgTXN, snp_mapTXN[p, 2], plTXN[p],
                length = 0, col = alpha('cyan', 0.8), lwd = 1.5)
        }
      }
    }
    mtext('pMV001', side = 3, line = -1.5, adj = 0.05)
    axis(side = 1, at = seq(start - floor(start / 50) * 50, floor(end / 50) * 50 + 10 * floor((end - start) / 50) + 1, 50),
          labels = rep('', length(positions)), cex.axis = 1)
    if (cond == names(cond.ls)[1]) {
      title(ylab = 'Modal fluorescence (log10, a.u.)', line = 2.5, cex.lab = 1.2)
      par(mar = c(3, 4, 0, 0))
    } else if (cond == names(cond.ls)[2]) {
      par(mar = c(3, 3, 0, 1))
    } else {
      par(mar = c(3, 2, 0, 2))
    }
    plot(3, 3, type = 'n', xlim = c(1, snp_mapTOT[lenTOT, 2]),
          ylim = c(2, 5), xaxt = 'n',
          cex.main = 1.5, xlab = '', ylab = '')
    for (an in 1:length(anns[, 1])) {
      rect(anns[an, 3], 2, anns[an, 4], 5, col = alpha(anns[an, 5], 0.25), border = NA)
    }
    arrows(1, mgTOT, snp_mapTOT[lenTOT, 2], mgTOT, length = 0)
    for (p in 1:length(plTOT)) {
      if (!is.na(plTOT[p])) {
        arrows(snp_mapTOT[p, 2], mgTOT, snp_mapTOT[p, 2], plTOT[p],
              length = 0, col = alpha('black', 0.8), lwd = 1.5)
      }
    }
    if (cond == names(cond.ls)[1]) {
      title(ylab = 'Modal fluorescence (log10, a.u.)', line = 2.5, cex.lab = 1.2)
    }
    mtext('pUA66', side = 3, line = -1.5, adj = 0.05)
    axis(side = 1, at = seq(start - floor(start / 50) * 50, floor(end / 50) * 50 + 10 * floor((end - start) / 50) + 1, 50),
          labels = positions, cex.axis = 1)
  }

dev.off()

### produce legend for Supp Figure 3a (mapping the effect of random SNPs
### relative to MG1655 variant to position within promoter sequence)
cat(paste('Producing Supplementary Figure 3a.2\n'))
pdf(file = 'SupplementaryFigure_3a.2.pdf', width = 12, height = 3)

  plot(NULL, axes = F, ann = F, xlim = c(0, 1), ylim = c(0, 1))
  legend('top', legend = c('open reading frame', '-35 or -10 element', 'inducer TF', 'repressor TF'),
          pch = 15, col = c(alpha('yellow', 0.25), alpha('grey', 0.25), alpha('green', 0.25), alpha('red', 0.25)),
          title = 'Annotations', bg = 'white', cex = 1.54, horiz = T)

dev.off()

### plot Supp Figure 3b (comparing differences in fluorescence
### relative to MG1655 between plasmids)
cat(paste('Producing Supplementary Figure 3b\n'))
pdf(file = 'SupplementaryFigure_3b.pdf', width = 8, height = 3)
  par(mfrow = c(1, 3),
      mar = c(4, 3.5, 2, 0.5),
      las = 1)

      for (cond in names(cond.ls)) {
        ### save modal fluorescence values for random variants
        muTXN <- as.numeric(mutTXN.ls[[cond]][, 3][which(as.numeric(mutTXN.ls[[cond]][, 2]) > 0 & as.numeric(mutTXN.ls[[cond]][, 2]) < 4)])
        names(muTXN) <- mutTXN.ls[[cond]][, 1][which(as.numeric(mutTXN.ls[[cond]][, 2]) > 0 & as.numeric(mutTXN.ls[[cond]][, 2]) < 4)]
        muTXN <- muTXN - offsetTXN.m[cond]
        muTOT <- as.numeric(mutTOT.ls[[cond]][, 3][which(as.numeric(mutTOT.ls[[cond]][, 2]) > 0 & as.numeric(mutTOT.ls[[cond]][, 2]) < 4)])
        names(muTOT) <- mutTOT.ls[[cond]][, 1][which(as.numeric(mutTOT.ls[[cond]][, 2]) > 0 & as.numeric(mutTOT.ls[[cond]][, 2]) < 4)]
        muTXN <- muTXN[order(names(muTXN))]
        muTOT <- muTOT[order(names(muTOT))]
        ### save modal fluorescence values for segregating variants in the MG1655 background
        segTXN <- as.numeric(segTXN.ls[[cond]][, 2][which(grepl('SC', segTXN.ls[[cond]][, 1]))])
        segTXN <- c(segTXN, as.numeric(segTXN.ls[[cond]][, 2][which(grepl('MG1655.TXN', segTXN.ls[[cond]][, 1]))]))
        names(segTXN) <- c(segTXN.ls[[cond]][, 1][which(grepl('SC', segTXN.ls[[cond]][, 1]))], 'MG1655')
        segTOT <- as.numeric(segTOT.ls[[cond]][, 2][which(grepl('SC', segTOT.ls[[cond]][, 1]))])
        segTOT <- c(segTOT, as.numeric(segTOT.ls[[cond]][, 2][which(grepl('MG1655.TOT', segTOT.ls[[cond]][, 1]))]))
        names(segTOT) <- c(segTOT.ls[[cond]][, 1][which(grepl('SC', segTOT.ls[[cond]][, 1]))], 'MG1655')
        ### combine segregating and random variants into single variables
        TXN <- c(muTXN, segTXN)
        TOT <- c(muTOT, segTOT)
        ### filter out variants not shared by both plasmid systems
        txn <- c()
        for (p in 1:length(TXN)) {
          if (names(TXN[p]) %in% names(TOT) && !(names(TXN[p]) %in% mut.bad)) {
            txn <- c(txn, TXN[p])
          }
        }

        tot <- c()
        nam <- c()
        for (p in 1:length(TOT)) {
          if (names(TOT[p]) %in% names(TXN) && !(names(TOT[p]) %in% mut.bad)) {
            tot <- c(tot, TOT[p])
            nam <- c(nam, names(TOT[p]))
          }
        }

        ct <- cor.test(txn, tot, method = 'pearson')
        plot(txn, tot, pch = 16, xlab = '', ylab = '',
              col = alpha('black', 0.3),
              main = cond.ls[[cond]],
              xlim = c(2, 5.2), ylim = c(2, 5.2))
        text(x = 3.6, y = 5.2, labels = paste('rho =', round(ct$estimate, digits = 3)))
        text(x = 3.6, y = 5, labels = paste('p =', signif(ct$p.value, digits = 2)))
        abline(0, 1, col = 'blue', lty = 3)
        title(xlab = 'Modal fluorescence in pMV001', line = 2.5)
        if (cond == names(cond.ls)[1]) {
          title(ylab = 'Modal fluorescence in pUA66', line = 2.5)
          par(mar = c(4, 2.5, 2, 1.5))
        } else if (cond == names(cond.ls)[2]) {
          par(mar = c(4, 1.5, 2, 2.5))
        }

      }

dev.off()

### sorting effect sizes of changes in fluorescence
### based on whether a SNP is inside (annTrue) or
### or outside (annFalse) TF or RNAP binding site
annTrue <- c()
annFalse <- c()
### obtain data about SNP positions for each random variant
snp_map <- read.csv('PhenAnalysis/1SNPmapMV.csv', header = T)
### obtain info about promoter annotations (TF and so on)
anns <- read.csv('PhenAnalysis/AnnotationsBasic.csv', header = T)
### loop through all environments
### (except for the exclusion mentioned above)
for (cond in names(cond.ls)) {
  ### save modal fluorescence values for random variants and the MG1655 varaint
  pl <- as.numeric(mutTXN.ls[[cond]][, 3][which(as.numeric(mutTXN.ls[[cond]][, 2]) == 1)])
  mg <- as.numeric(mutTXN.ls[[cond]][, 3][which(as.numeric(mutTXN.ls[[cond]][, 2]) == 0)])
  ### calculate fold changes in fluorescence caused by individual SNPs
  sizes <- pl / mg
  ### loop through all random variants
  for (snp in 1:length(snp_map[, 2])) {
    hit <- 0
    for (ann in rownames(anns)) {
      if (grepl('TF', anns[ann, 2]) || grepl('sigma70', anns[ann, 2])) {
        if (snp_map[snp, 2] >= anns[ann, 3] && snp_map[snp, 2] <= anns[ann, 4]) {
          hit <- 1
        }
      }
    }
    ### if the variant has a SNP inside TF or RNAP binding site
    ###  assing the fluorescence fold change value to 'annTrue' variable
    if (hit == 1) {
      annTrue <- c(annTrue, sizes[snp])
    ### if the SNP is outside of any TF or RNAP binding site
    ### assign the fluorescence fold change value to 'annFalse' variable
    } else {
      annFalse <- c(annFalse, sizes[snp])
    }
  }
}

##################################################
############### FLUORESCENCE LEVELS ##############
#################### FIGURE 3 ####################
##################################################

### produce Figure 3a (cumulative fluorescence)
cat(paste('Producing Figure 3a\n'))
pdf(file = 'Figure_3a.pdf', width = 6, height = 2)
  par(mfrow = c(1, 3),
      las = 1)

  cat(paste('p-value indicating significant differences in transcription rate variance between random and segregating variants in:\n'))
  ### loop through all environments
  for (n in 1:3) {
    cond <- names(cond.ls)[n]
    if (n == 1) {
      par(mar = c(3.5, 3.5, 1, 0.5))
    } else {
      if (n == 2) {
        par(mar = c(3.5, 2.5, 1, 1.5))
      } else {
        par(mar = c(3.5, 1.5, 1, 2.5))
      }
    }
    ### get modal fluorescence values for random variants
    ### and use offset so they are comparable to seg. variants
    mut1 <- as.numeric(mutTXN.ls[[cond]][, 3][which(as.numeric(mutTXN.ls[[cond]][, 2]) == 1)])
    mut1 <- mut1 - offsetTXN.m[cond]
    mut2 <- as.numeric(mutTXN.ls[[cond]][, 3][which(as.numeric(mutTXN.ls[[cond]][, 2]) == 2)])
    mut2 <- mut2 - offsetTXN.m[cond]
    mut3 <- as.numeric(mutTXN.ls[[cond]][, 3][which(as.numeric(mutTXN.ls[[cond]][, 2]) == 3)])
    mut3 <- mut3 - offsetTXN.m[cond]
    ### get modal fluorescence values for seg. variants in MG1655
    seg <- as.numeric(segTXN.ls[[cond]][, 2][which(grepl('SC', segTXN.ls[[cond]][, 1]))])
    seg <- c(seg, as.numeric(segTXN.ls[[cond]][, 2][which(grepl('MG1655.TXN', segTXN.ls[[cond]][, 1]))]))
    ### check whether increase and decrease in fluorescence
    ### as compared to MG1655 variant are equally probable
    mg1655 <- seg[length(seg)]
    m <- c(mut1, mut2, mut3)
    bin <- binom.test(x = length(m[which(m > mg1655)]), n = length(m))
    ### compare values from random and segregating variants
    dat <- c(m, seg)
    grp <- c(rep('A', length(m)), rep('B', length(seg)))
    flig <- fligner.test(x = dat, g = grp)
    cat(paste('\t', cond.ls[[cond]], ':', signif(flig$p.value, digits = 3), '\n'))
    ### plotting
    plot(sort(mut1), (1:length(mut1) / length(mut1)),
          pch = 16, type = 'o', cex = 0.75,
          xlim = c(2, 5.5), ylim = c(0, 1),
          col = alpha('blue', 0.75),
          xlab = '', ylab = '', main = '')
    points(sort(mut2), (1:length(mut2) / length(mut2)),
            col = alpha('gold', 0.75), pch = 16, type = 'o', cex = 0.75)
    points(sort(mut3), (1:length(mut3) / length(mut3)),
            col = alpha('red', 0.75), pch = 16, type = 'o', cex = 0.75)
    points(sort(seg), (1:length(seg) / length(seg)),
            col = alpha('black', 0.75), pch = 16, type = 'o', cex = 0.75)
    points(seg[length(seg)], (which(sort(seg) == seg[length(seg)]) / length(seg)),
            bg = 'white', pch = 21, cex = 0.8)
    if (bin$p.value <= 0.05) {
      mtext(signif(bin$p.value, digits = 2), side = 3, line = -1, adj = 0.05, cex = 0.5, font = 4)
    } else {
      mtext(signif(bin$p.value, digits = 2), side = 3, line = -1, adj = 0.05, cex = 0.5)
    }
    mtext(cond.ls[[cond]], side = 1, line = -1, adj = 0.95, cex = 0.75, font = 4)
    if (cond == names(cond.ls)[1]) {
      title(ylab = 'Cumulative density', line = 2.5)
    }
    title(xlab = 'Modal fluorescence (log10, a.u.)', line = 2.25)
  }

dev.off()

cols <- c('black', 'blue', 'gold', 'red')
lty.cond <- c(3, 2, 6)
### produce Figure 3b (fluorescence stdev)
cat(paste('Producing Figure 3b\n'))
pdf(file = 'Figure_3b.pdf', width = 6, height = 3)
  par(mar = c(1, 4, 1, 1),
      mfrow = c(1, 2),
      las = 1)

  for (n in 1:3) {
    cond <- names(cond.ls)[n]
    ### get modal fluorescence values for random variants
    ### and use offset so they are comparable to seg. variants
    mut1 <- as.numeric(mutTXN.ls[[cond]][, 3][which(as.numeric(mutTXN.ls[[cond]][, 2]) == 1)])
    mut2 <- as.numeric(mutTXN.ls[[cond]][, 3][which(as.numeric(mutTXN.ls[[cond]][, 2]) == 2)])
    mut3 <- as.numeric(mutTXN.ls[[cond]][, 3][which(as.numeric(mutTXN.ls[[cond]][, 2]) == 3)])
    ### get modal fluorescence values for seg. variants in MG1655
    seg <- as.numeric(segTXN.ls[[cond]][, 2][which(grepl('SC', segTXN.ls[[cond]][, 1]))])
    seg <- c(seg, as.numeric(segTXN.ls[[cond]][, 2][which(grepl('MG1655.TXN', segTXN.ls[[cond]][, 1]))]))
    ### get stdev values for all groups
    sdsm <- c(sd(seg), sd(mut1), sd(mut2), sd(mut3))
    ### save number of samples in each group
    nsam <- c(length(seg), length(mut1), length(mut2), length(mut3))
    ### plotting
    if (cond == names(cond.ls)[1]) {
      plot(1:4, sdsm, col = 'grey', pch = 16,
          type = 'b', xlim = c(0.5, 4.5), ylim = c(0, 0.6),
          xaxt = 'n', ylab = '', lty = lty.cond[n])
    } else {
      points(1:4, sdsm, col = 'grey', pch = 16, type = 'b', lty = lty.cond[n])
      text(x = 1:4, y = rep(0.6, 4), labels = nsam, cex = 0.75)
    }
    points(1:4, sdsm, col = cols, pch = 16)
  }
  title(ylab = 'Stdev of modal fluorescence', line = 2.5)

  par(mar = c(0.5, 4.5, 0.5, 3.5))
  plot(NULL, axes = F, ann = F, xlim = c(0, 1), ylim = c(0, 1))
  legend('top', legend = c('Segregating', '1 mutation', '2 mutations', '3 mutations'), pch = 16, col = cols, box.col = 'white')
  legend('bottom', legend = c('Glucose', 'Galactose', 'Lactose'), lty = lty.cond, col = 'grey', box.col = 'white')

dev.off()

##################################################
################ PLASTICITY LEVELS ###############
############# FIGURE 4 & SUPP FIG. 4 #############
##################################################

### set values for ylim and xlim in plotting
lims <- c(2, 5.5)
### produce Figure 4a (2D fluorescence levels)
cat(paste('Producing Figure 4a\n'))
pdf(file = "Figure_4a.pdf", width = 6, height = 2)
  par(mfrow = c(1, 3),
      las = 1)

  ### get modal fluorescence values for all variants in MG1655
  TXN <- list()
  for (cond in names(cond.ls)) {
    for (n in 1:3) {
      name <- paste0(cond, ':random_', n)
      hold <- as.numeric(mutTXN.ls[[cond]][, 3][which(as.numeric(mutTXN.ls[[cond]][, 2]) == n)])
      hold <- hold - offsetTXN.m[cond]
      TXN[[name]] <- hold
    }
    name <- paste0(cond, ':segregating')
    TXN[[name]] <- c(as.numeric(segTXN.ls[[cond]][, 2][which(grepl('SC', segTXN.ls[[cond]][, 1]))]),
                      as.numeric(segTXN.ls[[cond]][, 2][which(grepl('MG1655.TXN', segTXN.ls[[cond]][, 1]))]))
  }
  ### define fluorescence comparisons
  plasts <- c(paste0(names(cond.ls)[1], ':', names(cond.ls)[2]),
              paste0(names(cond.ls)[1], ':', names(cond.ls)[3]),
              paste0(names(cond.ls)[2], ':', names(cond.ls)[3]))
  for (plast in plasts) {
    conds <- unlist(strsplit(plast, ':'))
    ### get fluorescence values for environment 1
    cond1.mut1 <- unlist(TXN[which(grepl(conds[1], names(TXN)) & grepl('random_1', names(TXN)))])
    cond1.mut2 <- unlist(TXN[which(grepl(conds[1], names(TXN)) & grepl('random_2', names(TXN)))])
    cond1.mut3 <- unlist(TXN[which(grepl(conds[1], names(TXN)) & grepl('random_3', names(TXN)))])
    cond1.seg <- unlist(TXN[which(grepl(conds[1], names(TXN)) & grepl('segregating', names(TXN)))])
    ### get fluorescence levels for environment 2
    cond2.mut1 <- unlist(TXN[which(grepl(conds[2], names(TXN)) & grepl('random_1', names(TXN)))])
    cond2.mut2 <- unlist(TXN[which(grepl(conds[2], names(TXN)) & grepl('random_2', names(TXN)))])
    cond2.mut3 <- unlist(TXN[which(grepl(conds[2], names(TXN)) & grepl('random_3', names(TXN)))])
    cond2.seg <- unlist(TXN[which(grepl(conds[2], names(TXN)) & grepl('segregating', names(TXN)))])
    if (plast == plasts[1]) {
      par(mar = c(3, 3.5, 1, 0.5))
    } else {
      if (plast == plasts[2]) {
        par(mar = c(3, 3.5, 1, 0.5))
      } else {
        par(mar = c(3, 3.5, 1, 0.5))
      }
    }
    ### plotting
    plot(cond1.mut1, cond2.mut1, xlim = lims, ylim = lims,
          pch = 16, col = alpha('blue', 0.75), xlab = '',
          ylab = '', main = '', cex = 0.75)
    points(cond1.seg, cond2.seg, pch = 16, col = alpha('black', 0.75), cex = 0.75)
    points(cond1.mut2, cond2.mut2, pch = 16, col = alpha('gold', 0.75), cex = 0.75)
    points(cond1.mut3, cond2.mut3, pch = 16, col = alpha('red', 0.75), cex = 0.75)
    points(cond1.seg[length(cond1.seg)], cond2.seg[length(cond2.seg)], pch = 21, bg = 'white', cex = 0.75)
    abline(0, 1, lty = 3, col = 'blue')
    abline(abs(cond1.seg[length(cond1.seg)] - cond2.seg[length(cond2.seg)]),
            1, lty = 3, col = 'red')
    title(xlab = cond.ls[[conds[1]]], line = 2)
    title(ylab = cond.ls[[conds[2]]], line = 2.5)
  }

dev.off()

### produce Figure 4b (cumulative 2D plasticity)
cat(paste('Producing Figure 4b\n'))
pdf(file = "Figure_4b.pdf", width = 6, height = 2)
  par(mfrow = c(1, 3),
      las = 1)

  cat(paste('p-value indicating significant differences in plasticity between random and segregating variants in:\n'))
  ### calculate 2D and 3D plasticity for all variants
  TXN.P <- list()
  for (name in names(TXN)[1:4]) {
    group <- unlist(strsplit(name, ':'))[2]
    glu <- TXN[[paste0('_glucose_:', group)]]
    gal <- TXN[[paste0('_galactose_:', group)]]
    lac <- TXN[[paste0('_lactose_:', group)]]
    glu.gal <- c()
    glu.lac <- c()
    gal.lac <- c()
    all <- c()
    for (i in 1:length(glu)) {
      x1 <- rep(2, 2)
      x2 <- rep(6, 2)
      x0 <- c(glu[i], gal[i])
      glu.gal <- c(glu.gal, dist2d(x0, x1, x2))
      x0 <- c(glu[i], lac[i])
      glu.lac <- c(glu.lac, dist2d(x0, x1, x2))
      x0 <- c(gal[i], lac[i])
      gal.lac <- c(gal.lac, dist2d(x0, x1, x2))
      x1 <- rep(2, 3)
      x2 <- rep(6, 3)
      x0 <- c(glu[i], gal[i], lac[i])
      all <- c(all, dist3d(x0, x1, x2))
    }
    TXN.P[[paste0(names(cond.ls)[1], '-', names(cond.ls)[2], ':', group)]] <- glu.gal
    TXN.P[[paste0(names(cond.ls)[1], '-', names(cond.ls)[3], ':', group)]] <- glu.lac
    TXN.P[[paste0(names(cond.ls)[2], '-', names(cond.ls)[3], ':', group)]] <- gal.lac
    TXN.P[[paste0('all_', group)]] <- all
  }
  ### prepare dataset for plotting
  for (name in names(TXN.P)[1:3]) {
    if (name == names(TXN.P)[1]) {
      par(mar = c(3.5, 3.5, 1, 0.5))
    } else {
      if (name == names(TXN.P)[2]) {
        par(mar = c(3.5, 2.5, 1, 1.5))
      } else {
        par(mar = c(3.5, 1.5, 1, 2.5))
      }
    }
    group <- unlist(strsplit(name, ':'))[1]
    mut1 <- unlist(TXN.P[which(grepl(group, names(TXN.P)) & grepl('random_1', names(TXN.P)))])
    mut2 <- unlist(TXN.P[which(grepl(group, names(TXN.P)) & grepl('random_2', names(TXN.P)))])
    mut3 <- unlist(TXN.P[which(grepl(group, names(TXN.P)) & grepl('random_3', names(TXN.P)))])
    seg <- unlist(TXN.P[which(grepl(group, names(TXN.P)) & grepl('segregating', names(TXN.P)))])
    ### check whether increase and decrease in fluorescence
    ### as compared to MG1655 variant are equally probable
    mg1655 <- seg[length(seg)]
    m <- c(mut1, mut2, mut3)
    bin <- binom.test(x = length(m[which(m > mg1655)]), n = length(m))
    ### compare values from random and segregating variants
    wil <- wilcox.test(x = m, y = seg, exact = F)
    conds <- unlist(strsplit(group, '-'))
    cat(paste('\t', cond.ls[[conds[1]]], 'vs.', cond.ls[[conds[2]]], ':', signif(wil$p.value, digits = 3), '\n'))
    ### plotting
    plot(sort(mut1), (1:length(mut1) / length(mut1)),
          xlim = c(0, 2), ylim = c(0, 1), pch = 16,
          type = 'o', col = alpha('blue', 0.75), xlab = '',
          ylab = '', main = '', cex = 0.75)
    points(sort(mut2), (1:length(mut2) / length(mut2)), pch = 16,
          type = 'o', col = alpha('gold', 0.75), cex = 0.75)
    points(sort(mut3), (1:length(mut3) / length(mut3)), pch = 16,
          type = 'o', col = alpha('red', 0.75), cex = 0.75)
    points(sort(seg), (1:length(seg) / length(seg)), pch = 16,
          type = 'o', col = alpha('black', 0.75), cex = 0.75)
    points(seg[length(seg)], (which(sort(seg) == seg[length(seg)]) / length(seg)),
            bg = 'white', pch = 21, cex = 0.8)
    if (bin$p.value <= 0.05) {
      mtext(signif(bin$p.value, digits = 2), side = 3, line = -1, adj = 0.05, cex = 0.5, font = 4)
    } else {
      mtext(signif(bin$p.value, digits = 2), side = 3, line = -1, adj = 0.05, cex = 0.5)
    }
    mtext(paste0(short.ls[[conds[1]]], ':', short.ls[[conds[2]]]), side = 1, line = -1, adj = 0.95, cex = 0.75, font = 4)
    if (name == names(TXN.P)[1]) {
      title(ylab = 'Cumulative density', line = 2.5)
    }
    title(xlab = 'Plasticity', line = 2.25)
  }

dev.off()

### produce Figure 4c (comparing expression values of aceB
### promoter among all three environments -> 3D)
cat(paste('Producing Figure 4c\n'))
pdf(file = 'Figure_4c.pdf', width = 5, height = 5)
  par(las = 1)

  ### get fluorescence values for environment 1
  cond1.mut1 <- unlist(TXN[which(grepl('_glucose_', names(TXN)) & grepl('random_1', names(TXN)))])
  cond1.mut2 <- unlist(TXN[which(grepl('_glucose_', names(TXN)) & grepl('random_2', names(TXN)))])
  cond1.mut3 <- unlist(TXN[which(grepl('_glucose_', names(TXN)) & grepl('random_3', names(TXN)))])
  cond1.seg <- unlist(TXN[which(grepl('_glucose_', names(TXN)) & grepl('segregating', names(TXN)))])
  ### get fluorescence levels for environment 2
  cond2.mut1 <- unlist(TXN[which(grepl('_galactose_', names(TXN)) & grepl('random_1', names(TXN)))])
  cond2.mut2 <- unlist(TXN[which(grepl('_galactose_', names(TXN)) & grepl('random_2', names(TXN)))])
  cond2.mut3 <- unlist(TXN[which(grepl('_galactose_', names(TXN)) & grepl('random_3', names(TXN)))])
  cond2.seg <- unlist(TXN[which(grepl('_galactose_', names(TXN)) & grepl('segregating', names(TXN)))])
  ### get fluorescence levels for environment 3
  cond3.mut1 <- unlist(TXN[which(grepl('_lactose_', names(TXN)) & grepl('random_1', names(TXN)))])
  cond3.mut2 <- unlist(TXN[which(grepl('_lactose_', names(TXN)) & grepl('random_2', names(TXN)))])
  cond3.mut3 <- unlist(TXN[which(grepl('_lactose_', names(TXN)) & grepl('random_3', names(TXN)))])
  cond3.seg <- unlist(TXN[which(grepl('_lactose_', names(TXN)) & grepl('segregating', names(TXN)))])
  source(paste0(root.path, '/addgrids3d.r'))
  sp <- scatterplot3d(cond1.mut1, cond2.mut1, cond3.mut1,
        pch = '', grid = F, box = F,
        xlab = cond.ls[1], xlim = lims,
        zlab = cond.ls[2], ylim = lims,
        ylab = cond.ls[3], zlim = lims,
        angle = 55)
  addgrids3d(lims, lims, lims, grid = c('xy', 'xz', 'yz'), angle = 55)
  sp$points(cond1.mut1, cond2.mut1, cond3.mut1,
        col = alpha('blue', 0.75), pch = 16)
  sp$points3d(cond1.mut2, cond2.mut2, cond3.mut2,
        col = alpha('gold', 0.75), pch = 16)
  sp$points3d(cond1.mut3, cond2.mut3, cond3.mut3,
        col = alpha('red', 0.75), pch = 16)
  sp$points3d(cond1.seg, cond2.seg, cond3.seg,
        col = alpha('black', 0.75), pch = 16)
  sp$points3d(x = lims, y = lims, z = lims, type = 'l', lty = 3, col = 'blue')

dev.off()

### produce Figure 4d (cumulative 3D plasticity)
cat(paste('Producing Figure 4d\n'))
pdf(file = "Figure_4d.pdf", width = 6, height = 3)
  par(mar = c(3.5, 3.5, 1, 0.5),
      mfrow = c(1, 2),
      las = 1)

  cat(paste('p-value indicating significant differences in plasticity between random and segregating variants in:\n'))
  ### prepare dataset for plotting
  mut1 <- unlist(TXN.P[which(grepl('all_', names(TXN.P)) & grepl('random_1', names(TXN.P)))])
  mut2 <- unlist(TXN.P[which(grepl('all_', names(TXN.P)) & grepl('random_2', names(TXN.P)))])
  mut3 <- unlist(TXN.P[which(grepl('all_', names(TXN.P)) & grepl('random_3', names(TXN.P)))])
  seg <- unlist(TXN.P[which(grepl('all_', names(TXN.P)) & grepl('segregating', names(TXN.P)))])
  ### check whether increase and decrease in fluorescence
  ### as compared to MG1655 variant are equally probable
  mg1655 <- seg[length(seg)]
  m <- c(mut1, mut2, mut3)
  bin <- binom.test(x = length(m[which(m > mg1655)]), n = length(m))
  ### compare values from random and segregating variants
  wil <- wilcox.test(x = m, y = seg, exact = F)
  cat(paste('\t all three environments:', signif(wil$p.value, digits = 3), '\n'))
  ### plotting
  plot(sort(mut1), (1:length(mut1) / length(mut1)),
        xlim = c(0, 2), ylim = c(0, 1), pch = 16,
        type = 'o', col = alpha('blue', 0.75), xlab = '',
        ylab = '', main = '', cex = 0.75)
  points(sort(mut2), (1:length(mut2) / length(mut2)), pch = 16,
        type = 'o', col = alpha('gold', 0.75), cex = 0.75)
  points(sort(mut3), (1:length(mut3) / length(mut3)), pch = 16,
        type = 'o', col = alpha('red', 0.75), cex = 0.75)
  points(sort(seg), (1:length(seg) / length(seg)), pch = 16,
        type = 'o', col = alpha('black', 0.75), cex = 0.75)
  points(seg[length(seg)], (which(sort(seg) == seg[length(seg)]) / length(seg)),
          bg = 'white', pch = 21, cex = 0.8)
  if (bin$p.value <= 0.05) {
    mtext(signif(bin$p.value, digits = 2), side = 3, line = -1, adj = 0.05, cex = 0.75, font = 4)
  } else {
    mtext(signif(bin$p.value, digits = 2), side = 3, line = -1, adj = 0.05, cex = 0.75)
  }
  title(ylab = 'Cumulative density', line = 2.5)
  title(xlab = 'Plasticity', line = 2.25)

  par(mar = c(3.5, 1.5, 1, 3.5))
  plot(NULL, axes = F, ann = F, xlim = c(0, 1), ylim = c(0, 1))
  legend('left', legend = c('Segregating', '1 mutation', '2 mutations', '3 mutations'), pch = 16, col = cols, box.col = 'white', cex = 0.85)

dev.off()

##################################################
################## NOISE PLOTS ###################
############ FIGURE 5, 6c & SUPP FIG 5 ###########
##################################################

### define colors for ploting
cols <- c('magenta', 'cyan', 'gold')
names(cols) <- names(cond.ls)
noiTXN <- list()

### produce Figure 5a (fitting smooth spline
### to coefficient of variation and modal fluorescence data)
cat(paste('Producing Figure 5a\n'))
pdf(file = 'Figure_5a.pdf', width = 6, height = 3)
  par(mar = c(3.5, 4, 1, 1),
      mfrow = c(1, 2),
      las = 1)

  ### loop through all the environments
  for (cond in names(cond.ls)) {
    ### get and correct fluorescence level and
    ### coefficient of variation values for random variants
    mut.m <- mutTXN.ls[[cond]][, 3][which(as.numeric(mutTXN.ls[[cond]][, 2]) > 0 & as.numeric(mutTXN.ls[[cond]][, 2]) < 4)]
    mut.n <- mutTXN.ls[[cond]][, 4][which(as.numeric(mutTXN.ls[[cond]][, 2]) > 0 & as.numeric(mutTXN.ls[[cond]][, 2]) < 4)]
    mut.m <- as.numeric(mut.m) - offsetTXN.m[cond]
    mut.n <- as.numeric(mut.n) - offsetTXN.n[cond]
    mutNo <- mutTXN.ls[[cond]][, 2][which(as.numeric(mutTXN.ls[[cond]][, 2]) > 0 & as.numeric(mutTXN.ls[[cond]][, 2]) < 4)]
    ### get fluorescence level and coefficient of variation
    ### values for seg. variants in MG1655 background
    seg.m <- as.numeric(segTXN.ls[[cond]][, 2][which(grepl('SC', segTXN.ls[[cond]][, 1]))])
    seg.n <- as.numeric(segTXN.ls[[cond]][, 3][which(grepl('SC', segTXN.ls[[cond]][, 1]))])
    ### get fluorescence level and coefficient of variation
    ### values for seg. variants in native background
    nat.m <- as.numeric(natTXN.ls[[cond]][1:4, 2])
    nat.m <- c(nat.m, as.numeric(segTXN.ls[[cond]][, 2][which(grepl('MG1655.TXN', segTXN.ls[[cond]][, 1]))]))
    nat.n <- as.numeric(natTXN.ls[[cond]][1:4, 3])
    nat.n <- c(nat.n, as.numeric(segTXN.ls[[cond]][, 3][which(grepl('MG1655.TXN', segTXN.ls[[cond]][, 1]))]))
    ### calculate the smooth spline
    li <- smooth.spline(c(mut.m, seg.m, nat.m), c(mut.n, seg.n, nat.n), lambda = 0.01)
    ### predict coefficient of variation levels for each observed
    ### fluorescence level from the calculated smooth spline
    fit <- predict(li, c(mut.m, seg.m, nat.m))
    noise <- c(mut.n, seg.n, nat.n) - fit$y
    noise <- cbind(c(mutNo, rep(NA, length(seg.n) + length(nat.n))), noise)
    rownames(noise) <- c(mutTXN.ls[[cond]][, 1][which(as.numeric(mutTXN.ls[[cond]][, 2]) > 0 & as.numeric(mutTXN.ls[[cond]][, 2]) < 4)],
                segTXN.ls[[cond]][, 1][which(grepl('SC', segTXN.ls[[cond]][, 1]))],
                natTXN.ls[[cond]][1:4, 1], 'MG1655')
    noiTXN[[cond]] <- noise
    ### plotting
    if (cond == names(cond.ls)[1]) {
      plot(mut.m, mut.n, pch = 16, cex = 0.85,
            col = alpha(cols[cond], 0.3), xlim = c(2, 5.25),
            ylim = c(0.02, 0.2), log = 'y',
            xlab = '', ylab = '', main = '')
    } else {
      points(mut.m, mut.n, pch = 16, col = alpha(cols[cond], 0.3), cex = 0.85)
    }
    points(seg.m, seg.n, pch = 21, col = 'black',
          bg = alpha(cols[cond], 0.3), cex = 0.85)
    points(nat.m, nat.n, pch = 21, col = cols[cond],
          bg = alpha('black', 0.3), cex = 0.85)
    lines(predict(li), col = cols[cond], pch = '.')
  }
  title(xlab = 'Modal fluorescence (log10, a.u.)', line = 2)
  title(ylab = 'mCV (stdev/mode)', line = 3)
  par(mar = c(1, 1, 1, 1))
  plot(NULL, axes = F, ann = F, xlim = c(0, 1), ylim = c(0, 1))
  legend('top', legend = c('Random in MG1655 background', 'Segregating in MG1655 background', 'Segregating in native background'),
          pch = 21, title = 'Variant sets', col = c(alpha('magenta', 0.3), 'black', 'magenta'),
          pt.bg = c(alpha('magenta', 0.3), alpha('magenta', 0.3), alpha('black', 0.3)), box.col = 'white', cex = 0.85)
  legend('bottomleft', legend = unlist(cond.ls), pch = 16,
          col = c(alpha(cols[1], 0.3), alpha(cols[2], 0.3), alpha(cols[3], 0.3)),
          box.col = 'white', cex = 0.85)
  legend('bottomright', legend = c('Segregating', '1 mutation', '2 mutations', '3 mutations'),
          pch = 16, col = c('black', 'blue', 'gold', 'red'), box.col = 'white', cex = 0.85)

dev.off()

### produce Figure 5b (comparing noise
### between seg. and random variants)
cat(paste('Producing Figure 5b\n'))
pdf(file = 'Figure_5b.pdf', width = 6, height = 2)
  par(mfrow = c(1, 3),
      las = 1)

  cat(paste('p-value indicating significant differences in noise between random and segregating variants in:\n'))
  ### loop through all environments
  for (n in 1:3) {
    cond <- names(cond.ls)[n]
    if (n == 1) {
      par(mar = c(3.5, 3.5, 1, 0.5))
    } else {
      if (n == 2) {
        par(mar = c(3.5, 2.5, 1, 1.5))
      } else {
        par(mar = c(3.5, 1.5, 1, 2.5))
      }
    }
    ### get noise values for random variants
    mut1 <- as.numeric(noiTXN[[cond]][, 2][which(as.numeric(noiTXN[[cond]][, 1]) == 1)])
    mut2 <- as.numeric(noiTXN[[cond]][, 2][which(as.numeric(noiTXN[[cond]][, 1]) == 2)])
    mut3 <- as.numeric(noiTXN[[cond]][, 2][which(as.numeric(noiTXN[[cond]][, 1]) == 3)])
    ### get noise values for seg. variants in MG1655
    seg <- as.numeric(noiTXN[[cond]][, 2][which(is.na(noiTXN[[cond]][, 1]))])
    seg <- c(seg[1:(length(seg) - 5)], seg[length(seg)])
    ### check whether increase and decrease in fluorescence
    ### as compared to MG1655 variant are equally probable
    mg1655 <- seg[length(seg)]
    m <- c(mut1, mut2, mut3)
    bin <- binom.test(x = length(m[which(m > mg1655)]), n = length(m))
    ### compare values from random and segregating variants
    wil <- wilcox.test(x = m, y = seg, exact = F)
    cat(paste('\t', cond.ls[[cond]], ':', signif(wil$p.value, digits = 3),
              '\tmedian random =', signif(median(m), digits = 3), 'median segregating =', signif(median(seg), digits = 3), '\n'))
    ### plotting
    plot(sort(mut1), (1:length(mut1) / length(mut1)),
          pch = 16, type = 'o', cex = 0.75,
          xlim = c(-0.035, 0.025), ylim = c(0, 1),
          col = alpha('blue', 0.75),
          xlab = '', ylab = '', main = '')
    points(sort(mut2), (1:length(mut2) / length(mut2)),
            col = alpha('gold', 0.75), pch = 16, type = 'o', cex = 0.75)
    points(sort(mut3), (1:length(mut3) / length(mut3)),
            col = alpha('red', 0.75), pch = 16, type = 'o', cex = 0.75)
    points(sort(seg), (1:length(seg) / length(seg)),
            col = alpha('black', 0.75), pch = 16, type = 'o', cex = 0.75)
    points(seg[length(seg)], (which(sort(seg) == seg[length(seg)]) / length(seg)),
            bg = 'white', pch = 21, cex = 0.8)
    if (bin$p.value <= 0.05) {
      mtext(signif(bin$p.value, digits = 2), side = 3, line = -1, adj = 0.05, cex = 0.5, font = 4)
    } else {
      mtext(signif(bin$p.value, digits = 2), side = 3, line = -1, adj = 0.05, cex = 0.5)
    }
    mtext(cond.ls[[cond]], side = 1, line = -1, adj = 0.95, cex = 0.75, font = 4)
    if (cond == names(cond.ls)[1]) {
      title(ylab = 'Cumulative density', line = 2.5)
    }
    title(xlab = 'Noise', line = 2.25)
  }

dev.off()

### produce Figure 6b (comparing noise
### in seg. variants between environments)
cat(paste('Producing Figure 6b\n'))
pdf(file = "Figure_6b.pdf", width = 4, height = 8)
  par(mar = c(2, 4.5, 3, 1),
      mfrow = c(2, 1),
      las = 1)

  ### obtain noise values seg. variants in glucose
  cond <- names(cond.ls)[1]
  gl <- noiTXN[[cond]][, 2][which(is.na(noiTXN[[cond]][, 1]))]
  glu <- as.numeric(c(gl[1:(length(gl) - 5)], gl[length(gl)]))
  ### obtain noise values seg. variants in galactose
  cond <- names(cond.ls)[2]
  ga <- noiTXN[[cond]][, 2][which(is.na(noiTXN[[cond]][, 1]))]
  gal <-  as.numeric(c(ga[1:(length(ga) - 5)], ga[length(ga)]))
  ### obtain noise values seg. variants in lactose
  cond <- names(cond.ls)[3]
  la <- noiTXN[[cond]][, 2][which(is.na(noiTXN[[cond]][, 1]))]
  lac <-  as.numeric(c(la[1:(length(la) - 5)], la[length(la)]))
  ### test for correlations in noise between environments
  cor1 <- cor.test(glu, gal, method = 'spearman')
  cor2 <- cor.test(glu, lac, method = 'spearman')
  ### plotting
  plot(glu, gal, pch = 16, col = alpha('black', 0.2),
          xlim = c(-0.035, 0.015), ylim = c(-0.015, 0.01),
          xlab = '', ylab = '', yaxt = 'n')
  axis(side = 2, at = c(-0.01, 0, 0.01), labels = c(-0.01, '0.00', 0.01))
  title(ylab = 'Noise in Galactose', line = 3.5)
  fit <- lm(gal ~ glu)
  abline(fit, lwd = 0.5)
  if (cor1$p.value <= 0.05) {
    text(x = -0.01, y = 0.01, labels = paste('rho =', signif(cor1$estimate, digits = 2)), font = 4)
    text(x = -0.01, y = 0.008, labels = paste('p =', signif(cor1$p.value, digits = 2)), font = 4)
  } else {
    text(x = -0.01, y = 0.01, labels = paste('rho =', signif(cor1$estimate, digits = 2)))
    text(x = -0.01, y = 0.008, labels = paste('p =', signif(cor1$p.value, digits = 2)))
  }

  par(mar = c(4, 4.5, 1, 1))
  plot(glu, lac, pch = 16, col = alpha('black', 0.2),
          xlim = c(-0.035, 0.015), ylim = c(-0.015, 0.01),
          xlab = '', ylab = '', yaxt = 'n')
  axis(side = 2, at = c(-0.01, 0, 0.01), labels = c(-0.01, '0.00', 0.01))
  title(ylab = 'Noise in Lactose', line = 3.5)
  title(xlab = 'Noise in Glucose', line = 2)
  fit <- lm(lac ~ glu)
  abline(fit, lwd = 0.5)
  if (cor2$p.value <= 0.05) {
    text(x = -0.01, y = 0.01, labels = paste('rho =', signif(cor2$estimate, digits = 2)), font = 4)
    text(x = -0.01, y = 0.008, labels = paste('p =', signif(cor2$p.value, digits = 2)), font = 4)
  } else {
    text(x = -0.01, y = 0.01, labels = paste('rho =', signif(cor2$estimate, digits = 2)))
    text(x = -0.01, y = 0.008, labels = paste('p =', signif(cor2$p.value, digits = 2)))
  }

dev.off()

##################################################
############### NATIVE BACKGROUND ################
################### FIGURE 7 #####################
##################################################

setwd(root.path)
### sort content of mixed libraries into separate variables
segTXN.ls2 <- list()
natTXN.ls2 <- list()
for (path in paths.noRJ) {
  ### extract info about last directory
  end.dir <- unlist(strsplit(path, split = '/', fixed = T))[ndir]
  cond <- unlist(strsplit(end.dir, split = '_', fixed = T))[3]
  ### load modal fluorescence and coefficient of variation data
  data.m <- read.csv(paste0(path, 'AllPeaks.csv'), header = T)
  data.m <- data.m[, -1]
  data.n <- read.csv(paste0(path, 'Allcvs.csv'), header = T)
  data.n <- data.n[, -1]
  ### extract fluorescence data from segr. variants
  ### in MG1655 with RiboJ shared in chrosome & native bakcgrounds
  seg <- cbind(rbind(data.m[76,], data.m[85:86,], data.m[91,]),
              rbind(data.n[76,], data.n[85:86,], data.n[91,]))
  names <- c('SC312', 'SC400', 'SC418', 'MG1655')
  df.seg <- cbind(names, seg)
  colnames(df.seg) <- c('name', 'mode1', 'mode2', 'mode3', 'sd1', 'sd2', 'sd3')
  segTXN.ls2[[paste0('_', cond, '_')]] <- as.data.frame(df.seg)
  ### extract fluorescence data from segr. variants
  ### in native background with RiboJ including controls
  nat <- cbind(rbind(data.m[93:96,], data.m[91,]),
              rbind(data.n[93:96,], data.n[91,]))
  names <- c('SC312', 'SC358', 'SC400', 'SC418', 'MG1655')
  df.nat <- cbind(names, nat)
  colnames(df.nat) <- c('name', 'mode1', 'mode2', 'mode3', 'sd1', 'sd2', 'sd3')
  natTXN.ls2[[paste0('_', cond, '_')]] <- as.data.frame(df.nat)
}

### produce Figure 7a (modal fluorescence in native bbbackground)
cols <- c('black', 'red', 'blue', 'gold')
cat(paste('Producing Figure 7a\n'))
pdf(file = 'Figure_7a.pdf', width = 4, height = 3)
par(mar = c(3.5, 4, 1, 1),
    las = 1)

    ### loop through all environments
    for (n in 1:3) {
      cond <- names(cond.ls)[n]
      ### save fluorescence values from MG1655 background
      mg.all <- segTXN.ls2[[cond]][, 2:4]
      rownames(mg.all) <- segTXN.ls2[[cond]][, 1]
      mg <- as.numeric(segTXN.ls[[cond]][, 2][which(segTXN.ls[[cond]][, 1] %in% segTXN.ls2[[cond]][, 1])])
      mg <- c(mg, as.numeric(segTXN.ls[[cond]][, 2][which(grepl('MG1655.TXN', segTXN.ls[[cond]][, 1]))]))
      names(mg) <- c(segTXN.ls[[cond]][, 1][which(segTXN.ls[[cond]][, 1] %in% segTXN.ls2[[cond]][, 1])], 'MG1655')
      ### save fluorescence values from native backgrounds
      nat.all <- natTXN.ls2[[cond]][, 2:4]
      rownames(nat.all) <- natTXN.ls2[[cond]][, 1]
      nat <- as.numeric(c(natTXN.ls[[cond]][1:4, 2], natTXN.ls[[cond]][10, 2]))
      names(nat) <- c(natTXN.ls[[cond]][1:4, 1], 'MG1655')
      ### plotting
      for (i in 1:length(names(mg))) {
        name <- names(mg)[i]
        pl.line <- c(mg[name], nat[name])
        pl.point <- rbind(as.numeric(mg.all[name,]), as.numeric(nat.all[name,]))
        if (n == 1 && i == 1) {
          plot(rep(n - 0.25, 3), pl.point[1,], pch = 16, xaxt = 'n',
              ylim = c(2, 5), xlim = c(0.5, 3.5), main = '',
              col = alpha(cols[i], 0.5), xlab = '', ylab = '', cex = 0.75)
        } else {
          points(rep(n - 0.25, 3), pl.point[1,], pch = 16, col = alpha(cols[i], 0.5), cex = 0.75)
        }
        points(rep(n + 0.25, 3), pl.point[2,], pch = 16, col = alpha(cols[i], 0.5), cex = 0.75)
        lines(c(n - 0.25, n + 0.25), pl.line, col = alpha(cols[i], 0.5))
      }
    }
    text(x = c(0.75, 1.25, 1.75, 2.25, 2.75, 3.25), y = rep(2, 6), labels = rep(c('MG', 'NAT'), 3), cex = 0.5)
    axis(side = 1, at = c(1:3), labels = unlist(cond.ls))
    title(ylab = 'Modal fluorescence (log10, a.u.)', line = 2.5)

dev.off()

### produce Figure 7b (plasticity in native background)
cat(paste('Producing Figure 7b\n'))
pdf(file = 'Figure_7b.pdf', width = 4.5, height = 3)
par(mar = c(3.5, 4, 1, 1),
    las = 1)

    ### save fluorescence values from MG1655 background
    glu.mg <- as.numeric(segTXN.ls[['_glucose_']][, 2][which(segTXN.ls[['_glucose_']][, 1] %in% segTXN.ls2[['_glucose_']][, 1])])
    glu.mg <- c(glu.mg, as.numeric(segTXN.ls[['_glucose_']][, 2][which(grepl('MG1655.TXN', segTXN.ls[['_glucose_']][, 1]))]))
    gal.mg <- as.numeric(segTXN.ls[['_galactose_']][, 2][which(segTXN.ls[['_glucose_']][, 1] %in% segTXN.ls2[['_glucose_']][, 1])])
    gal.mg <- c(gal.mg, as.numeric(segTXN.ls[['_galactose_']][, 2][which(grepl('MG1655.TXN', segTXN.ls[['_glucose_']][, 1]))]))
    lac.mg <- as.numeric(segTXN.ls[['_lactose_']][, 2][which(segTXN.ls[['_glucose_']][, 1] %in% segTXN.ls2[['_glucose_']][, 1])])
    lac.mg <- c(lac.mg, as.numeric(segTXN.ls[['_lactose_']][, 2][which(grepl('MG1655.TXN', segTXN.ls[['_glucose_']][, 1]))]))
    mg <- cbind(glu.mg, gal.mg, lac.mg)
    rownames(mg) <- c(segTXN.ls[['_glucose_']][, 1][which(segTXN.ls[['_glucose_']][, 1] %in% segTXN.ls2[['_glucose_']][, 1])], 'MG1655')
    ### save fluorescence values from native backgrounds
    glu.nat <- as.numeric(c(natTXN.ls[['_glucose_']][1:4, 2], natTXN.ls[['_glucose_']][10, 2]))
    gal.nat <- as.numeric(c(natTXN.ls[['_galactose_']][1:4, 2], natTXN.ls[['_galactose_']][10, 2]))
    lac.nat <- as.numeric(c(natTXN.ls[['_lactose_']][1:4, 2], natTXN.ls[['_lactose_']][10, 2]))
    nat <- cbind(glu.nat, gal.nat, lac.nat)
    rownames(nat) <- c(natTXN.ls[['_glucose_']][1:4, 1], 'MG1655')
    ### save environment combinations for future reference in plasticity calculation
    combs <- rbind(c(1, 2), c(1, 3), c(2, 3))
    ### loop through all environmental combinations (including 3D)
    ### and calculate plasticity of each group separately
    for (clm in 1:4) {
      plast.mg <- c()
      plast.nat <- c()
      if (clm < 4) {
        x1 <- rep(0, 2)
        x2 <- rep(6, 2)
        for (i in 1:5) {
          if (i < 5) {
            x0 <- c(mg[i, combs[clm, 1]], mg[i, combs[clm, 2]])
            plast.mg <- c(plast.mg, dist2d(x0, x1, x2))
          }
          x0 <- c(nat[i, combs[clm, 1]], nat[i, combs[clm, 2]])
          plast.nat <- c(plast.nat, dist2d(x0, x1, x2))
        }
      } else {
        x1 <- rep(0, 3)
        x2 <- rep(6, 3)
        for (i in 1:5) {
          if (i < 5) {
            x0 <- c(mg[i, 1], mg[i, 2], mg[i, 3])
            plast.mg <- c(plast.mg, dist3d(x0, x1, x2))
          }
          x0 <- c(nat[i, 1], nat[i, 2], nat[i, 3])
          plast.nat <- c(plast.nat, dist3d(x0, x1, x2))
        }
      }
      ### add plasticity values to the matrix with fluorescence values
      mg <- cbind(mg, plast.mg)
      nat <- cbind(nat, plast.nat)
    }
    cnames <- c('Glu', 'Gal', 'Lac', 'Glu:Gal', 'Glu:Lac', 'Gal:Lac', 'All')
    colnames(mg) <- cnames
    colnames(nat) <- cnames
    ### plotting
    for (n in 1:4) {
      for (i in 1:length(rownames(mg))) {
        name <- rownames(mg)[i]
        ypl <- c(mg[name, n + 3], nat[name, n + 3])
        xpl <- c(n - 0.25, n + 0.25)
        if (n == 1 && i == 1) {
          plot(xpl, ypl, pch = 16, xaxt = 'n', type = 'o',
              ylim = c(-0.1, 1.6), xlim = c(0.5, 4.5), main = '',
              col = alpha(cols[i], 0.75), xlab = '', ylab = '', cex = 0.75)
        } else {
          points(xpl, ypl, pch = 16, type = 'o', col = alpha(cols[i], 0.75), cex = 0.75)
        }
      }
    }
    text(x = c(0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25), y = rep(-0.1, 8), labels = rep(c('MG', 'NAT'), 4), cex = 0.5)
    axis(side = 1, at = c(1:4), labels = colnames(mg)[4:7])
    title(ylab = 'Plasticity', line = 2.5)

dev.off()

### produce Figure 7c (noise native)
cat(paste('Producing Figure 7c\n'))
pdf(file = 'Figure_7c.pdf', width = 4, height = 3)
par(mar = c(3.5, 4.5, 1, 0.5),
    las = 1)

    ### loop through all environments
    for (n in 1:3) {
      cond <- names(cond.ls)[n]
      ### save noise for native backgrounds
      nat <- as.numeric(noiTXN[[cond]][110:114, 2])
      names(nat) <- names(noiTXN[[cond]][110:114, 2])
      ### save noise for MG1655 background
      mg <- as.numeric(noiTXN[[cond]][93:109, 2][which(names(noiTXN[[cond]][93:109, 2]) %in% names(nat))])
      mg <- c(mg, as.numeric(noiTXN[[cond]][114, 2]))
      names(mg) <- c(names(noiTXN[[cond]][93:109, 2][which(names(noiTXN[[cond]][93:109, 2]) %in% names(nat))]), 'MG1655')
      ### plotting
      for (i in 1:length(names(mg))) {
        name <- names(mg)[i]
        ypl <- c(mg[name], nat[name])
        xpl <- c(n - 0.25, n + 0.25)
        if (n == 1 && i == 1) {
          plot(xpl, ypl, pch = 16, xaxt = 'n',
              ylim = c(-0.04, 0.015), xlim = c(0.5, 3.5), main = '',
              col = alpha(cols[i], 0.75), xlab = '', ylab = '', cex = 0.75)
        } else {
          points(xpl, ypl, pch = 16, col = alpha(cols[i], 0.75), cex = 0.75)
        }
        lines(xpl, ypl, col = alpha(cols[i], 0.5))
      }
    }
    text(x = c(0.75, 1.25, 1.75, 2.25, 2.75, 3.25), y = rep(-0.04, 6), labels = rep(c('MG', 'NAT'), 3), cex = 0.5)
    axis(side = 1, at = c(1:3), labels = unlist(cond.ls))
    title(ylab = 'Noise', line = 3.5)

dev.off()

### produce Figure 7d (legend for Figure 7)
cat(paste('Producing Figure 7d\n'))
pdf(file = 'Figure_7d.pdf', width = 3, height = 3)
par(mar = c(3.5, 4, 1, 1),
    las = 1)

    plot(NULL, axes = F, ann = F, xlim = c(0, 1), ylim = c(0, 1))
    legend('left', legend = names(mg), pch = 16,
            col = c(alpha(cols, 0.75)), title = 'Variants')

dev.off()

##################################################
################## CHROMOSOME ####################
################# SUPP FIGURE 4 ##################
##################################################

### extract fluorescence data from all replicates
### of seg. variants in MG1655 chromosome
chrom.ls2 <- list()
data.m <- read.csv(paste0(path.chrom, 'AllPeaks.csv'), header = T)
data.m <- data.m[, -1]
data.n <- read.csv(paste0(path.chrom, 'Allcvs.csv'), header = T)
data.n <- data.n[, -1]
glu <- cbind(rbind(data.m[1,], data.m[4,], data.m[7,], data.m[10,], data.m[13,]),
            rbind(data.n[1,], data.n[4,], data.n[7,], data.n[10,], data.n[13,]))
colnames(glu) <- c('mode1', 'mode2', 'mode3', 'sd1', 'sd2', 'sd3')
gal <- cbind(rbind(data.m[2,], data.m[5,], data.m[8,], data.m[11,], data.m[14,]),
            rbind(data.n[2,], data.n[5,], data.n[8,], data.n[11,], data.n[14,]))
colnames(gal) <- c('mode1', 'mode2', 'mode3', 'sd1', 'sd2', 'sd3')
lac <- cbind(rbind(data.m[3,], data.m[6,], data.m[9,], data.m[12,], data.m[15,]),
            rbind(data.n[3,], data.n[6,], data.n[9,], data.n[12,], data.n[15,]))
colnames(lac) <- c('mode1', 'mode2', 'mode3', 'sd1', 'sd2', 'sd3')
names <- c('MG1655', 'SC312', 'SC400', 'SC418', 'SC358')
chrom.ls2[[paste0('_glucose_')]] <- as.data.frame(cbind(names, glu))
chrom.ls2[[paste0('_galactose_')]] <- as.data.frame(cbind(names, gal))
chrom.ls2[[paste0('_lactose_')]] <- as.data.frame(cbind(names, lac))

### produce Supplementary Figure 4a (modal fluorescence in chromosome)
cat(paste('Producing Supplementary Figure 4a\n'))
pdf(file = 'SupplementaryFigure_4a.pdf', width = 4, height = 3)
par(mar = c(3.5, 4, 1, 1),
    las = 1)

    ### loop through all environments
    for (n in 1:3) {
      cond <- names(cond.ls)[n]
      ### save fluorescence values from plasmid (in MG1655)
      mg.all <- segTXN.ls2[[cond]][, 2:4]
      rownames(mg.all) <- segTXN.ls2[[cond]][, 1]
      mg <- as.numeric(segTXN.ls[[cond]][, 2][which(segTXN.ls[[cond]][, 1] %in% segTXN.ls2[[cond]][, 1])])
      mg <- c(mg, as.numeric(segTXN.ls[[cond]][, 2][which(grepl('MG1655.TXN', segTXN.ls[[cond]][, 1]))]))
      names(mg) <- c(segTXN.ls[[cond]][, 1][which(segTXN.ls[[cond]][, 1] %in% segTXN.ls2[[cond]][, 1])], 'MG1655')
      ### save fluorescence values from chromosome
      ch.all <- chrom.ls2[[cond]][, 2:4]
      ch.all <- ch.all - offsetCH.m[cond]
      rownames(ch.all) <- chrom.ls2[[cond]][, 1]
      ch <- as.numeric(chrom.ls[[cond]][1:5, 2])
      ch <- ch - offsetCH.m[cond]
      names(ch) <- c(chrom.ls[[cond]][1:5, 1])
      ### plotting
      for (i in 1:length(names(mg))) {
        name <- names(mg)[i]
        ypl <- c(mg[name], ch[name])
        xpl <- c(n - 0.25, n + 0.25)
        if (n == 1 && i == 1) {
          plot(xpl, ypl, xaxt = 'n', type = 'l',
              ylim = c(2, 5), xlim = c(0.5, 3.5), main = '',
              col = alpha(cols[i], 0.75), xlab = '', ylab = '', cex = 0.75)
        } else {
          lines(xpl, ypl, col = alpha(cols[i], 0.75), cex = 0.75)
        }
        points(rep(n - 0.25, 3), mg.all[name, ], pch = 16, col = alpha(cols[i], 0.5), cex = 0.75)
        points(rep(n + 0.25, 3), ch.all[name, ], pch = 16, col = alpha(cols[i], 0.5), cex = 0.75)
      }
    }
    text(x = c(0.75, 1.25, 1.75, 2.25, 2.75, 3.25), y = rep(2, 6), labels = rep(c('MG', 'CH'), 3), cex = 0.5)
    axis(side = 1, at = c(1:3), labels = unlist(cond.ls))
    title(ylab = 'Modal fluorescence (log10, a.u.)', line = 2.5)

dev.off()

### produce Supplementary Figure 4b (relative modal fluorescence in chromosome)
cat(paste('Producing Supplementary Figure 4b\n'))
pdf(file = 'SupplementaryFigure_4b.pdf', width = 4, height = 3)
par(mar = c(3.5, 4, 1, 1),
    las = 1)

    ### loop through all environments
    for (n in 1:3) {
      cond <- names(cond.ls)[n]
      ### save fluorescence values from plasmid (in MG1655)
      mg <- as.numeric(segTXN.ls[[cond]][, 2][which(segTXN.ls[[cond]][, 1] %in% segTXN.ls2[[cond]][, 1])])
      mg <- c(mg, as.numeric(segTXN.ls[[cond]][, 2][which(grepl('MG1655.TXN', segTXN.ls[[cond]][, 1]))]))
      mg <- mg / mg[length(mg)]
      names(mg) <- c(segTXN.ls[[cond]][, 1][which(segTXN.ls[[cond]][, 1] %in% segTXN.ls2[[cond]][, 1])], 'MG1655')
      ### save fluorescence values from chromosome
      ch <- as.numeric(chrom.ls[[cond]][1:5, 2])
      ch <- ch - offsetCH.m[cond]
      ch <- ch / ch[1]
      names(ch) <- c(chrom.ls[[cond]][1:5, 1])
      ### plotting
      for (i in 1:length(names(mg))) {
        name <- names(mg)[i]
        ypl <- c(mg[name], ch[name])
        xpl <- c(n - 0.25, n + 0.25)
        if (n == 1 && i == 1) {
          plot(xpl, ypl, pch = 16, xaxt = 'n', type = 'o',
              ylim = c(0.9, 1.3), xlim = c(0.5, 3.5), main = '',
              col = alpha(cols[i], 0.75), xlab = '', ylab = '', cex = 0.75)
        } else {
          points(xpl, ypl, pch = 16, type = 'o', col = alpha(cols[i], 0.75), cex = 0.75)
        }
        abline(h = 1, lty = 3, lwd = 0.5, col = 'grey')
      }
    }
    text(x = c(0.75, 1.25, 1.75, 2.25, 2.75, 3.25), y = rep(0.9, 6), labels = rep(c('MG', 'CH'), 3), cex = 0.5)
    axis(side = 1, at = c(1:3), labels = unlist(cond.ls))
    title(ylab = 'Relative modal fluorescence', line = 2.5)

dev.off()

### produce Supplementary Figure 4c (plasticity in native background)
cat(paste('Producing Supplementary Figure 4c\n'))
pdf(file = 'SupplementaryFigure_4c.pdf', width = 4.5, height = 3)
par(mar = c(3.5, 4, 1, 1),
    las = 1)

    ### save fluorescence values from MG1655 background
    mg <- cbind(glu.mg, gal.mg, lac.mg)
    rownames(mg) <- c(segTXN.ls[['_glucose_']][, 1][which(segTXN.ls[['_glucose_']][, 1] %in% segTXN.ls2[['_glucose_']][, 1])], 'MG1655')
    ### save fluorescence values from native backgrounds
    glu.ch <- as.numeric(chrom.ls[['_glucose_']][1:5, 2])
    gal.ch <- as.numeric(chrom.ls[['_galactose_']][1:5, 2])
    lac.ch <- as.numeric(chrom.ls[['_lactose_']][1:5, 2])
    ch <- cbind(glu.ch, gal.ch, lac.ch)
    rownames(ch) <- c(chrom.ls[['_glucose_']][1:5, 1])
    ### save environment combinations for future reference in plasticity calculation
    combs <- rbind(c(1, 2), c(1, 3), c(2, 3))
    ### loop through all environmental combinations (including 3D)
    ### and calculate plasticity of each group separately
    for (clm in 1:4) {
      plast.mg <- c()
      plast.ch <- c()
      if (clm < 4) {
        x1 <- rep(0, 2)
        x2 <- rep(6, 2)
        for (i in 1:5) {
          if (i < 5) {
            x0 <- c(mg[i, combs[clm, 1]], mg[i, combs[clm, 2]])
            plast.mg <- c(plast.mg, dist2d(x0, x1, x2))
          }
          x0 <- c(ch[i, combs[clm, 1]], ch[i, combs[clm, 2]])
          plast.ch <- c(plast.ch, dist2d(x0, x1, x2))
        }
      } else {
        x1 <- rep(0, 3)
        x2 <- rep(6, 3)
        for (i in 1:5) {
          if (i < 5) {
            x0 <- c(mg[i, 1], mg[i, 2], mg[i, 3])
            plast.mg <- c(plast.mg, dist3d(x0, x1, x2))
          }
          x0 <- c(ch[i, 1], ch[i, 2], ch[i, 3])
          plast.ch <- c(plast.ch, dist3d(x0, x1, x2))
        }
      }
      ### add plasticity values to the matrix with fluorescence values
      mg <- cbind(mg, plast.mg)
      ch <- cbind(ch, plast.ch)
    }
    cnames <- c('Glu', 'Gal', 'Lac', 'Glu:Gal', 'Glu:Lac', 'Gal:Lac', 'All')
    colnames(mg) <- cnames
    colnames(ch) <- cnames
    ### plotting
    for (n in 1:4) {
      for (i in 1:length(rownames(mg))) {
        name <- rownames(mg)[i]
        ypl <- c(mg[name, n + 3], ch[name, n + 3])
        xpl <- c(n - 0.25, n + 0.25)
        if (n == 1 && i == 1) {
          plot(xpl, ypl, pch = 16, xaxt = 'n', type = 'o',
              ylim = c(-0.1, 1.6), xlim = c(0.5, 4.5), main = '',
              col = alpha(cols[i], 0.75), xlab = '', ylab = '', cex = 0.75)
        } else {
          points(xpl, ypl, pch = 16, type = 'o', col = alpha(cols[i], 0.75), cex = 0.75)
        }
      }
    }
    text(x = c(0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25), y = rep(-0.1, 8), labels = rep(c('MG', 'CH'), 4), cex = 0.5)
    axis(side = 1, at = c(1:4), labels = colnames(mg)[4:7])
    title(ylab = 'Plasticity', line = 2.5)

dev.off()
