library(AlphaSimR)
library(dplyr)
library(MCMCpack)
library(tictoc)
library(Matrix)
#library(snpStats)
'InvLogit'  = function(x){
    exp(x)/(1 + exp(x))
}
'LogitNoise' = function(x, sigma2){
    x   = ifelse(x > 0.99, 0.99, ifelse(x < 0.01, 0.01, x))
    y   = log(x/(1 - x))
    tmp = rnorm(length(x))*sqrt(sigma2) + y
    z   = ifelse(tmp > 4.0, 4.0 - (tmp - 4.0),
                 ifelse(tmp < -4.0, -4.0 + abs(tmp + 4.0), tmp))
    return(InvLogit(z))
}
'BuildMME'  = function(npool, nsires, sigma2e, sigma2a, G, y, nvec){
    #sigma2e       = 1 - h2
    #sigma2a       = h2
    X             = as.matrix(rep(1, npool))
    Z             = cbind(matrix(0, ncol = nsires, nrow = npool), diag(npool))
    W             = as(cbind(X, Z), 'dgCMatrix')#matrix(0, ncol = nsires + npool + 1, nrow = npool + nsires + 1)
    Rinv          = diag(nvec)/sigma2e
    LHS           = crossprod(W, crossprod(Rinv, W))
    rhs           = crossprod(W, crossprod(Rinv, y))
    #k             = sigma2e/sigma2a
    Ginv          = MASS::ginv(G)/sigma2a
    idx           = 2:(npool + nsires + 1)
    LHS[idx, idx] = LHS[idx, idx] + Ginv
    return(list(LHS = LHS, rhs = rhs))
}
'BuildMME_i' = function(nind, nsires, sigma2e, sigma2a, y, G){
    X           = as.matrix(rep(1, nind))
    Z           = cbind(matrix(0, nrow = nind, ncol = nsires), Diagonal(nind))
    Rinv        = Diagonal(nind)/sigma2e
    W           = rbind(cbind(crossprod(X, crossprod(Rinv, X)), crossprod(X, crossprod(Rinv, Z))),
                        cbind(crossprod(Z, crossprod(Rinv, X)), crossprod(Z, crossprod(Rinv, Z))))
    #W           = cbind(X, Z)
    #C           = as(crossprod(W), 'dgCMatrix')
    C           = as(W, 'dgCMatrix')
    idx         = 2:(nsires + nind + 1)
    C[idx, idx] = C[idx, idx] + MASS::ginv(G)*(1/sigma2a)
    rhs         = rbind(crossprod(X, crossprod(Rinv, y)), crossprod(Z, crossprod(Rinv, y)))
    return(list(LHS = C, rhs = rhs))
}
'PEVpool'   = function(npool, nsires, sigma2e, sigma2a, G, n_k){
    ZRZ         = as(matrix(0, ncol = npool + nsires, nrow = npool + nsires), 'dgCMatrix')
    diag(ZRZ)   = c(rep(0, nsires), n_k)
    Ginv        = MASS::ginv(G)
    C           = solve(ZRZ + sigma2e/sigma2a*Ginv)
    return(C)
}
'CorRel'    = function(npool, nsires, G, sigma2e, sigma2a, n_k){
    Z   = as(cbind(matrix(0, ncol = nsires, nrow = npool), diag(npool)), 'dgCMatrix')
    Va  = (G*sigma2a)%*%t(Z)%*%solve(Z%*%(G*sigma2a)%*%t(Z) + diag(n_k)*sigma2e, Z%*%(G*sigma2a))
    return(Va)
}
'buildPed'  = function (ID, Par1, Par2, gener = NULL, sex = NULL, add.ancestors = FALSE, unknown = 0) {
    ID <- as.character(ID)
    Par1 <- as.character(Par1)
    Par2 <- as.character(Par2)
    n <- length(ID)
    if (length(unique(ID)) != n) 
        warning("ID is not unique, removing duplicated individuals")
    if (length(Par1) != n) 
        stop("Par1 must have same length as ID")
    if (length(Par2) != n) 
        stop("Par2 must have same length as ID")
    if (!(0 %in% ID | "0" %in% ID | any(unknown %in% ID))) {
        Par1[Par1 == "0" | Par1 %in% unknown | Par1 == 
                 0] <- NA
        Par2[Par2 == "0" | Par2 %in% unknown | Par2 == 
                 0] <- NA
    }
    else if (as.character(ID[ID == 0 | ID == "0"]) == as.character(Par1[ID == 
                                                                        0 | ID == "0"]) & as.character(ID[ID == 0 | ID == 
                                                                                                          "0"]) == as.character(Par2[ID == 0 | ID == "0"])) {
        Par1 <- Par1[ID != 0 | ID != "0"]
        Par2 <- Par2[ID != 0 | ID != "0"]
        ID <- ID[ID != 0 | ID != "0"]
        Par1[Par1 == "0" | Par1 %in% unknown | Par1 == 
                 0] <- NA
        Par2[Par2 == "0" | Par2 %in% unknown | Par2 == 
                 0] <- NA
    }
    Pars <- unique(c(Par1[!is.na(Par1)], Par2[!is.na(Par2)]))
    ancestors <- Pars[!Pars %in% ID]
    ID <- c(ancestors, ID)
    Par1 <- c(rep(NA, length(ancestors)), Par1)
    Par2 <- c(rep(NA, length(ancestors)), Par2)
    if (!is.null(gener)) 
        gener <- c(rep(NA, length(ancestors)), gener)
    if (!is.null(sex)) 
        sex <- c(rep(NA, length(ancestors)), sex)
    n <- length(ID)
    if (is.null(gener)) {
        generOld <- gener <- rep(n + 100, n)
        gener[is.na(Par1) & is.na(Par2)] <- 0
        i <- 0
        while (!all(generOld == gener)) {
            generOld <- gener
            gener[Par1 %in% ID[gener == i]] <- i + 1
            gener[Par2 %in% ID[gener == i]] <- i + 1
            i <- i + 1
            if (i > n + 10) 
                break
        }
    }
    if (add.ancestors) 
        ancestors <- FALSE
    if (!is.null(sex)) 
        pedigree <- data.frame(ID = ID[!ID %in% ancestors], Par1 = Par1[!ID %in% 
                                                                            ancestors], Par2 = Par2[!ID %in% ancestors], gener = gener[!ID %in% 
                                                                                                                                           ancestors], sex = sex[!ID %in% ancestors], stringsAsFactors = FALSE)
    else pedigree <- data.frame(ID = ID[!ID %in% ancestors], 
                                Par1 = Par1[!ID %in% ancestors], Par2 = Par2[!ID %in% 
                                                                                 ancestors], gener = gener[!ID %in% ancestors], stringsAsFactors = FALSE)
    pedigree <- pedigree[!duplicated(pedigree), ]
    pedigree <- pedigree[order(pedigree$gener, partial = pedigree$ID), 
    ]
    class(pedigree) <- c("pedigree", "data.frame")
    pedigree[is.na(pedigree)] <- 0
    return(pedigree)
}
'renumped'  = function(Pedigree){
    ID              = match(Pedigree[, 1], Pedigree[, 1])
    SID             = match(Pedigree[, 2], Pedigree[, 1])
    SID[is.na(SID)] = 0
    DID             = match(Pedigree[, 3], Pedigree[, 1])
    DID[is.na(DID)] = 0
    ped             = data.frame(ID = ID, SID = SID, DID = DID, Old = Pedigree[, 1], stringsAsFactors = F)
    return(ped)
}
'PBLUP' = function(X, Z, Ainv, Y, alpha){
    xtx     = crossprod(X)
    xtz     = crossprod(X, Z)
    ztx     = crossprod(Z, X)
    ztzainv = crossprod(Z) + alpha * Ainv
    nr      = nrow(xtx) + nrow(ztx)
    LHS     = Matrix(0, ncol = nr, nrow = nr)
    LHS     = Matrix(0, ncol = nr, nrow = nr)
    LHS[1:nrow(xtx), 1:ncol(xtx)]                   = xtx
    LHS[nrow(xtx) + seq(1:nrow(ztx)), 1:ncol(ztx)]  = ztx
    LHS[1:nrow(xtz), ncol(xtx) + seq(1, ncol(xtz))] = xtz
    LHS[nrow(xtz) + seq(1:nrow(ztzainv)), ncol(ztx) + seq(1, ncol(ztzainv))] = ztzainv
    RHS                 = Matrix(0, ncol = 1, nrow = nr)
    xty                 = crossprod(X, Y)
    RHS[1:nrow(xty), 1] = xty
    zty                 = crossprod(Z, Y)
    RHS[nrow(xty) + seq(1, nrow(zty))] = zty
    sol                 = solve(LHS, RHS)
    return(list(sol = sol, C = LHS))
}

Rcpp::sourceCpp('/work/rmlewis/napov/PoolingSimulation/G3/XXtEig.cpp')
ngen                          = 15
nchr                          = 26
sigmavec                      = c(0.125, 0.625, 1.125)
poolsize                      = c(5, 10, 25, 50, 100)
NS                            = length(sigmavec)
NP                            = length(poolsize)
nrep                          = 100 ### make sure this is set to 100
ngen                          = 15
chrlen                        = c(275406953, 248966461, 223996068, 119216639, 
                                  107836144, 116888256, 100009711, 90615088,
                                  94583238, 86377204, 62170480, 79028859, 
                                  83079144, 62568341, 80783214, 71693149, 
                                  72251135, 68494538, 60445663, 51176841, 
                                  49987992, 50780147, 62282865, 42034648, 
                                  45223504, 44047080)
nmarkers                      = round(chrlen*2.2/1e5)
nqtl                          = round(chrlen*0.1/1e5)
nsnp                          = sum(nmarkers)
cat('Accuracy', 'Replication', 'Poolsize_Pools', 'Error', 'Design', 'Scaling', 'Sires_Pool', 'Prog_sire',
    'Sire_rep', 'MeanRel', 'SDRel', 'ConstError' , sep = " ", 
    file = '/work/rmlewis/napov/PoolingSimulation/G3/SelSim1.txt', append = T)
cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/SelSim1.txt', append = T)
cat('Mean', 'SD', 'Min', 'Max', 'PoolSize1', 'Npooled', 'Replication', 'Pools', 'Error' , sep = " ", 
    file = '/work/rmlewis/napov/PoolingSimulation/G3/SelClus1.txt', append = T)
cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/SelClus1.txt', append = T)
#cat('Scenario', 'Replication', 'Mean', 'SD', 'Max', 'Min',
#    file = '/work/rmlewis/napov/PoolingSimulation/Selection/LDsummary2.txt', append = T)
#cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/Selection/LDsummary2.txt', append = T)
cat('PoolSize', 'Error', 'Design', 'Scaling', 'ConstError', 'Replication', 'PopTBV', 
    'T10', 'T20', 'T40', 'B40', 'B20', 'B10', sep = ' ', 
    file = '/work/rmlewis/napov/PoolingSimulation/G3/TopBotSel.txt', append = T)
cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/TopBotSel.txt', append = T)
t1 = proc.time()
for(irep in 1:nrep){
    set.seed(202006 + irep)
    tic()
    founderPop                    = runMacs2(nInd = 1000, nChr = nchr, Ne = 350, 
                                             segSites = (nmarkers + nqtl)*1.1, bp = chrlen)
    toc()
    SP                            = SimParam$new(founderPop)
    SP$restrSegSites(minQtlPerChr = nqtl,
                     minSnpPerChr = nmarkers,
                     overlap = FALSE) 
    SP$addTraitA(nqtl) 
    SP$setVarE(h2 = 0.3) 
    SP$addSnpChip(nmarkers) 
    SP$setSexes("yes_sys")
    SP$setTrackPed(TRUE)
    nfemales                      = 500
    lsize                         = 2
    n                             = ceiling(nfemales*lsize)
    repfemales                    = 0.90
    pop                           = newPop(founderPop) #Create initial population
    rowid                         = 1:n
    i                             = 1
    while(nfemales < 5000){
        if(i == 1){
            selsire = pop[pop@sex == 'M']
            ids     = sort(sample(length(selsire@id), 25))
            sirepop = selsire[ids]
            seldam  = pop[pop@sex == 'F']
            idd     = 1:length(seldam@id)
            dampop  = seldam[idd]
        } else if (i > 1){
            selsire1 = which(pop1@sex == 'M')
            selsire2 = which(pop2@sex == 'M')
            ids1     = sample(selsire1, floor(nfemales*0.02))
            ids2     = sample(selsire2, floor(nfemales*0.03))
            sirepop  = c(pop1[ids1], pop2[ids2])
            seldam1  = which(pop1@sex == 'F')
            seldam2  = which(pop2@sex == 'F')
            idd1     = sample(seldam1, floor(nfemales*0.8))#seldam1
            idd2     = sample(seldam2, floor(nfemales*0.3))
            dampop   = c(pop1[idd1], pop2[idd2])
        }
        pop1        = c(dampop, sirepop)
        pop2        = randCross(pop1, nCrosses = sum(pop1@sex == 'F'), nProgeny = 2)
        nfemales    = sum(pop2@sex == 'F')
        cat('Increasing population size. Generation: ', i, '\n')
        i = i + 1
    }
    for(gen in 1:ngen){
        selsire1 = pop1[which(pop1@sex == 'M')]
        selsire2 = pop2[which(pop2@sex == 'M')]
        dfsire1  = data.frame(ID = selsire1@id, Pheno = selsire1@pheno)
        dfsire2  = data.frame(ID = selsire2@id, Pheno = selsire2@pheno)
        ids1     = head(arrange(dfsire1, -Pheno), 100)
        ids2     = head(arrange(dfsire2, -Pheno), 100)
        sirepop  = c(selsire1[which(selsire1@id%in%ids1$ID)], selsire2[which(selsire2@id%in%ids2$ID)])
        seldam1  = pop1[which(pop1@sex == 'F')]
        seldam2  = pop2[which(pop2@sex == 'F')]
        dfdam1   = data.frame(ID = seldam1@id, Pheno = seldam1@pheno)
        dfdam2   = data.frame(ID = seldam2@id, Pheno = seldam2@pheno)
        idd1     = head(arrange(dfdam1, -Pheno), floor(nfemales*0.8))
        idd2     = head(arrange(dfdam2, -Pheno), floor(nfemales*0.2))
        dampop   = c(seldam1[which(seldam1@id%in%idd1$ID)], seldam2[which(seldam2@id%in%idd2$ID)])
        pop1     = c(dampop, sirepop)
        pop2     = selectCross(pop1, nFemale = length(dampop@id), nMale = length(sirepop@id), nCrosses = nfemales,
                               use = 'pheno', selectTop = T, balance = T, nProgeny = 2)
        nfemales = sum(pop2@sex == 'F')        
        if(gen == (ngen - 1)){
            SireGeno = pullSnpGeno(c(pop1, pop2))
            SireTBV  = data.frame(ID = c(pop1@id, pop2@id), TBV = c(pop1@gv, pop2@gv))
        }
        if(gen == ngen){
            Genotypes = pullSnpGeno(pop2)
            ProgTBV   = data.frame(ID = pop2@id, TBV = pop2@gv)
        }
        cat('Under selection: ', gen, '\n')
    }
    cat('\n')
    DF          = data.frame(ID = pop2@id, Sire = pop2@father, Dam = pop2@mother, Sex = pop2@sex,
                             TBV = pop2@gv, Pheno = pop2@pheno, stringsAsFactors = F)
    DF          = DF %>% dplyr::mutate(ID = as.numeric(ID)) %>% dplyr::mutate(Sire = as.numeric(Sire)) %>%
        dplyr::mutate(Dam = as.numeric(Dam))
    nsires      = 30
    #LDmat       = matrix(0, nchr, 1)
    #c2          = cumsum(nmarkers)
    #last        = 0
    #for(c1 in 1:nchr){
    #    first       = last + 1
    #    last        = c2[c1]
    #    X           = Genotypes[, first:last]
    #    Y           = as(X, 'SnpMatrix')
    #    A           = ld(x = Y, stats = 'R.squared', depth = 1)
    #    LDmat[c1]   = mean(A@x, na.rm = T)
    #}
    #--------------#
    # Sample sires #
    #--------------#
    rsires      = sample(unique(DF$Sire), nsires)
    if(sum(rsires%in%rownames(SireGeno)) != 30){
        cat('Not all sires in Genotype file. Check population!\n')
        stop()
    }
    #-----------------#
    # Extreme Pooling #
    #-----------------#
    np          = dim(DF)[1]
    subid       = sort(sample(nrow(DF), 3000))
    nindpool    = length(subid)
    DFF         = arrange(DF, desc(Pheno))
    SDF         = c(head(DFF$ID, nindpool/2), tail(DFF$ID, nindpool/2))
    for(ps in 1:NP){
        for(ss in 1:NS){
            K                   = length(subid)/poolsize[ps]
            W                   = matrix(0, nrow = K, ncol = nsnp)
            n_k                 = rep(floor(nindpool/K), K)
            y_k                 = c()
            idq                 = 1:n_k[1]
            siresperpool        = c()
            progpersire         = c()
            SireRep             = array(0, dim = c(nsires, K))
            for(k in 1:K){
                pivec           = c(rdirichlet(1, rep(10, n_k[k])))
                idvec           = SDF[idq]
                siresperpool[k] = sum(rsires%in%DFF$Sire[which(DFF$ID%in%idvec)])
                progpersire[k]  = sum(DFF$ID[DFF$Sire%in%rsires[rsires%in%DFF$Sire[which(DFF$ID%in%idvec)]]]%in%idvec)
                SireRep[, k]    = (rsires%in%DFF$Sire[which(DFF$ID%in%idvec)])*1
                wtmp            = Genotypes[which(DF$ID%in%idvec), ]
                W[k, ]          = LogitNoise(0.5*colSums(pivec*wtmp), sigmavec[ss]^2)
                y_k[k]          = mean(DFF$Pheno[which(DFF$ID%in%idvec)])
                idq             = idq + rep(1, n_k[1])*n_k[1]
                cat('Pool: ', k, ' constructed\n')
            }
            cat("\n")
            idx                 = which(rownames(SireGeno)%in%rsires)
            X                   = SireGeno[idx, ]
            phat                = colMeans(X)/2
            muhat               = colMeans(W)
            indx                = which(phat < 0.01 | phat > 0.99)
            nrem                = nsnp - length(indx)
            M                   = rbind(X[, -indx]/2, W[, -indx])
            muvec               = colMeans(M)
            G                   = GetGPool(M, muvec, dim(M)[1], dim(M)[2])
            diag(G)             = diag(G)*1.01
            tic()
            MME                 = BuildMME(npool = K, nsires = nsires, sigma2e = 0.7, sigma2a = 0.3, 
                                           G = G, y_k, rep(n_k[k], K))
            ans                 = solve(MME$LHS, MME$rhs)
            C22                 = diag(ginv(as.matrix(MME$LHS)))[2:(nsires + 1)]
            rel                 = sqrt(pmax(1 - C22*(0.7/0.3), 0)) 
            toc()
            c1                  = cor(SireTBV$TBV[which(SireTBV$ID %in% rsires)], ans[2:(nsires + 1)])
            cat(c1, irep, poolsize[ps], sigmavec[ss], 'Extremes', 'Cov',
                mean(siresperpool), mean(progpersire), mean(rowSums(SireRep)), mean(rel), sd(rel), 'Unequal',
                sep = " ", file = '/work/rmlewis/napov/PoolingSimulation/G3/SelSim1.txt', append = T)
            cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/SelSim1.txt', append = T)
            d1                  = data.frame(TBV = SireTBV$TBV[which(SireTBV$ID %in% rsires)], 
                                             EBV = ans[2:(nsires + 1)])
            topvals             = rep(0, 3)
            topidx              = list(19:30, 25:30, 28:30)
            botvals             = rep(0, 3)
            botidx              = list(1:12, 1:6, 1:3) 
            for(i in 1:3){
                topvals[i] = d1 %>% arrange(EBV) %>% slice(topidx[[i]]) %>% summarise(x = mean(TBV)) %>% unlist()
                botvals[i] = d1 %>% arrange(EBV) %>% slice(botidx[[i]]) %>% summarise(x = mean(TBV)) %>% unlist()
            }
            cat(poolsize[ps], sigmavec[ss], 'Extremes', 'Cov', 'Unequal', irep, mean(DF$TBV), 
                c(rev(topvals), botvals),
                sep = " ", file = '/work/rmlewis/napov/PoolingSimulation/G3/TopBotSel.txt', append = T)
            cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/TopBotSel.txt', append = T)
            tic()
            MME                 = BuildMME(npool = K, nsires = nsires, sigma2e = 0.7, sigma2a = 0.3,
                                           G = cov2cor(G), y_k, rep(n_k[k], K))
            ans                 = solve(MME$LHS, MME$rhs)
            C22                 = diag(ginv(as.matrix(MME$LHS)))[2:(nsires + 1)]
            rel                 = sqrt(pmax(1 - C22*(0.7/0.3), 0))
            toc()
            c1                  = cor(SireTBV$TBV[which(SireTBV$ID %in% rsires)], ans[2:(nsires + 1)])
            cat(c1, irep, poolsize[ps], sigmavec[ss], 'Extremes', 'Cor',
                mean(siresperpool), mean(progpersire), mean(rowSums(SireRep)), mean(rel), sd(rel), 'Unequal',
                sep = " ", file = '/work/rmlewis/napov/PoolingSimulation/G3/SelSim1.txt', append = T)
            cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/SelSim1.txt', append = T)
            d1                  = data.frame(TBV = SireTBV$TBV[which(SireTBV$ID %in% rsires)], 
                                             EBV = ans[2:(nsires + 1)])
            topvals             = rep(0, 3)
            topidx              = list(19:30, 25:30, 28:30)
            botvals             = rep(0, 3)
            botidx              = list(1:12, 1:6, 1:3) 
            for(i in 1:3){
                topvals[i] = d1 %>% arrange(EBV) %>% slice(topidx[[i]]) %>% summarise(x = mean(TBV)) %>% unlist()
                botvals[i] = d1 %>% arrange(EBV) %>% slice(botidx[[i]]) %>% summarise(x = mean(TBV)) %>% unlist()
            }
            cat(poolsize[ps], sigmavec[ss], 'Extremes', 'Cor', 'Unequal', irep, mean(DF$TBV), 
                c(rev(topvals), botvals),
                sep = " ", file = '/work/rmlewis/napov/PoolingSimulation/G3/TopBotSel.txt', append = T)
            cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/TopBotSel.txt', append = T)
        }
    }
    #-----------------------#
    # No construction error #
    #-----------------------#
    for(ps in 1:NP){
        for(ss in 1:NS){
            K                   = length(subid)/poolsize[ps]
            W                   = matrix(0, nrow = K, ncol = nsnp)
            n_k                 = rep(floor(nindpool/K), K)
            y_k                 = c()
            idq                 = 1:n_k[1]
            siresperpool        = c()
            progpersire         = c()
            SireRep             = array(0, dim = c(nsires, K))
            for(k in 1:K){
                pivec           = rep(1/n_k[k], n_k[k])
                idvec           = SDF[idq]
                siresperpool[k] = sum(rsires%in%DFF$Sire[which(DFF$ID%in%idvec)])
                progpersire[k]  = sum(DFF$ID[DFF$Sire%in%rsires[rsires%in%DFF$Sire[which(DFF$ID%in%idvec)]]]%in%idvec)
                SireRep[, k]    = (rsires%in%DFF$Sire[which(DFF$ID%in%idvec)])*1
                wtmp            = Genotypes[which(DF$ID%in%idvec), ]
                W[k, ]          = LogitNoise(0.5*colSums(pivec*wtmp), sigmavec[ss]^2)
                y_k[k]          = mean(DFF$Pheno[which(DFF$ID%in%idvec)])
                idq             = idq + rep(1, n_k[1])*n_k[1]
                cat('Pool: ', k, ' constructed\n')
            }
            cat("\n")
            idx                 = which(rownames(SireGeno)%in%rsires)
            X                   = SireGeno[idx, ]
            phat                = colMeans(X)/2
            muhat               = colMeans(W)
            indx                = which(phat < 0.01 | phat > 0.99)
            nrem                = nsnp - length(indx)
            M                   = rbind(X[, -indx]/2, W[, -indx])
            muvec               = colMeans(M)
            G                   = GetGPool(M, muvec, dim(M)[1], dim(M)[2])
            diag(G)             = diag(G)*1.01
            tic()
            MME                 = BuildMME(npool = K, nsires = nsires, sigma2e = 0.7, sigma2a = 0.3, 
                                           G = G, y_k, rep(n_k[k], K))
            ans                 = solve(MME$LHS, MME$rhs)
            C22                 = diag(ginv(as.matrix(MME$LHS)))[2:(nsires + 1)]
            rel                 = sqrt(pmax(1 - C22*(0.7/0.3), 0)) 
            toc()
            c1                  = cor(SireTBV$TBV[which(SireTBV$ID %in% rsires)], ans[2:(nsires + 1)])
            cat(c1, irep, poolsize[ps], sigmavec[ss], 'Extremes', 'Cov',
                mean(siresperpool), mean(progpersire), mean(rowSums(SireRep)), mean(rel), sd(rel), 'Equal',
                sep = " ", file = '/work/rmlewis/napov/PoolingSimulation/G3/SelSim1.txt', append = T)
            cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/SelSim1.txt', append = T)
            d1                  = data.frame(TBV = SireTBV$TBV[which(SireTBV$ID %in% rsires)], 
                                             EBV = ans[2:(nsires + 1)])
            topvals             = rep(0, 3)
            topidx              = list(19:30, 25:30, 28:30)
            botvals             = rep(0, 3)
            botidx              = list(1:12, 1:6, 1:3) 
            for(i in 1:3){
                topvals[i] = d1 %>% arrange(EBV) %>% slice(topidx[[i]]) %>% summarise(x = mean(TBV)) %>% unlist()
                botvals[i] = d1 %>% arrange(EBV) %>% slice(botidx[[i]]) %>% summarise(x = mean(TBV)) %>% unlist()
            }
            cat(poolsize[ps], sigmavec[ss], 'Extremes', 'Cov', 'Equal', irep, mean(DF$TBV), 
                c(rev(topvals), botvals),
                sep = " ", file = '/work/rmlewis/napov/PoolingSimulation/G3/TopBotSel.txt', append = T)
            cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/TopBotSel.txt', append = T)
            tic()
            MME                 = BuildMME(npool = K, nsires = nsires, sigma2e = 0.7, sigma2a = 0.3, 
                                           G = cov2cor(G), y_k, rep(n_k[k], K))
            ans                 = solve(MME$LHS, MME$rhs)
            C22                 = diag(ginv(as.matrix(MME$LHS)))[2:(nsires + 1)]
            rel                 = sqrt(pmax(1 - C22*(0.7/0.3), 0))
            toc()
            c1                  = cor(SireTBV$TBV[which(SireTBV$ID %in% rsires)], ans[2:(nsires + 1)])
            cat(c1, irep, poolsize[ps], sigmavec[ss], 'Extremes', 'Cor',
                mean(siresperpool), mean(progpersire), mean(rowSums(SireRep)), mean(rel), sd(rel), 'Equal',
                sep = " ", file = '/work/rmlewis/napov/PoolingSimulation/G3/SelSim1.txt', append = T)
            cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/SelSim1.txt', append = T)
            d1                  = data.frame(TBV = SireTBV$TBV[which(SireTBV$ID %in% rsires)], 
                                             EBV = ans[2:(nsires + 1)])
            topvals             = rep(0, 3)
            topidx              = list(19:30, 25:30, 28:30)
            botvals             = rep(0, 3)
            botidx              = list(1:12, 1:6, 1:3) 
            for(i in 1:3){
                topvals[i] = d1 %>% arrange(EBV) %>% slice(topidx[[i]]) %>% summarise(x = mean(TBV)) %>% unlist()
                botvals[i] = d1 %>% arrange(EBV) %>% slice(botidx[[i]]) %>% summarise(x = mean(TBV)) %>% unlist()
            }
            cat(poolsize[ps], sigmavec[ss], 'Extremes', 'Cor', 'Equal', irep, mean(DF$TBV), 
                c(rev(topvals), botvals),
                sep = " ", file = '/work/rmlewis/napov/PoolingSimulation/G3/TopBotSel.txt', append = T)
            cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/TopBotSel.txt', append = T)
        }
    }
    #----------------#
    # Random Pooling #
    #----------------#
    subid       = sample(nrow(DF), 3000)
    DFF         = DF[subid, ]
    SDF         = DFF$ID
    for(ps in 1:NP){
        for(ss in 1:NS){
            DFF                 = DF[subid, ]
            SDF                 = DFF$ID
            K                   = length(subid)/poolsize[ps]
            W                   = matrix(0, nrow = K, ncol = nsnp)
            n_k                 = rep(floor(nindpool/K), K)
            y_k                 = c()
            idq                 = 1:n_k[1]
            siresperpool        = c()
            progpersire         = c()
            SireRep             = array(0, dim = c(nsires, K))
            for(k in 1:K){
                pivec           = c(rdirichlet(1, rep(10, n_k[k])))
                tmp             = SDF
                idvec           = sort(sample(SDF, n_k[k]))
                SDF             = sort(setdiff(tmp, idvec))
                siresperpool[k] = sum(rsires%in%DFF$Sire[which(DFF$ID%in%idvec)])
                progpersire[k]  = sum(DFF$ID[DFF$Sire%in%rsires[rsires%in%DFF$Sire[which(DFF$ID%in%idvec)]]]%in%idvec)
                wtmp            = Genotypes[which(DF$ID%in%idvec), ]
                W[k, ]          = LogitNoise(0.5*colSums(pivec*wtmp), sigmavec[ss]^2)
                y_k[k]          = mean(DFF$Pheno[which(DFF$ID%in%idvec)])
                cat('Pool: ', k, ' constructed\n')
            }
            cat("\n")
            idx                 = which(rownames(SireGeno)%in%rsires)
            X                   = SireGeno[idx, ]
            phat                = colMeans(X)/2
            muhat               = colMeans(W)
            indx                = which(phat < 0.01 | phat > 0.99)
            nrem                = nsnp - length(indx)
            M                   = rbind(X[, -indx]/2, W[, -indx])
            muvec               = colMeans(M)
            G                   = GetGPool(M, muvec, dim(M)[1], dim(M)[2])
            diag(G)             = diag(G)*1.01
            tic()
            MME                 = BuildMME(npool = K, nsires = nsires, sigma2e = 0.7, sigma2a = 0.3, 
                                           G = G, y_k, rep(n_k[k], K))
            ans                 = solve(MME$LHS, MME$rhs)
            C22                 = diag(ginv(as.matrix(MME$LHS)))[2:(nsires + 1)]
            rel                 = sqrt(pmax(1 - C22*(0.7/0.3), 0))
            toc()
            c1                  = cor(SireTBV$TBV[which(SireTBV$ID %in% rsires)], ans[2:(nsires + 1)])
            cat(c1, irep, poolsize[ps], sigmavec[ss], 'Random', 'Cov',
                mean(siresperpool), mean(progpersire), mean(rowSums(SireRep)), mean(rel), sd(rel), 'Unequal',
                sep = " ", file = '/work/rmlewis/napov/PoolingSimulation/G3/SelSim1.txt', append = T)
            cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/SelSim1.txt', append = T)
            d1                  = data.frame(TBV = SireTBV$TBV[which(SireTBV$ID %in% rsires)], 
                                             EBV = ans[2:(nsires + 1)])
            topvals             = rep(0, 3)
            topidx              = list(19:30, 25:30, 28:30)
            botvals             = rep(0, 3)
            botidx              = list(1:12, 1:6, 1:3) 
            for(i in 1:3){
                topvals[i] = d1 %>% arrange(EBV) %>% slice(topidx[[i]]) %>% summarise(x = mean(TBV)) %>% unlist()
                botvals[i] = d1 %>% arrange(EBV) %>% slice(botidx[[i]]) %>% summarise(x = mean(TBV)) %>% unlist()
            }
            cat(poolsize[ps], sigmavec[ss], 'Random', 'Cov', 'Unequal', irep, mean(DF$TBV), 
                c(rev(topvals), botvals),
                sep = " ", file = '/work/rmlewis/napov/PoolingSimulation/G3/TopBotSel.txt', append = T)
            cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/TopBotSel.txt', append = T)
            tic()
            MME                 = BuildMME(npool = K, nsires = nsires, sigma2e = 0.7, sigma2a = 0.3, 
                                           G = cov2cor(G), y_k, rep(n_k[k], K))
            ans                 = solve(MME$LHS, MME$rhs)
            C22                 = diag(ginv(as.matrix(MME$LHS)))[2:(nsires + 1)]
            rel                 = sqrt(pmax(1 - C22*(0.7/0.3), 0))
            toc()
            c1                  = cor(SireTBV$TBV[which(SireTBV$ID %in% rsires)], ans[2:(nsires + 1)])
            cat(c1, irep, poolsize[ps], sigmavec[ss], 'Random', 'Cor',
                mean(siresperpool), mean(progpersire), mean(rowSums(SireRep)), mean(rel), sd(rel), 'Unequal',
                sep = " ", file = '/work/rmlewis/napov/PoolingSimulation/G3/SelSim1.txt', append = T)
            cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/SelSim1.txt', append = T)
            d1                  = data.frame(TBV = SireTBV$TBV[which(SireTBV$ID %in% rsires)], 
                                             EBV = ans[2:(nsires + 1)])
            topvals             = rep(0, 3)
            topidx              = list(19:30, 25:30, 28:30)
            botvals             = rep(0, 3)
            botidx              = list(1:12, 1:6, 1:3) 
            for(i in 1:3){
                topvals[i] = d1 %>% arrange(EBV) %>% slice(topidx[[i]]) %>% summarise(x = mean(TBV)) %>% unlist()
                botvals[i] = d1 %>% arrange(EBV) %>% slice(botidx[[i]]) %>% summarise(x = mean(TBV)) %>% unlist()
            }
            cat(poolsize[ps], sigmavec[ss], 'Random', 'Cor', 'Unequal', irep, mean(DF$TBV), 
                c(rev(topvals), botvals),
                sep = " ", file = '/work/rmlewis/napov/PoolingSimulation/G3/TopBotSel.txt', append = T)
            cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/TopBotSel.txt', append = T)
        }
    }
    #-----------------------#
    # No construction error #
    #-----------------------#
    for(ps in 1:NP){
        for(ss in 1:NS){
            DFF                 = DF[subid, ]
            SDF                 = DFF$ID
            K                   = length(subid)/poolsize[ps]
            W                   = matrix(0, nrow = K, ncol = nsnp)
            n_k                 = rep(floor(nindpool/K), K)
            y_k                 = c()
            idq                 = 1:n_k[1]
            siresperpool        = c()
            progpersire         = c()
            SireRep             = array(0, dim = c(nsires, K))
            for(k in 1:K){
                pivec           = rep(1/n_k[k], n_k[k])
                tmp             = SDF
                idvec           = sort(sample(SDF, n_k[k]))
                SDF             = sort(setdiff(tmp, idvec))
                siresperpool[k] = sum(rsires%in%DFF$Sire[which(DFF$ID%in%idvec)])
                progpersire[k]  = sum(DFF$ID[DFF$Sire%in%rsires[rsires%in%DFF$Sire[which(DFF$ID%in%idvec)]]]%in%idvec)
                wtmp            = Genotypes[which(DF$ID%in%idvec), ]
                W[k, ]          = LogitNoise(0.5*colSums(pivec*wtmp), sigmavec[ss]^2)
                y_k[k]          = mean(DFF$Pheno[which(DFF$ID%in%idvec)])
                cat('Pool: ', k, ' constructed\n')
            }
            cat("\n")
            idx                 = which(rownames(SireGeno)%in%rsires)
            X                   = SireGeno[idx, ]
            phat                = colMeans(X)/2
            muhat               = colMeans(W)
            indx                = which(phat < 0.01 | phat > 0.99)
            nrem                = nsnp - length(indx)
            M                   = rbind(X[, -indx]/2, W[, -indx])
            muvec               = colMeans(M)
            G                   = GetGPool(M, muvec, dim(M)[1], dim(M)[2])
            diag(G)             = diag(G)*1.01
            tic()
            MME                 = BuildMME(npool = K, nsires = nsires, sigma2e = 0.7, sigma2a = 0.3, 
                                           G = G, y_k, rep(n_k[k], K))
            ans                 = solve(MME$LHS, MME$rhs)
            C22                 = diag(ginv(as.matrix(MME$LHS)))[2:(nsires + 1)]
            rel                 = sqrt(pmax(1 - C22*(0.7/0.3), 0))
            toc()
            c1                  = cor(SireTBV$TBV[which(SireTBV$ID %in% rsires)], ans[2:(nsires + 1)])
            cat(c1, irep, poolsize[ps], sigmavec[ss], 'Random', 'Cov',
                mean(siresperpool), mean(progpersire), mean(rowSums(SireRep)), mean(rel), sd(rel), 'Equal',
                sep = " ", file = '/work/rmlewis/napov/PoolingSimulation/G3/SelSim1.txt', append = T)
            cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/SelSim1.txt', append = T)
            d1                  = data.frame(TBV = SireTBV$TBV[which(SireTBV$ID %in% rsires)], 
                                             EBV = ans[2:(nsires + 1)])
            topvals             = rep(0, 3)
            topidx              = list(19:30, 25:30, 28:30)
            botvals             = rep(0, 3)
            botidx              = list(1:12, 1:6, 1:3) 
            for(i in 1:3){
                topvals[i] = d1 %>% arrange(EBV) %>% slice(topidx[[i]]) %>% summarise(x = mean(TBV)) %>% unlist()
                botvals[i] = d1 %>% arrange(EBV) %>% slice(botidx[[i]]) %>% summarise(x = mean(TBV)) %>% unlist()
            }
            cat(poolsize[ps], sigmavec[ss], 'Random', 'Cov', 'Equal', irep, mean(DF$TBV), 
                c(rev(topvals), botvals),
                sep = " ", file = '/work/rmlewis/napov/PoolingSimulation/G3/TopBotSel.txt', append = T)
            cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/TopBotSel.txt', append = T)
            tic()
            MME                 = BuildMME(npool = K, nsires = nsires, sigma2e = 0.7, sigma2a = 0.3, 
                                           G = cov2cor(G), y_k, rep(n_k[k], K))
            ans                 = solve(MME$LHS, MME$rhs)
            C22                 = diag(ginv(as.matrix(MME$LHS)))[2:(nsires + 1)]
            rel                 = sqrt(pmax(1 - C22*(0.7/0.3), 0))
            toc()
            c1                  = cor(SireTBV$TBV[which(SireTBV$ID %in% rsires)], ans[2:(nsires + 1)])
            cat(c1, irep, poolsize[ps], sigmavec[ss], 'Random', 'Cor',
                mean(siresperpool), mean(progpersire), mean(rowSums(SireRep)), mean(rel), sd(rel), 'Equal',
                sep = " ", file = '/work/rmlewis/napov/PoolingSimulation/G3/SelSim1.txt', append = T)
            cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/SelSim1.txt', append = T)
            d1                  = data.frame(TBV = SireTBV$TBV[which(SireTBV$ID %in% rsires)], 
                                             EBV = ans[2:(nsires + 1)])
            topvals             = rep(0, 3)
            topidx              = list(19:30, 25:30, 28:30)
            botvals             = rep(0, 3)
            botidx              = list(1:12, 1:6, 1:3) 
            for(i in 1:3){
                topvals[i] = d1 %>% arrange(EBV) %>% slice(topidx[[i]]) %>% summarise(x = mean(TBV)) %>% unlist()
                botvals[i] = d1 %>% arrange(EBV) %>% slice(botidx[[i]]) %>% summarise(x = mean(TBV)) %>% unlist()
            }
            cat(poolsize[ps], sigmavec[ss], 'Random', 'Cor', 'Equal', irep, mean(DF$TBV), 
                c(rev(topvals), botvals),
                sep = " ", file = '/work/rmlewis/napov/PoolingSimulation/G3/TopBotSel.txt', append = T)
            cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/TopBotSel.txt', append = T)
        }
    }
    #-----------------------------------------#
    # K-means pooling with construction error #
    #-----------------------------------------#
    DFF         = DF
    DFF         = arrange(DFF, desc(Pheno))
    for(ps in 1:NP){
        for(ss in 1:NS){
            K                   = length(subid)/poolsize[ps]
            SDF                 = DFF$ID
            W                   = matrix(0, nrow = K, ncol = nsnp)
            y_k                 = c()
            m1                  = kmeans(DFF$Pheno, centers = 3*K)
            DFclus              = data.frame(Pool = 1:(3*K), Val = m1$centers)
            DFC                 = arrange(DFclus, desc(Val))
            cvec                = m1$cluster
            DFC$Size            = sapply(DFC$Pool, function(x) sum(cvec == x))
            m2                  = rbind(head(DFC, K/2), tail(DFC, K/2))
            y_k                 = c()
            siresperpool        = c()
            progpersire         = c()
            SireRep             = array(0, dim = c(nsires, K))
            nvec                = c()
            for(k in 1:K){
                idq             = cvec == m2$Pool[k]
                nvec[k]         = sum(idq)
                idvec           = SDF[idq]
                pivec           = c(rdirichlet(1, rep(10, nvec[k])))
                siresperpool[k] = sum(rsires%in%DFF$Sire[which(DFF$ID%in%idvec)])
                progpersire[k]  = sum(DFF$ID[DFF$Sire%in%rsires[rsires%in%DFF$Sire[which(DFF$ID%in%idvec)]]]%in%idvec)
                wtmp            = Genotypes[which(DF$ID%in%idvec), ]
                if(sum(idq) == 1) {
                    W[k, ]      = LogitNoise(0.5*wtmp, sigmavec[ss]^2)
                } else {
                    W[k, ]      = LogitNoise(0.5*colSums(pivec*wtmp), sigmavec[ss]^2)
                }
                y_k[k]          = m2$Val[k]
                cat('Pool: ', k, ' constructed\n')
            }
            cat("\n")
            idx                 = which(rownames(SireGeno)%in%rsires)
            X                   = SireGeno[idx, ]
            phat                = colMeans(X)/2
            muhat               = colMeans(W)
            indx                = which(phat < 0.01 | phat > 0.99)
            nrem                = nsnp - length(indx)
            M                   = rbind(X[, -indx]/2, W[, -indx])
            muvec               = colMeans(M)
            G                   = GetGPool(M, muvec, dim(M)[1], dim(M)[2])
            diag(G)             = diag(G)*1.01
            tic()
            MME                 = BuildMME(npool = K, nsires = nsires, sigma2e = 0.7, sigma2a = 0.3, 
                                           G = G, y_k, nvec)
            ans                 = solve(MME$LHS, MME$rhs)
            C22                 = diag(ginv(as.matrix(MME$LHS)))[2:(nsires + 1)]
            rel                 = sqrt(pmax(1 - C22*(0.7/0.3), 0))
            toc()
            c1                  = cor(SireTBV$TBV[which(SireTBV$ID %in% rsires)], ans[2:(nsires + 1)])
            cat(c1, irep, poolsize[ps], sigmavec[ss], 'Kmeans', 'Cov',
                mean(siresperpool), mean(progpersire), mean(rowSums(SireRep)), mean(rel), sd(rel), 'Unequal',
                sep = " ", file = '/work/rmlewis/napov/PoolingSimulation/G3/SelSim1.txt', append = T)
            cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/SelSim1.txt', append = T)
            d1                  = data.frame(TBV = SireTBV$TBV[which(SireTBV$ID %in% rsires)], 
                                             EBV = ans[2:(nsires + 1)])
            topvals             = rep(0, 3)
            topidx              = list(19:30, 25:30, 28:30)
            botvals             = rep(0, 3)
            botidx              = list(1:12, 1:6, 1:3) 
            for(i in 1:3){
                topvals[i] = d1 %>% arrange(EBV) %>% slice(topidx[[i]]) %>% summarise(x = mean(TBV)) %>% unlist()
                botvals[i] = d1 %>% arrange(EBV) %>% slice(botidx[[i]]) %>% summarise(x = mean(TBV)) %>% unlist()
            }
            cat(poolsize[ps], sigmavec[ss], 'Kmeans', 'Cov', 'Unequal', irep, mean(DF$TBV), 
                c(rev(topvals), botvals),
                sep = " ", file = '/work/rmlewis/napov/PoolingSimulation/G3/TopBotSel.txt', append = T)
            cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/TopBotSel.txt', append = T)
            tic()
            MME                 = BuildMME(npool = K, nsires = nsires, sigma2e = 0.7, sigma2a = 0.3, 
                                           G = cov2cor(G), y_k, nvec)
            ans                 = solve(MME$LHS, MME$rhs)
            C22                 = diag(ginv(as.matrix(MME$LHS)))[2:(nsires + 1)]
            rel                 = sqrt(pmax(1 - C22*(0.7/0.3), 0))
            toc()
            c1                  = cor(SireTBV$TBV[which(SireTBV$ID %in% rsires)], ans[2:(nsires + 1)])
            cat(c1, irep, poolsize[ps], sigmavec[ss], 'Kmeans', 'Cor',
                mean(siresperpool), mean(progpersire), mean(rowSums(SireRep)), mean(rel), sd(rel), 'Unequal',
                sep = " ", file = '/work/rmlewis/napov/PoolingSimulation/G3/SelSim1.txt', append = T)
            cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/SelSim1.txt', append = T)
            cat(mean(m2$Size), sd(m2$Size), range(m2$Size), sum(m2$Size == 1), sum(m2$Size),
                irep, K, sigmavec[ss], sep = " ", 
                file = '/work/rmlewis/napov/PoolingSimulation/G3/SelClus1.txt', append = T)
            cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/SelClus1.txt', append = T)
            d1                  = data.frame(TBV = SireTBV$TBV[which(SireTBV$ID %in% rsires)], 
                                             EBV = ans[2:(nsires + 1)])
            topvals             = rep(0, 3)
            topidx              = list(19:30, 25:30, 28:30)
            botvals             = rep(0, 3)
            botidx              = list(1:12, 1:6, 1:3) 
            for(i in 1:3){
                topvals[i] = d1 %>% arrange(EBV) %>% slice(topidx[[i]]) %>% summarise(x = mean(TBV)) %>% unlist()
                botvals[i] = d1 %>% arrange(EBV) %>% slice(botidx[[i]]) %>% summarise(x = mean(TBV)) %>% unlist()
            }
            cat(poolsize[ps], sigmavec[ss], 'Kmeans', 'Cor', 'Unequal', irep, mean(DF$TBV), 
                c(rev(topvals), botvals),
                sep = " ", file = '/work/rmlewis/napov/PoolingSimulation/G3/TopBotSel.txt', append = T)
            cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/TopBotSel.txt', append = T)
        }
    }
    #-----------------------#
    # No construction error #
    #-----------------------#
    for(ps in 1:NP){
        for(ss in 1:NS){
            K                   = length(subid)/poolsize[ps]
            SDF                 = DFF$ID
            W                   = matrix(0, nrow = K, ncol = nsnp)
            y_k                 = c()
            m1                  = kmeans(DFF$Pheno, centers = 3*K)
            DFclus              = data.frame(Pool = 1:(3*K), Val = m1$centers)
            DFC                 = arrange(DFclus, desc(Val))
            cvec                = m1$cluster
            DFC$Size            = sapply(DFC$Pool, function(x) sum(cvec == x))
            m2                  = rbind(head(DFC, K/2), tail(DFC, K/2))
            y_k                 = c()
            siresperpool        = c()
            progpersire         = c()
            SireRep             = array(0, dim = c(nsires, K))
            nvec                = c()
            for(k in 1:K){
                idq             = cvec == m2$Pool[k]
                nvec[k]         = sum(idq)
                idvec           = SDF[idq]
                pivec           = rep(1/nvec[k], nvec[k])
                siresperpool[k] = sum(rsires%in%DFF$Sire[which(DFF$ID%in%idvec)])
                progpersire[k]  = sum(DFF$ID[DFF$Sire%in%rsires[rsires%in%DFF$Sire[which(DFF$ID%in%idvec)]]]%in%idvec)
                wtmp            = Genotypes[which(DF$ID%in%idvec), ]
                if(sum(idq) == 1) {
                    W[k, ]      = LogitNoise(0.5*wtmp, sigmavec[ss]^2)
                } else {
                    W[k, ]      = LogitNoise(0.5*colSums(pivec*wtmp), sigmavec[ss]^2)
                }
                y_k[k]          = m2$Val[k]
                cat('Pool: ', k, ' constructed\n')
            }
            cat("\n")
            idx                 = which(rownames(SireGeno)%in%rsires)
            X                   = SireGeno[idx, ]
            phat                = colMeans(X)/2
            muhat               = colMeans(W)
            indx                = which(phat < 0.01 | phat > 0.99)
            nrem                = nsnp - length(indx)
            M                   = rbind(X[, -indx]/2, W[, -indx])
            muvec               = colMeans(M)
            G                   = GetGPool(M, muvec, dim(M)[1], dim(M)[2])
            diag(G)             = diag(G)*1.01
            tic()
            MME                 = BuildMME(npool = K, nsires = nsires, sigma2e = 0.7, sigma2a = 0.3, 
                                           G = G, y_k, nvec)
            ans                 = solve(MME$LHS, MME$rhs)
            C22                 = diag(ginv(as.matrix(MME$LHS)))[2:(nsires + 1)]
            rel                 = sqrt(pmax(1 - C22*(0.7/0.3), 0))
            toc()
            c1                  = cor(SireTBV$TBV[which(SireTBV$ID %in% rsires)], ans[2:(nsires + 1)])
            cat(c1, irep, poolsize[ps], sigmavec[ss], 'Kmeans', 'Cov',
                mean(siresperpool), mean(progpersire), mean(rowSums(SireRep)), mean(rel), sd(rel), 'Equal',
                sep = " ", file = '/work/rmlewis/napov/PoolingSimulation/G3/SelSim1.txt', append = T)
            cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/SelSim1.txt', append = T)
            d1                  = data.frame(TBV = SireTBV$TBV[which(SireTBV$ID %in% rsires)], 
                                             EBV = ans[2:(nsires + 1)])
            topvals             = rep(0, 3)
            topidx              = list(19:30, 25:30, 28:30)
            botvals             = rep(0, 3)
            botidx              = list(1:12, 1:6, 1:3) 
            for(i in 1:3){
                topvals[i] = d1 %>% arrange(EBV) %>% slice(topidx[[i]]) %>% summarise(x = mean(TBV)) %>% unlist()
                botvals[i] = d1 %>% arrange(EBV) %>% slice(botidx[[i]]) %>% summarise(x = mean(TBV)) %>% unlist()
            }
            cat(poolsize[ps], sigmavec[ss], 'Kmeans', 'Cov', 'Equal', irep, mean(DF$TBV), 
                c(rev(topvals), botvals),
                sep = " ", file = '/work/rmlewis/napov/PoolingSimulation/G3/TopBotSel.txt', append = T)
            cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/TopBotSel.txt', append = T)
            tic()
            MME                 = BuildMME(npool = K, nsires = nsires, sigma2e = 0.7, sigma2a = 0.3,
                                           G = cov2cor(G), y_k, nvec)
            ans                 = solve(MME$LHS, MME$rhs)
            C22                 = diag(ginv(as.matrix(MME$LHS)))[2:(nsires + 1)]
            rel                 = sqrt(pmax(1 - C22*(0.7/0.3), 0))
            toc()
            c1                  = cor(SireTBV$TBV[which(SireTBV$ID %in% rsires)], ans[2:(nsires + 1)])
            cat(c1, irep, poolsize[ps], sigmavec[ss], 'Kmeans', 'Cor',
                mean(siresperpool), mean(progpersire), mean(rowSums(SireRep)), mean(rel), sd(rel), 'Equal',
                sep = " ", file = '/work/rmlewis/napov/PoolingSimulation/G3/SelSim1.txt', append = T)
            cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/SelSim1.txt', append = T)
            cat(mean(m2$Size), sd(m2$Size), range(m2$Size), sum(m2$Size == 1), sum(m2$Size),
                irep, K, sigmavec[ss], sep = " ", 
                file = '/work/rmlewis/napov/PoolingSimulation/G3/SelClus1.txt', append = T)
            cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/SelClus1.txt', append = T)
            d1                  = data.frame(TBV = SireTBV$TBV[which(SireTBV$ID %in% rsires)], 
                                             EBV = ans[2:(nsires + 1)])
            topvals             = rep(0, 3)
            topidx              = list(19:30, 25:30, 28:30)
            botvals             = rep(0, 3)
            botidx              = list(1:12, 1:6, 1:3) 
            for(i in 1:3){
                topvals[i] = d1 %>% arrange(EBV) %>% slice(topidx[[i]]) %>% summarise(x = mean(TBV)) %>% unlist()
                botvals[i] = d1 %>% arrange(EBV) %>% slice(botidx[[i]]) %>% summarise(x = mean(TBV)) %>% unlist()
            }
            cat(poolsize[ps], sigmavec[ss], 'Kmeans', 'Cor', 'Equal', irep, mean(DF$TBV), 
                c(rev(topvals), botvals),
                sep = " ", file = '/work/rmlewis/napov/PoolingSimulation/G3/TopBotSel.txt', append = T)
            cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/TopBotSel.txt', append = T)
        }
    }
    #-------------------------------------------#
    # GBLUP with all (subset) progeny genotyped #
    #-------------------------------------------#
    idx     = which(rownames(SireGeno)%in%rsires)
    DFF     = DF[subid, ]
    idi     = which(rownames(Genotypes)%in%DFF$ID)
    phat    = colMeans(Genotypes)/2
    indx    = which(phat < 0.01 | phat > 0.99)
    nrem    = nsnp - length(indx)
    M       = rbind(SireGeno[idx, -indx], Genotypes[idi, -indx])
    muvec   = colMeans(M)/2
    G       = GetGInt(M, muvec, dim(M)[1], dim(M)[2])
    MME     = BuildMME_i(length(subid), nsires, 0.7, 0.3, arrange(DFF, ID)$Pheno, G)
    ans     = solve(MME$LHS, MME$rhs)
    C22     = diag(ginv(as.matrix(MME$LHS)))[2:(nsires + 1)]
    rel     = sqrt(pmax(1 - C22*(0.7/0.3), 0))
    c1      = cor(SireTBV$TBV[which(SireTBV$ID %in% rsires)], ans[2:(nsires + 1)])
    cat(c1, irep, '--', '--', 'GBLUP', '--',
        '--', '--', '--', mean(rel), sd(rel), '--',
        sep = " ", file = '/work/rmlewis/napov/PoolingSimulation/G3/SelSim1.txt', append = T)
    cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/SelSim1.txt', append = T)
    d1                  = data.frame(TBV = SireTBV$TBV[which(SireTBV$ID %in% rsires)], 
                                     EBV = ans[2:(nsires + 1)])
    topvals             = rep(0, 3)
    topidx              = list(19:30, 25:30, 28:30)
    botvals             = rep(0, 3)
    botidx              = list(1:12, 1:6, 1:3) 
    for(i in 1:3){
        topvals[i] = d1 %>% arrange(EBV) %>% slice(topidx[[i]]) %>% summarise(x = mean(TBV)) %>% unlist()
        botvals[i] = d1 %>% arrange(EBV) %>% slice(botidx[[i]]) %>% summarise(x = mean(TBV)) %>% unlist()
    }
    cat('--', '--', 'GBLUP', '--', '--', irep, mean(DF$TBV), 
        c(rev(topvals), botvals),
        sep = " ", file = '/work/rmlewis/napov/PoolingSimulation/G3/TopBotSel.txt', append = T)
    cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/TopBotSel.txt', append = T)
    rm(G)
    rm(MME)
    rm(SireGeno)
    rm(Genotypes)
    #---------------------------#
    # PBLLUP with full pedigree #
    #---------------------------#
    ped1    = buildPed(DF$ID, DF$Sire, DF$Dam, add.ancestors = T)[, 1:3]
    ped2    = renumped(ped1)
    Ainv    = nadiv::makeAinv(ped1)
    xx1     = as(Ainv$Ainv, 'dgTMatrix')
    xmat    = rep(1, nrow(DF))
    zmat    = cbind(Matrix(0, ncol = nrow(ped2) - nrow(DF), nrow = nrow(DF)), Diagonal(nrow(DF)))
    tmp     = PBLUP(X = xmat, Z = zmat, Ainv = xx1, Y = DF$Pheno, alpha = 0.7/0.3)
    ans     = tmp$sol[-1,]
    C22     = diag(Matrix::solve(tmp$C, Diagonal(ncol(tmp$C))))[-1]
    c22     = C22[which(ped2$Old%in%rsires)]
    rel     = sqrt(pmax(1 - c22*(0.7/0.3), 0))
    d1      = data.frame(ID = ped2$Old[ped2$Old%in%rsires], EBV = ans[which(ped2$Old%in%rsires)])
    d2      = rbind(SireTBV, ProgTBV)
    df      = merge(d1, d2, by = 'ID')
    c1      = cor(df$TBV, df$EBV)
    cat(c1, irep, '--', '--', 'PBLUP', '--',
        '--', '--', '--', mean(rel), sd(rel), '--',
        sep = " ", file = '/work/rmlewis/napov/PoolingSimulation/G3/SelSim1.txt', append = T)
    cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/SelSim1.txt', append = T)
    #cat('Selection', irep, mean(LDmat), sd(LDmat), max(LDmat), min(LDmat), sep = ' ', 
    #    file = '/work/rmlewis/napov/PoolingSimulation/Selection/LDsummary2.txt', append = T)
    #cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/Selection/LDsummary2.txt', append = T)
    d1                  = data.frame(TBV = SireTBV$TBV[which(SireTBV$ID %in% rsires)], 
                                     EBV = ans[2:(nsires + 1)])
    topvals             = rep(0, 3)
    topidx              = list(19:30, 25:30, 28:30)
    botvals             = rep(0, 3)
    botidx              = list(1:12, 1:6, 1:3) 
    for(i in 1:3){
        topvals[i] = d1 %>% arrange(EBV) %>% slice(topidx[[i]]) %>% summarise(x = mean(TBV)) %>% unlist()
        botvals[i] = d1 %>% arrange(EBV) %>% slice(botidx[[i]]) %>% summarise(x = mean(TBV)) %>% unlist()
    }
    cat('--', '--', 'PBLUP', '--', '--', irep, mean(DF$TBV), 
        c(rev(topvals), botvals),
        sep = " ", file = '/work/rmlewis/napov/PoolingSimulation/G3/TopBotSel.txt', append = T)
    cat('\n', file = '/work/rmlewis/napov/PoolingSimulation/G3/TopBotSel.txt', append = T)
    cat('#----------------------------------#\n')
    cat('# Replication :', irep, ' completed.\n')
    cat('#----------------------------------#\n')
    rm(tmp)
    rm(C22)
    gc()
}
t2 = proc.time()
tt = t2 - t1; tt