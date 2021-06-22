library(AlphaSimR)
library(dplyr)
library(MCMCpack)
library(tictoc)
library(Matrix)
# Inverse Logit function
'InvLogit'  = function(x){
    exp(x)/(1 + exp(x))
}
# Function to add gaussian noise
'LogitNoise' = function(x, sigma2){
    x   = ifelse(x > 0.99, 0.99, ifelse(x < 0.01, 0.01, x))
    y   = log(x/(1 - x))
    tmp = rnorm(length(x))*sqrt(sigma2) + y
    z   = ifelse(tmp > 4.0, 4.0 - (tmp - 4.0),
                 ifelse(tmp < -4.0, -4.0 + abs(tmp + 4.0), tmp))
    return(InvLogit(z))
}
#Function to build MME
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

# Change this to folder where cpp function was saved
Rcpp::sourceCpp('~/XXtEig.cpp')
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
# Change the folders where these are saved
cat('Accuracy', 'Replication', 'Poolsize_Pools', 'Error', 'Design', 'Scaling', 'Sires_Pool', 'Prog_sire',
    'Sire_rep', 'MeanRel', 'SDRel', 'ConstError' , sep = " ", 
    file = '~/SelSim1.txt', append = T)
cat('\n', file = '~/SelSim1.txt', append = T)
cat('Mean', 'SD', 'Min', 'Max', 'PoolSize1', 'Npooled', 'Replication', 'Pools', 'Error' , sep = " ", 
    file = '~/SelClus1.txt', append = T)
cat('\n', file = '~/SelClus1.txt', append = T)
cat('PoolSize', 'Error', 'Design', 'Scaling', 'ConstError', 'Replication', 'PopTBV', 
    'T10', 'T20', 'T40', 'B40', 'B20', 'B10', sep = ' ', 
    file = '~/TopBotSel.txt', append = T)
cat('\n', file = '~/TopBotSel.txt', append = T)
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
                sep = " ", file = '~/SelSim1.txt', append = T)
            cat('\n', file = '~/SelSim1.txt', append = T)
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
                sep = " ", file = '~/TopBotSel.txt', append = T)
            cat('\n', file = '~/TopBotSel.txt', append = T)
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
                sep = " ", file = '~/SelSim1.txt', append = T)
            cat('\n', file = '~/SelSim1.txt', append = T)
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
                sep = " ", file = '~/TopBotSel.txt', append = T)
            cat('\n', file = '~/TopBotSel.txt', append = T)
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
                sep = " ", file = '~/SelSim1.txt', append = T)
            cat('\n', file = '~/SelSim1.txt', append = T)
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
                sep = " ", file = '~/TopBotSel.txt', append = T)
            cat('\n', file = '~/TopBotSel.txt', append = T)
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
                sep = " ", file = '~/SelSim1.txt', append = T)
            cat('\n', file = '~/SelSim1.txt', append = T)
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
                sep = " ", file = '~/TopBotSel.txt', append = T)
            cat('\n', file = '~/TopBotSel.txt', append = T)
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
                sep = " ", file = '~/SelSim1.txt', append = T)
            cat('\n', file = '~/SelSim1.txt', append = T)
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
                sep = " ", file = '~/TopBotSel.txt', append = T)
            cat('\n', file = '~/TopBotSel.txt', append = T)
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
                sep = " ", file = '~/SelSim1.txt', append = T)
            cat('\n', file = '~/SelSim1.txt', append = T)
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
                sep = " ", file = '~/TopBotSel.txt', append = T)
            cat('\n', file = '~/TopBotSel.txt', append = T)
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
                sep = " ", file = '~/SelSim1.txt', append = T)
            cat('\n', file = '~/SelSim1.txt', append = T)
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
                sep = " ", file = '~/TopBotSel.txt', append = T)
            cat('\n', file = '~/TopBotSel.txt', append = T)
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
                sep = " ", file = '~/SelSim1.txt', append = T)
            cat('\n', file = '~/SelSim1.txt', append = T)
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
                sep = " ", file = '~/TopBotSel.txt', append = T)
            cat('\n', file = '~/TopBotSel.txt', append = T)
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
                sep = " ", file = '~/SelSim1.txt', append = T)
            cat('\n', file = '~/SelSim1.txt', append = T)
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
                sep = " ", file = '~/TopBotSel.txt', append = T)
            cat('\n', file = '~/TopBotSel.txt', append = T)
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
                sep = " ", file = '~/SelSim1.txt', append = T)
            cat('\n', file = '~/SelSim1.txt', append = T)
            cat(mean(m2$Size), sd(m2$Size), range(m2$Size), sum(m2$Size == 1), sum(m2$Size),
                irep, K, sigmavec[ss], sep = " ", 
                file = '~/SelClus1.txt', append = T)
            cat('\n', file = '~/SelClus1.txt', append = T)
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
                sep = " ", file = '~/TopBotSel.txt', append = T)
            cat('\n', file = '~/TopBotSel.txt', append = T)
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
                sep = " ", file = '~/SelSim1.txt', append = T)
            cat('\n', file = '~/SelSim1.txt', append = T)
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
                sep = " ", file = '~/TopBotSel.txt', append = T)
            cat('\n', file = '~/TopBotSel.txt', append = T)
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
                sep = " ", file = '~/SelSim1.txt', append = T)
            cat('\n', file = '~/SelSim1.txt', append = T)
            cat(mean(m2$Size), sd(m2$Size), range(m2$Size), sum(m2$Size == 1), sum(m2$Size),
                irep, K, sigmavec[ss], sep = " ", 
                file = '~/SelClus1.txt', append = T)
            cat('\n', file = '~/SelClus1.txt', append = T)
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
                sep = " ", file = '~/TopBotSel.txt', append = T)
            cat('\n', file = '~/TopBotSel.txt', append = T)
        }
    }
}
t2 = proc.time()
tt = t2 - t1; tt
