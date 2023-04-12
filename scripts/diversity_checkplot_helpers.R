# helper functions for generating diversity checkplots
#MR function to take abundance vector and "l" and return the quantile of B bootstrap iterations of Chao technique that true value falls on. Could be extended to include sample Hill as Chao seems to have
# x<-1:5
# B<-10
# l<-1
# truediv<-5

checkchao <- function(x, B, l, truediv, conf = 0.95){ #, truemu_n
    n<-sum(x)
    #columns of this matrix are replicate boostraps
    data.bt = rmultinom(B,n,Bt_prob_abu(x))
    #get estimator for each bootstrapped sample
    pro = apply(data.bt,2,function(boot)Chao_Hill_abu(boot, l))
    #mean correction
    pro<-pro-mean(pro)+Chao_Hill_abu(x, l)
    
    #break ties
    less<-sum(pro<truediv)/length(pro)
    more<-(length(pro)-sum(pro>truediv))/length(pro)
    p<-runif(1, min(less, more), max(less, more))
    
    lower<-max(pro[which(min_rank(pro)<=max(floor(B*(1-conf)/2),1))])
    upper<-min(pro[which(min_rank(-pro)<=max(floor(B*(1-conf)/2),1))])
    
    return(data.frame(p=p
                      , lower=lower
                      , upper=upper
                      , truediv=truediv
                      , "chaoest"=Chao_Hill_abu(x, l)
                      , "obsD"=rarity(x,l)
    )
    
    )}

################
#truemu computes the emprical average sample diveristy under sampling without replacement
truemu<-function(comm, size, reps, l,...){
    sam<-replicate(reps, subsam(comm, size))
    return(
        mean(
            apply(sam,2, function(x){dfun(x, l)})
        )
    )
}


truemu_inf<-function(comm, size, reps, l,...){ #comm is abundance vector; size, reps, l all constants
    sam<-replicate(reps, sample_infinite(comm, size=size))
    
    return(
        mean(
            apply(sam,2, function(x){dfun(x, l)})
        )
    )
}


checkplot_inf <- function(SAD, B = 999, l, inds, reps){
    hillname <- ifelse(l == -1
                       , "Hill-Simpson"
                       , ifelse(l == 0
                                , "Hill-Shannon"
                                , "richness"))
    td <- SAD$community_info[hillname] #grab true diversity from SAD object
    
    future_map_dfr(1:reps, function(x){
        
        
        obs <- MeanRarity::sample_infinite(SAD$rel_abundances, size = inds) #subsample the whole community with # individuals=size
        
        chaotile <- checkchao(x = obs, B = B, l = l, truediv = td) #then do B bootstrap samples for the augmented community based on that sample
        return(myout = data.frame(p = chaotile$p
                                  , truediv = chaotile$truediv
                                  , chaoest = chaotile$chaoest
                                  , obsD = chaotile$obsD
                                  , upper = chaotile$upper
                                  , lower = chaotile$lower
                                  , l
                                  , inds
                                  , reps)
        )
        
    })
}
