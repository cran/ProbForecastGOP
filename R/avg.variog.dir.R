"avg.variog.dir" <-
function(day,coord1,coord2,id,gop.res,cut.points,tol.angle.rad1,tol.angle.rad2,type){
   unique.day <- unique(day)
   l.day <- length(unique.day)
   n.stations <- length(unique(id))
   n.cuts <- length(cut.points)
   W <- rep(0,(n.cuts-1))
   W.obs <- rep(0,(n.cuts-1))
   for(i in 1:l.day){
       distance.day <- NULL
       difference.day <- NULL
       new.dist.day <- NULL
       s.new.dist.day <- NULL
       o.new.dist.day <- NULL
       new.diff.day <- NULL
       s.new.diff.day <- NULL
       gop.res.day <- NULL
       id.day <- NULL
       coord1.day <- NULL
       coord2.day <- NULL
       gop.res.day <- gop.res[day==unique.day[i]]
       id.day <- id[day==unique.day[i]]
       coord1.day <- coord1[day==unique.day[i]]
       coord2.day <- coord2[day==unique.day[i]]
       distance.day <- 
calc.dist.dir(coord1.day,coord2.day,id.day,tol.angle.rad1,tol.angle.rad2,type)
       new.dist.day <- distance.day[lower.tri(distance.day)]
       s.new.dist.day <- sort(new.dist.day)
       o.new.dist.day <- order(new.dist.day)
       difference.day  <- calc.difference(gop.res.day)      
       new.diff.day <- difference.day[lower.tri(difference.day)] 
       s.new.diff.day <- new.diff.day[o.new.dist.day]
       W.new <- rep(0,(n.cuts-1))
       W.obs.new <- rep(0,(n.cuts-1))
   
       for(s in 1:(n.cuts-1)){
           low.bound <- cut.points[s]
           upp.bound <- cut.points[s+1]
           v.dist <- NULL
           v.dist1 <- NULL
           v.diff1 <- NULL
           index.v.dist <- NULL
           index.v.dist1 <- NULL

           if(s < (n.cuts-1)){
                v.dist <- s.new.dist.day[s.new.dist.day >= low.bound & 
s.new.dist.day < upp.bound]
                index.v.dist <- 
seq(1:length(s.new.dist.day))[s.new.dist.day >= low.bound & s.new.dist.day 
< upp.bound]
                v.dist1 <- v.dist[v.dist!=0]
                index.v.dist1 <- index.v.dist[v.dist!=0]

                if(length(v.dist1) >=1){  
                   v.diff1 <- s.new.diff.day[index.v.dist1]
                   W.new[s] <- sum(v.diff1)
                   W.obs.new[s] <- length(v.dist1)}
                if(length(v.dist1) ==0){  
                   W.new[s] <- 0
                   W.obs.new[s] <- 0}
           }
           if(s==(n.cuts-1)){
                 v.dist <- s.new.dist.day[s.new.dist.day >= low.bound & 
s.new.dist.day <= upp.bound]
                 index.v.dist <- 
seq(1:length(s.new.dist.day))[s.new.dist.day >= 
low.bound & s.new.dist.day <= upp.bound]
                 v.dist1 <- v.dist[v.dist!=0]
                 index.v.dist1 <- index.v.dist[v.dist!=0]

                 if(length(v.dist1) >=1){  
                   v.diff1 <- s.new.diff.day[index.v.dist1]
                   W.new[s] <- sum(v.diff1)
                   W.obs.new[s] <- length(v.dist1)}
                 if(length(v.dist1) ==0){  
                   W.new[s] <- 0
                   W.obs.new[s] <- 0}
           }
       }           
   W <- W+W.new
   W.obs <- W.obs + W.obs.new
   }         
   avg.variog <- round(W/(2*W.obs),2)  
   n.h <- W.obs
   x.vec <- NULL
   for(i in 1:(n.cuts-1)){
     x.vec[i] <- (cut.points[i]+cut.points[i+1])/2
   }
   fin.avg.variog <- c(0,avg.variog)   
   fin.x.vec <- c(0,x.vec)
   fin.n.h <- c(n.stations,n.h)
   B <- cbind(fin.x.vec,fin.n.h,fin.avg.variog)
   return(B)
}
