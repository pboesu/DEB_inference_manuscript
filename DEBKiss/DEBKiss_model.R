## DEBKiss Process Model
library(deSolve)

DEBKiss1<-function(time,y,params){

    with(as.list(c(y, params)),{

        ec<-0.000001 ## value of WB below which hatching occurs, this
                     ## is >0 for numerical reasons

        ## physical length, Lw, is the input
        L<-Lw*deltaM ## structural length        
        L3<-L^3 ## structural volume         

        ## calculating the length at puberty from the weight at
        ## puberty
        if(!is.null(Wp)) Lp<-(Wp/dV)^(1/3)
        
        JA<-f*JAM*L^2 ## assimilation
        JX<-JA/yAX ## feeding

        ## first check if still egg stage and modify if needed
        if(WB>ec){
            JA <- JA*fb/f ## this sets f to fb, if needed
            JX <- 0
            dWB <- - JA
        }else dWB<-0
        
        JM<-exp(logJMv)*L3 ## volume maintance costs
        JV<-yVA*(kappa*JA-JM) ## growth
        JR<- (1-kappa)*JA ## repo buffer
                

        ## Starvation conditions, after hatching
        if(WB<=ec){
            if(JA>JM){
                if( JV<0 ){
                    JV<-0
                    JR<-JA-JM
                }
            }else{
                JV<-(JA-JM)/yAV
                JR<-0
            }
        }

        if(L<Lp) JR<-0 ## check for being mature

       
        dLw <- JV/(3*dV*L^2*deltaM) ## change in length
        dWR <- JR ## change in reproductive buffer


        list(c(dWB,dLw,dWR))
    })

}

