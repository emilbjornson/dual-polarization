function [Rate_UatF] = functionDownlinkMMSESICbound_unipl(Hunipl,W,K,nbrOfRealizations,prelogFactor_unipl)

Term1_unipl=zeros(K,nbrOfRealizations);
Term2_unipl=zeros(K,K,nbrOfRealizations);


%The rate calculation Monte-Carlo
for k=1:K
    for nr=1:nbrOfRealizations
        
        
        Term1_unipl(k,nr)=Hunipl(:,nr,k)'*W(:,k,nr);
        
        
        for l=1:K
            
            Term2_unipl(k,l,nr)=Hunipl(:,nr,k)'*W(:,l,nr)*W(:,l,nr)'*Hunipl(:,nr,k);
            
        end
        
        
    end
    
end



meanTerm1_unipl=mean(Term1_unipl,2);
meanTerm2_unipl=sum(mean(Term2_unipl,3),2);
Omega_unipl = zeros(K,1);


for k=1:K
    
    
    
    %MRT uni-polarized antennas
    Omega_unipl(k)=(meanTerm2_unipl(k)  - abs(meanTerm1_unipl(k))^2 + 1);
    Rate_UatF(k) = prelogFactor_unipl*log2(1+ (abs(meanTerm1_unipl(k))^2)/Omega_unipl(k));
    
    
end

end

