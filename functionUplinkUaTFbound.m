function [Rate_UaTF] = functionUplinkUaTFbound(H,v_kV,v_kH,K,nbrOfRealizations,pV_ul,pH_ul,prelogFactor)
%This file calculates the uplink UaTF bound for given combining vectors

Term1V=zeros(K,nbrOfRealizations);
Term1H=zeros(K,nbrOfRealizations);
Term3V=zeros(K,nbrOfRealizations);
Term3H=zeros(K,nbrOfRealizations);

Term2kVV=zeros(K,K);
Term2kVH=zeros(K,K);
Term2kHV=zeros(K,K);
Term2kHH=zeros(K,K);
sumTerm2V=zeros(K,nbrOfRealizations);
sumTerm2H=zeros(K,nbrOfRealizations);


for nr=1:nbrOfRealizations
    
    for k=1:K
        
        
        Term1V(k,nr)= sqrt(pV_ul(k))*v_kV(:,k,nr)'*H(1,:,nr,k)';
        Term1H(k,nr)= sqrt(pH_ul(k))*v_kH(:,k,nr)'*H(2,:,nr,k)';
        
        Term3V(k,nr)= v_kV(:,k,nr)'*v_kV(:,k,nr);
        Term3H(k,nr)= v_kH(:,k,nr)'*v_kH(:,k,nr);
        
        for l=1:K
            %MMSE Part
            Term2kVV(k,l)= pV_ul(l)*abs(v_kV(:,k,nr)'*H(1,:,nr,l)')^2;
            Term2kVH(k,l)= pH_ul(l)*abs(v_kV(:,k,nr)'*H(2,:,nr,l)')^2;
            
            Term2kHV(k,l)= pV_ul(l)*abs(v_kH(:,k,nr)'*H(1,:,nr,l)')^2;
            Term2kHH(k,l)= pH_ul(l)*abs(v_kH(:,k,nr)'*H(2,:,nr,l)')^2;
            
            
        end
    end
    
    %MMSE part
    sumTerm2V(:,nr)=sum(Term2kVV,2) + sum(Term2kVH,2);
    sumTerm2H(:,nr)=sum(Term2kHV,2) + sum(Term2kHH,2);
   
end

%MMSE PART
meanTerm1V=mean(Term1V,2);
meanTerm3V=mean(Term3V,2);
meanTerm2V=mean(sumTerm2V,2);


meanTerm1H=mean(Term1H,2);
meanTerm3H=mean(Term3H,2);
meanTerm2H=mean(sumTerm2H,2);

%Calculate the rates
for k=1:K
    
    %v_kV MMSE
    rate_kV(k)=log2(1 + abs(meanTerm1V(k))^2 /(meanTerm2V(k) - abs(meanTerm1V(k))^2 + meanTerm3V(k) )  );
    rate_kH(k)=log2(1 + abs(meanTerm1H(k))^2 /(meanTerm2H(k) - abs(meanTerm1H(k))^2 + meanTerm3H(k) )  );
    Rate_UaTF(k)=prelogFactor*abs(rate_kV(k)+ rate_kH(k));

end
end

