function [Rate_SIC] = functionDownlinkMMSESICbound(H,W,K,nbrOfRealizations,prelogFactor)
%Downlink MMSE-SIC bound calculation

Term1=zeros(2,2,K,nbrOfRealizations);
Term2=zeros(2,2,K,K,nbrOfRealizations);


%Monte-Carlo rate calculation
for k=1:K
    for nr=1:nbrOfRealizations
        Term1(:,:,k,nr) = H(:,:,nr,k)*W(:,:,k,nr);
        for l=1:K
            Term2(:,:,k,l,nr)= H(:,:,nr,k)*W(:,:,l,nr)*W(:,:,l,nr)'*H(:,:,nr,k)';  
        end
    end
end

meanTerm1=mean(Term1,4);
meanTerm2=sum(mean(Term2,5),4);
Omega = zeros(2,2,K);



for k=1:K
    %MRT
    Omega(:,:,k)=(meanTerm2(:,:,k)  - meanTerm1(:,:,k)* meanTerm1(:,:,k)' + eye(2));
    Rate_SIC(k) = prelogFactor*log2(det( eye(2) + (meanTerm1(:,:,k)'/Omega(:,:,k))*meanTerm1(:,:,k)  ));
end
end

