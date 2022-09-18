function [Rate_woSIC] = functionDownlinkwithoutSICbound(H,W,K,nbrOfRealizations,prelogFactor)
%The rate calculation (without SIC bound)

Term1=zeros(2,2,K,nbrOfRealizations);
Term2=zeros(2,2,K,K,nbrOfRealizations);


%The rate calculation Monte-Carlo
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

for k=1:K
    
    OmegaBar(:,:,k)= meanTerm2(:,:,k) + eye(2);
    OmegaBarInv(:,:,k)=inv(OmegaBar(:,:,k));
    
    vDL_kV(:,k)=OmegaBarInv(:,:,k)*meanTerm1(:,1,k);
    vDL_kH(:,k)=OmegaBarInv(:,:,k)*meanTerm1(:,2,k);
    

    
    Rate_woSIC_kV(k)=prelogFactor*abs(log2(1 +(abs(vDL_kV(:,k)'*meanTerm1(:,1,k))^2)/...
        (vDL_kV(:,k)'*(OmegaBar(:,:,k)-meanTerm1(:,1,k)*meanTerm1(:,1,k)')*vDL_kV(:,k) )));
    
    Rate_woSIC_kH(k)=prelogFactor*abs(log2(1 +(abs(vDL_kH(:,k)'*meanTerm1(:,2,k))^2)/...
        (vDL_kH(:,k)'*(OmegaBar(:,:,k)-meanTerm1(:,2,k)*meanTerm1(:,2,k)')*vDL_kH(:,k) )));
    
    Rate_woSIC(k)= Rate_woSIC_kV(k)+   Rate_woSIC_kH(k);
    
    
end
end

