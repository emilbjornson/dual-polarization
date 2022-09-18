function [W_MMSE_unipl,W_ZF_unipl,W_MR_unipl] = functionPrecoders_unipl(hk_est_unipl,C_unipl,MMSE_Cov_unipl,Muni,K,nbrOfRealizations,p_unipl,rho_unipl)
%Precoding Part of uni-polarized antennas
Hestall_unipl=zeros(Muni,K,nbrOfRealizations);
W_ZF_unipl=zeros(Muni,K,nbrOfRealizations);

for nr=1:nbrOfRealizations
    
    %For MMSE precoder
    for k=1:K
        phhH(:,:,k)=p_unipl*hk_est_unipl(:,nr,k)*hk_est_unipl(:,nr,k)';
    end
    
    Upsilon_unipl = sum(phhH,3)+  p_unipl*sum(C_unipl,3) + eye(Muni);
    invUpsilon_unipl=inv(Upsilon_unipl);
    
    %For ZF precoder
    Hestall_unipl(:,:,nr)=reshape(hk_est_unipl(:,nr,:),Muni,K);
    W_ZF_unipl(:,:,nr) = Hestall_unipl(:,:,nr)/(Hestall_unipl(:,:,nr)'*Hestall_unipl(:,:,nr));
    for k=1:K
        
        vMMSE_unipl(:,k,nr)=sqrt(p_unipl)*invUpsilon_unipl*hk_est_unipl(:,nr,k);
        
        denMMSEunipl(k,nr)=norm(vMMSE_unipl(:,k,nr))^2 ;
        denZFunipl(k,nr)=norm(W_ZF_unipl(:,k,nr))^2 ;
        
    end
end

denMMSEmean=  mean(denMMSEunipl,2);
denZFmean=  mean(denZFunipl,2);

W_MMSE_unipl=zeros(Muni,K,nbrOfRealizations);
W_MR_unipl=zeros(Muni,K,nbrOfRealizations);

%Normalize
for nr=1:nbrOfRealizations
    for k=1:K
        
        W_MMSE_unipl(:,k,nr)= sqrt(rho_unipl(k))*vMMSE_unipl(:,k,nr)/sqrt(denMMSEmean(k));
        W_ZF_unipl(:,k,nr)= sqrt(rho_unipl(k))*W_ZF_unipl(:,k,nr)/sqrt(denZFmean(k));
        W_MR_unipl(:,k,nr)=sqrt(rho_unipl(k))*hk_est_unipl(:,nr,k)/sqrt(trace(MMSE_Cov_unipl(:,:,k)));
    end
end
% disp('MR and ZF are DONE.')
end

