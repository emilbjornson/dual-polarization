function [v_kV,v_kH, vZF_kV,vZF_kH,vMR_kV,vMR_kH] = functionCombiningVectors(Hk_est,C_kV,C_kH,M,K,nbrOfRealizations,pV_ul,pH_ul)
%This function calculate combining vectors

HPH=zeros(M,M,K);
Hestall=zeros(M,2*K,nbrOfRealizations);
W_ZF_all=zeros(M,2*K,nbrOfRealizations);


for nr=1:nbrOfRealizations
    
    for k=1:K
        HPH(:,:,k)=Hk_est(:,:,nr,k)*[pV_ul(k) 0; 0 pH_ul(k)]*Hk_est(:,:,nr,k)';
        UpsilonSICx(:,:,k)=pV_ul(k)*C_kV(:,:,k)+ pH_ul(k)*C_kH(:,:,k) ;
    end
    
    Upsilon = sum(HPH,3)+  sum(UpsilonSICx,3) +eye(M);
    invUpsilon=inv(Upsilon);
    
    Hestall(:,:,nr)=reshape(Hk_est(:,:,nr,:),M,2*K);
    W_ZF_all(:,:,nr) = Hestall(:,:,nr)/(Hestall(:,:,nr)'*Hestall(:,:,nr));
    W_ZF=reshape(W_ZF_all,M,2,K,nbrOfRealizations);
    
    for k=1:K
        
        %MMSE Filter
        v_kV(:,k,nr)=sqrt(pV_ul(k))*invUpsilon*Hk_est(:,1,nr,k);
        v_kH(:,k,nr)=sqrt(pH_ul(k))*invUpsilon*Hk_est(:,2,nr,k);
        
        %ZF Combining
        vZF_kV(:,k,nr)=W_ZF(:,1,k,nr);
        vZF_kH(:,k,nr)=W_ZF(:,2,k,nr);
        
        %MR Combining
        vMR_kV(:,k,nr)=Hk_est(:,1,nr,k);
        vMR_kH(:,k,nr)=Hk_est(:,2,nr,k);
        
    end
end
end

