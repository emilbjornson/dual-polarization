function [RateSICnr] = functionUplinkMMSESICbound(Hk_est,C_kV,C_kH,M,K,pV_ul,pH_ul,nbrOfRealizations,prelogFactor)
%Calculate MMSE-SIC Rate (for 2K stream decoding)

for nr=1:nbrOfRealizations
    
    
    for k=1:K
        UpsilonSICx(:,:,k)=pV_ul(k)*C_kV(:,:,k)+ pH_ul(k)*C_kH(:,:,k) ;
    end
    invUpsilonSIC=inv(sum(UpsilonSICx,3)  +eye(M));
    
    
    SNR_SIC=zeros(M,M,K);
    for l=1:K
        
        SNR_SIC(:,:,l)=Hk_est(:,:,nr,l)*[pV_ul(l) 0; 0 pH_ul(l)]*Hk_est(:,:,nr,l)'*invUpsilonSIC;
        
    end
    
    RateSICnr (nr) = prelogFactor*abs(log2(det(eye(M) + sum(SNR_SIC,3))));
    
end

end

