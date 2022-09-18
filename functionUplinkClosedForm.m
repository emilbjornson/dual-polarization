function [Rate_th_MR] = functionUplinkClosedForm(MMSEmatrixV,MMSEmatrixH,Rkv,Rkh,K,pV_ul,pH_ul,prelogFactor)
%Calculate the closed form rate for MR as given in the paper

%Theoretical rate calculation
part1V_th=zeros(K,1);
part1H_th=zeros(K,1);
part2V_th=zeros(K,K);
part2H_th=zeros(K,K);



for k=1:K
    
    part1V_th(k)=pV_ul(k)*abs(trace(MMSEmatrixV(:,:,k)));
    part1H_th(k)=pH_ul(k)*abs(trace(MMSEmatrixH(:,:,k)));
    
    for l=1:K
        
        part2V_th(k,l)= pV_ul(l)*trace(MMSEmatrixV(:,:,k)*Rkv(:,:,l))/abs(trace(MMSEmatrixV(:,:,k))) + pH_ul(l)*trace(MMSEmatrixV(:,:,k)*Rkh(:,:,l))/abs(trace(MMSEmatrixV(:,:,k)));
        
        part2H_th(k,l)= pH_ul(l)*trace(MMSEmatrixH(:,:,k)*Rkh(:,:,l))/abs(trace(MMSEmatrixH(:,:,k))) + pV_ul(l)*trace(MMSEmatrixH(:,:,k)*Rkv(:,:,l))/abs(trace(MMSEmatrixH(:,:,k)));
        
    end
    
    sumPart2V_th=sum(part2V_th(k,:));
    sumPart2H_th=sum(part2H_th(k,:));
    
    Rate_th_MR(k)=prelogFactor*log2(1 + part1V_th(k)/(sumPart2V_th +1)) + prelogFactor*log2(1 + part1H_th(k)/(sumPart2H_th + 1));
    
end
end

