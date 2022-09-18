function [Rate_th_MR] = functionDownlinkClosedForm(MMSEmatrixV,MMSEmatrixH,Rkv,Rkh,K,rho1,rho2,prelogFactor)
%Theoretical rate calculation as given in the paper

part1V_th=zeros(K,1);
part1H_th=zeros(K,1);
part2V_th=zeros(K,K);
part2H_th=zeros(K,K);

for k=1:K
    
    part1V_th(k)=rho1(k)*abs(trace(MMSEmatrixV(:,:,k)));
    part1H_th(k)=rho2(k)*abs(trace(MMSEmatrixH(:,:,k)));
    
    for l=1:K
        
        part2V_th(k,l)= rho1(l)*trace(MMSEmatrixV(:,:,l)*Rkv(:,:,k))/abs(trace(MMSEmatrixV(:,:,l))) + rho2(l)*trace(MMSEmatrixH(:,:,l)*Rkv(:,:,k))/abs(trace(MMSEmatrixH(:,:,l)));
        
        part2H_th(k,l)= rho2(l)*trace(MMSEmatrixH(:,:,l)*Rkh(:,:,k))/abs(trace(MMSEmatrixH(:,:,l))) + rho1(l)*trace(MMSEmatrixV(:,:,l)*Rkh(:,:,k))/abs(trace(MMSEmatrixV(:,:,l)));
        
    end
    
    sumPart2V_th=sum(part2V_th(k,:));
    sumPart2H_th=sum(part2H_th(k,:));
    
    Rate_th_MR(k)=prelogFactor*log2(1 + part1V_th(k)/(sumPart2V_th +1)) + prelogFactor*log2(1 + part1H_th(k)/(sumPart2H_th + 1));
    
end
end

