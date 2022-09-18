function [Rate_th_MR_unipl] = functionDownlinkClosedForm_unipl(MMSE_Cov_unipl,R_unipl,K,rho_unipl,prelogFactor_unipl)

%Theoretical uni-polarized rate calculation
part1MRuni_th=zeros(K,1);
part2MRuni_th=zeros(K,K);


for k=1:K
    
    
    part1MRuni_th(k)=rho_unipl(k)*abs(trace(MMSE_Cov_unipl(:,:,k)));
    
    for l=1:K
        
        
        part2MRuni_th(k,l)= rho_unipl(l)*trace(MMSE_Cov_unipl(:,:,l)*R_unipl(:,:,k))/abs(trace(MMSE_Cov_unipl(:,:,l)));
    end
    
    
    
    sumPart2MRuni_th=sum(part2MRuni_th(k,:));
    Rate_th_MR_unipl(k)=prelogFactor_unipl*log2(1 + part1MRuni_th(k)/(sumPart2MRuni_th +1));
end
