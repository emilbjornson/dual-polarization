function [W_MMSE,W_ZF,W_MR] = functionPrecoderMatrices(Hk_est,MMSEmatrixV,MMSEmatrixH,C_kV,C_kH,M,K,nbrOfRealizations,p1,p2,rho1,rho2)
%This file calculates the precoder matrices


[v_kV,v_kH, vZF_kV,vZF_kH,vMR_kV,vMR_kH] = functionCombiningVectors(Hk_est,C_kV,C_kH,M,K,nbrOfRealizations,p1,p2);


%Normalize the precoders
for nr=1:nbrOfRealizations
    for k=1:K
        
        denMMSEV(k,nr)=norm(v_kV(:,k,nr))^2;
        denZFV(k,nr)=norm(vZF_kV(:,k,nr))^2;
        
        denMMSEH(k,nr)=norm(v_kH(:,k,nr))^2;
        denZFH(k,nr)=norm(vZF_kH(:,k,nr))^2;
    end
end

denMMSEmeanV=  mean(denMMSEV,2);
denZFmeanV=  mean(denZFV,2);
denMMSEmeanH=  mean(denMMSEH,2);
denZFmeanH=  mean(denZFH,2);

for nr=1:nbrOfRealizations
    for k=1:K
        W_ZF(:,1,k,nr)= vZF_kV(:,k,nr)/sqrt(denZFmeanV(k));
        W_MMSE(:,1,k,nr)= v_kV(:,k,nr)/sqrt(denMMSEmeanV(k));
        
        W_ZF(:,2,k,nr)= vZF_kH(:,k,nr)/sqrt(denZFmeanH(k));
        W_MMSE(:,2,k,nr)= v_kH(:,k,nr)/sqrt(denMMSEmeanH(k));
        for na=1:2
            if na==1
                W_MR(:,na,k,nr)=vMR_kV(:,k,nr)/sqrt(trace(MMSEmatrixV(:,:,k)));
            else
                W_MR(:,na,k,nr)=vMR_kH(:,k,nr)/sqrt(trace(MMSEmatrixH(:,:,k)));
            end
        end
        
    end
end
% disp('Precoding is DONE.')



%Power allocation
for k=1:K
    Psqrt=diag([sqrt(rho1(k)),sqrt(rho2(k))]);
    for nr=1:nbrOfRealizations
        
        W_ZF(:,:,k,nr)=W_ZF(:,:,k,nr)*Psqrt;
        W_MR(:,:,k,nr)=W_MR(:,:,k,nr)*Psqrt;
        W_MMSE(:,:,k,nr)=W_MMSE(:,:,k,nr)*Psqrt;
    end
end
end

