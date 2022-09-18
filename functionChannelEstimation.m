function [vecHk_est,MMSEmatrixV, MMSEmatrixH,Rkv, Rkh] = functionChannelEstimation(H,Rk_sqrtm,q_XPD,M,K,nbrOfRealizations,tau_p,p1,p2)
%MMSE channel estimation function

L_k_sqrtm=diag([sqrt(p1), sqrt(p2)]);
%Generate realizations of normalized noise
Np = sqrt(0.5)*(randn(M,2,nbrOfRealizations,K) + 1i*randn(M,2,nbrOfRealizations,K));
Y_pk=zeros(M,2,nbrOfRealizations,K);
A=kron(tau_p*L_k_sqrtm,eye(M));
%Processed pilot signal Mx2
for k=1:K
    for nr=1:nbrOfRealizations
        Y_pk(:,:,nr,k)=tau_p*H(:,:,nr,k)'*L_k_sqrtm + sqrt(tau_p)*Np(:,:,nr,k);
    end
end


%MMSE channel estimation
diag1=repmat([(1-q_XPD), q_XPD],1,M/2);
diag2=repmat([q_XPD, (1-q_XPD) ],1,M/2);
Rkv=zeros(M,M,K);
Rkh=zeros(M,M,K);
R_bk=zeros(2*M,2*M,K);
C_MSE=zeros(K,1);
Cov_MMSE=zeros(K,1);
vecHk_est=zeros(2*M,nbrOfRealizations,K);
test1_normhest=zeros(K,1);
test2_normhest=zeros(K,1);
MMSEmatrixV= zeros(M,M,K);
MMSEmatrixH= zeros(M,M,K);

for k=1:K
    
    %Create the covariance matrices
    Rkv(:,:,k)= Rk_sqrtm(:,:,k)*diag(diag1)*Rk_sqrtm(:,:,k)';
    Rkh(:,:,k)= Rk_sqrtm(:,:,k)*diag(diag2)*Rk_sqrtm(:,:,k)';
    R_bk(:,:,k) = blkdiag(Rkv(:,:,k),Rkh(:,:,k));
    
    
    C_MSE(k) = trace(R_bk(:,:,k) - (R_bk(:,:,k)*A')/(A*R_bk(:,:,k)*A' + (tau_p)*eye(2*M))*(A*R_bk(:,:,k)));
    Cov_MMSE(k) = trace((R_bk(:,:,k)*A')/(A*R_bk(:,:,k)*A' + (tau_p)*eye(2*M))*(A*R_bk(:,:,k)));
    
    
    test1_normhest(k)=trace((p1*tau_p*Rkv(:,:,k)/(p1*tau_p*Rkv(:,:,k)+eye(M)))*Rkv(:,:,k));
    test2_normhest(k)=trace((p2*tau_p*Rkh(:,:,k)/(p2*tau_p*Rkh(:,:,k)+eye(M)))*Rkh(:,:,k));
    
    
    MMSEmatrixV(:,:,k)=(p1*tau_p*Rkv(:,:,k)/(p1*tau_p*Rkv(:,:,k)+eye(M)))*Rkv(:,:,k);
    MMSEmatrixH(:,:,k)=(p2*tau_p*Rkh(:,:,k)/(p2*tau_p*Rkh(:,:,k)+eye(M)))*Rkh(:,:,k);
    
    %MMSE estimator
    for nr=1:nbrOfRealizations
        
        vecY_pk=Y_pk(:,:,nr,k);
        vecHk_est(:,nr,k)=(R_bk(:,:,k)*A')/(A*R_bk(:,:,k)*A' + (tau_p)*eye(2*M))*vecY_pk(:);
        
    end
    
end
end

