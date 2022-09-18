function [H,Rk_sqrtm] = functionChannelGeneration(Sigma,R_BSk, C_BS,C_UE,F_BS,F_UE,M,K,nbrOfRealizations)
% This file generates the propagation channel

    Mdual=M/2;
    first_term=kron(ones(1,Mdual),Sigma);
    Rk=zeros(M,M,K);
    Rk_sqrtm=zeros(M,M,K);
    R_UE_sqrtm=zeros(2,2,K);
    %the propagation channel
    Z=zeros(2,M,nbrOfRealizations,K);
    M_BS = kron(eye(Mdual),F_BS);
    H=zeros(2,M,nbrOfRealizations,K);

    for k=1:K

        Rk(:,:,k)=kron(R_BSk(:,:,k),C_BS);
        Rk_sqrtm(:,:,k)=sqrtm(Rk(:,:,k));
        R_UE_sqrtm(:,:,k)=sqrtm(C_UE); % R_UE=1 for co-located DP at UE

        for nr=1:nbrOfRealizations

            S=sqrt(0.5)*(randn(2,M) + 1i*randn(2,M));
            %This is the propagation channel (without the XPI effects)
            Z(:,:,nr,k)=first_term.*(R_UE_sqrtm(:,:,k)*S*Rk_sqrtm(:,:,k));
            %The channel including antenna effects at both BS and
            %UE sides
            H(:,:,nr,k)=F_UE(:,:,k)*Z(:,:,nr,k)*M_BS;

        end
    end
end

