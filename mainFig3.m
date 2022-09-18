% This Matlab script generates Figure 3 in the paper:
%
%O. Ozdogan and E. Bjornson, "Massive MIMO with Dual-Polarized Antennas,"
%in IEEE Transactions on Wireless Communications, 2022,
%doi: 10.1109/TWC.2022.3205471
%
% This is version 1.0 (Last edited: 2021-11-11)
%
% License: This code is licensed under the GPLv2 license. If you in any way
% use this code for research that results in publications, please cite our
% paper as described above.

%M-unipolarized antennas
%Comparison of dual-polarized and uni-polarized setups
%Average downlink sum SE for 10 UEs with different preocoders
%as a function of the number of BS antennas for dual-polarized and
%uni-polarized setups
close all;clear;clc;

%Number of users
K=10;
L=1; %number of cells (keep this as it is)

%Communication bandwidth
B = 20e6;
%Noise figure at the BS (in dB)
noiseFigure = 7;
%Compute noise power
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;
%Select length of coherence block
tau_c = 200;


%Spatial correlation matrix parameters
ASDdeg = 5;
antennaSpacing = 1/2; %d_H

%%Create polarization-related matrices
%The deterministic polarization matrix of the BS
F_BS= eye(2);
%The polarization matrix of the UE
F_UE=reshape(repmat(eye(2),1,K),2,2,K);



%Polarization correlation matrices
r_p = 0;
t_p = 0;
C_UE = [1, r_p; r_p', 1];
C_BS = [1, t_p; t_p', 1];

%XPD: polar-discrimination (same for all users)
%0<= q_XPD <= 1
%Small  q_XPD means LoS-like channels 
XPDdB = 5;
q_XPD= 1/(1 + db2pow(XPDdB)); %keep 1 XPD value for the uni-polarized
%(uni-polarized setups does not depend on XPD)

%Prepare for the simulations
nbrOfRealizations=200;
nbrOfSetups=50;
%number of maximum antennas
Mmax=100;
Mrange=40:10:Mmax;

%Pilot transmit powers
p1=100; %mW
p2=100; %mW
pV_ul=100*ones(K,1);
pH_ul=100*ones(K,1);

%Unipolarized 
p_unipl=p1;%DO NOT FORGET TO CHANGE THIS FOR M/2 (2*p1) and M (p1)

%Dowlink transmit power
rho1=100*ones(K,1);
rho2=100*ones(K,1);
rho_unipl=2*rho1; %DO NOT FORGET TO CHANGE THIS FOR M/2 and M

%pilot length 
tau_p=2*K;
tau_p_unipl=K;
prelogFactor= (tau_c-tau_p)/tau_c;
prelogFactor_unipl= (tau_c-tau_p_unipl)/tau_c;

%Prepare to store the rates
Rate_th_MR=zeros(length(Mrange),length(q_XPD),nbrOfSetups,K);
Rate_th_MR_unipl=zeros(length(Mrange),length(q_XPD),nbrOfSetups,K);


Rate_SIC_MMSE= zeros(length(Mrange),length(q_XPD),nbrOfSetups,K);
Rate_SIC_ZF= zeros(length(Mrange),length(q_XPD),nbrOfSetups,K);
Rate_SIC_MR= zeros(length(Mrange),length(q_XPD),nbrOfSetups,K);

Rate_UatF_MMSEunipl= zeros(length(Mrange),length(q_XPD),nbrOfSetups,K);
Rate_UatF_ZFunipl= zeros(length(Mrange),length(q_XPD),nbrOfSetups,K);
Rate_UatF_MRunipl= zeros(length(Mrange),length(q_XPD),nbrOfSetups,K);

for n=1:nbrOfSetups
    
    %Compute spatial correlation matrices using the local scattering model
    [R_SS_normalized,R_SS_UE,channelGaindB] = functionExampleSetup(L,K,Mmax,ASDdeg);
    
    %Compute the normalized average channel gain, where the normalization
    %is based on the noise power
    channelGainOverNoise = channelGaindB - noiseVariancedBm;
    
    %Controlled UL channel gain over noise
    %channelGainOverNoise= functionPowerControl( channelGaindB,noiseVariancedBm,deltadB,L);
    betas=10.^(channelGainOverNoise./10);
    
    
    for m=1:length(Mrange)
        
        M=Mrange(m);
        Mdual=M/2;
        
        Muni=M; %Compare with unipolarized antennas, M or M/2
        
        R_unipl=zeros(Muni,Muni,K);
        R_BSk=zeros(Mdual,Mdual,K);
        for k=1:K
            R_unipl(:,:,k)= R_SS_normalized(1:Muni,1:Muni,k)*betas(k);
            R_BSk(:,:,k)=R_SS_normalized(1:Mdual,1:Mdual,k)*betas(k);
        end
        
        
        for xpd=1:length(q_XPD)
            
            %Generating the propagation channel
            %Dual-polarized antennas
            Sigma = [sqrt(1-q_XPD(xpd)), sqrt(q_XPD(xpd)); sqrt(q_XPD(xpd)), sqrt(1-q_XPD(xpd))];
            
            [H,Rk_sqrtm] = functionChannelGeneration(Sigma,R_BSk, C_BS,C_UE,F_BS,F_UE,M,K,nbrOfRealizations);
            
            %Unipolarized antennas
            Rk_unipl_sqrtm=zeros(Muni,Muni,K);
            Hunipl=zeros(Muni,nbrOfRealizations,K);
            
            for k=1:K
                
                Rk_unipl_sqrtm(:,:,k)=sqrtm(R_unipl(:,:,k));
                
                for nr=1:nbrOfRealizations
                    
                    %uni-polarized antennas
                    %both BS and UE has uni-polarized antennas
                    S_unipl=sqrt(0.5)*(randn(Muni,1) + 1i*randn(Muni,1));
                    Hunipl(:,nr,k)=Rk_unipl_sqrtm(:,:,k)*S_unipl;
                    
                end
            end
            
            
            
            
            %Channel estimation part
            %Dual-polarized antennas
            [vecHk_est,MMSEmatrixV, MMSEmatrixH,Rkv, Rkh] = functionChannelEstimation(H,Rk_sqrtm,q_XPD(xpd),M,K,nbrOfRealizations,tau_p,p1,p2);
            %Here are the estimated channels
            %Channel Estimation Error Matrices
            C_kV=Rkv-MMSEmatrixV;
            C_kH=Rkh-MMSEmatrixH;
            Hk_est=reshape(vecHk_est, M,2,nbrOfRealizations,K);
            
            %channel estimation - unipolarized antennas
            Np_unipl = sqrt(0.5)*(randn(Muni,nbrOfRealizations,K) + 1i*randn(Muni,nbrOfRealizations,K));
            Y_pk_unipl=zeros(Muni,nbrOfRealizations,K);
            %Received Pilot Signal
            for k=1:K
                
                for nr=1:nbrOfRealizations
                    %processed pilot signal
                    Y_pk_unipl(:,nr,k)=sqrt(p_unipl)*tau_p_unipl*Hunipl(:,nr,k) + sqrt(tau_p_unipl)*Np_unipl(:,nr,k);
                end
            end
            
            
            %MMSE estimation
            hk_est_unipl=zeros(Muni,nbrOfRealizations,K);
            MMSE_Cov_unipl=zeros(Muni,Muni,K);
            %Define other covariance matrices
            for k=1:K
                
                %Covariance matrix of uni-polarized antennas
                MMSE_Cov_unipl(:,:,k)=(p_unipl*tau_p_unipl*R_unipl(:,:,k)/(p_unipl*tau_p_unipl*R_unipl(:,:,k) + eye(Muni)))*R_unipl(:,:,k);
                %MMSE estimator
                for nr=1:nbrOfRealizations
                    
                    hk_est_unipl(:,nr,k)=sqrt(p_unipl)*R_unipl(:,:,k)/(p_unipl*tau_p_unipl*R_unipl(:,:,k) + eye(Muni))*Y_pk_unipl(:,nr,k);
                    
                end
                
            end
            
            C_unipl=R_unipl-MMSE_Cov_unipl;
            % disp('Channel Estimation is DONE.')
            
            
            
            
            %Precoding part
            %Dual-polarized
            [W_MMSE,W_ZF,W_MR] = functionPrecoderMatrices(Hk_est,MMSEmatrixV,MMSEmatrixH,C_kV,C_kH,M,K,nbrOfRealizations,pV_ul,pH_ul,rho1,rho2);
            
            %Unipolarized
            [W_MMSE_unipl,W_ZF_unipl,W_MR_unipl] = functionPrecoders_unipl(hk_est_unipl,C_unipl,MMSE_Cov_unipl,Muni,K,nbrOfRealizations,p_unipl,rho_unipl);
            

            
            %Rate of dual-polarized antennas with SIC
            Rate_SIC_MMSE(m,xpd,n,:) = functionDownlinkMMSESICbound(H,W_MMSE,K,nbrOfRealizations,prelogFactor);
            Rate_SIC_ZF(m,xpd,n,:) = functionDownlinkMMSESICbound(H,W_ZF,K,nbrOfRealizations,prelogFactor);
            Rate_SIC_MR(m,xpd,n,:) = functionDownlinkMMSESICbound(H,W_MR,K,nbrOfRealizations,prelogFactor);
            
            
            %Rate expressions of uni-polarized arrays(UaTF downlink bound)
            Rate_UatF_MMSEunipl(m,xpd,n,:)= functionDownlinkMMSESICbound_unipl(Hunipl,W_MMSE_unipl,K,nbrOfRealizations,prelogFactor_unipl);
            Rate_UatF_ZFunipl(m,xpd,n,:)= functionDownlinkMMSESICbound_unipl(Hunipl,W_ZF_unipl,K,nbrOfRealizations,prelogFactor_unipl);
            Rate_UatF_MRunipl(m,xpd,n,:)= functionDownlinkMMSESICbound_unipl(Hunipl,W_MR_unipl,K,nbrOfRealizations,prelogFactor_unipl);
            
            
            
            %Theoretical rate calculation
            %Dual-polarized closed-form MR precoding
            Rate_th_MR(m,xpd,n,:) = functionDownlinkClosedForm(MMSEmatrixV,MMSEmatrixH,Rkv,Rkh,K,rho1,rho2,prelogFactor);
            %Uni-polarized closed-form MR precoding
            Rate_th_MR_unipl(m,xpd,n,:) = functionDownlinkClosedForm_unipl(MMSE_Cov_unipl,R_unipl,K,rho_unipl,prelogFactor_unipl);
            

            
            %Output simulation progress
            %disp([num2str(xpd) ' XPDs out of ' num2str(length(q_XPD))]);
        end
        %Output simulation progress
        disp([num2str(M) ' antennas of ' num2str(Mmax)]);
        
    end
    
    %Output simulation progress
    disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
end

Rate_th_MRsum=mean(sum(abs(Rate_th_MR),4),3);
Rate_SIC_MMSEsum=mean(sum(abs(Rate_SIC_MMSE),4),3);
Rate_SIC_ZFsum=mean(sum(abs(Rate_SIC_ZF),4),3);
Rate_SIC_MRsum=mean(sum(abs(Rate_SIC_MR),4),3);

%uni-polarized
Rate_th_MRunipl_sum=mean(sum(abs(Rate_th_MR_unipl),4),3);
Rate_UaTF_MMSEunipl_sum=mean(sum(abs(Rate_UatF_MMSEunipl),4),3);
Rate_UaTF_ZFunipl_sum=mean(sum(abs(Rate_UatF_ZFunipl),4),3);
Rate_UaTF_MRunipl_sum=mean(sum(abs(Rate_UatF_MRunipl),4),3);

ratioMMSE=mean(Rate_SIC_MMSEsum./Rate_UaTF_MMSEunipl_sum)
ratioZF=mean(Rate_SIC_ZFsum./Rate_UaTF_ZFunipl_sum)
ratioMR=mean(Rate_SIC_MRsum./Rate_UaTF_MRunipl_sum)



figure;
dp1=plot(Mrange, Rate_SIC_MMSEsum,'r-o');
hold on
dp2=plot(Mrange, Rate_SIC_ZFsum,'r->');
hold on
dp3=plot(Mrange, Rate_SIC_MRsum,'r-+');
hold on
plot(Mrange,Rate_th_MRsum,'bs');
hold on
up1=plot(Mrange,Rate_UaTF_MMSEunipl_sum,'k--o');
up2=plot(Mrange,Rate_UaTF_ZFunipl_sum,'k-->');
up3=plot(Mrange,Rate_UaTF_MRunipl_sum,'k--+');
m=plot(Mrange,Rate_th_MRunipl_sum, 'bs');
legend([dp1,dp2,dp3,up1,up2,up3,m],{'Dual-Polarized MMSE','Dual-Polarized ZF','Dual-Polarized MR',...
    'Uni-Polarized MMSE','Uni-Polarized ZF','Uni-Polarized MR','Closed-form MR'},'Interpreter','latex','Location','NorthWest');
xlabel('Number of Antennas (M)','Interpreter','latex');
ylabel('Average sum DL SE [bps/Hz]','Interpreter','latex');
set(gca,'fontsize',18);




%save the figure
paperPos = [0.361111111111111,2.25694444444445,7.77777777777778,5];

fig = gcf;
fig.PaperPosition = paperPos;
print(fig,'figDPvsUP','-depsc'); %This is an alternative way to write the same thing
saveas(fig,'figDPvsUP.fig');

