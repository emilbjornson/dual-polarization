% This Matlab script generates Figure 7 in the paper:
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

%XPC comparison
%Average downlink sum SE for 10 UEs with different preocoders
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
%The antenna
%The deterministic polarization matrix of the BS
F_BS= eye(2);
%The polarization matrix of the UE
F_UE=reshape(repmat(eye(2),1,K),2,2,K);



%Polarization correlation matrices
r_p = 0.8;
t_p = 0.8;
C_UE = [1, r_p; r_p', 1];
C_BS = [1, t_p; t_p', 1];

%XPD: polar-discrimination (same for all users)
%0<= q_XPD <= 1
%Small  q_XPD means LoS-like channels
XPDdB = 5;
q_XPD=  1/(1 + db2pow(XPDdB));


%Prepare for the simulations
nbrOfRealizations=100;
nbrOfSetups=1;
%number of maximum antennas
Mmax=100;
Mrange=40:10:Mmax;

%Pilot transmit powers
p1=100; %mW
p2=100; %mW
pV_UL= 100*ones(K,1);
pH_UL= 100*ones(K,1);


%Dowlink transmit power
rho1=100*ones(K,1);
rho2=100*ones(K,1);

%pilot length
tau_p=2*K;
prelogFactor= (tau_c-tau_p)/tau_c;


%Prepare to store the rates
%UE Side with SIC
Rate_SIC_MMSE=zeros(length(Mrange),length(q_XPD),nbrOfSetups,K);
Rate_SIC_ZF= zeros(length(Mrange),length(q_XPD),nbrOfSetups,K);
Rate_SIC_MR= zeros(length(Mrange),length(q_XPD),nbrOfSetups,K);
Rate_th_MR=zeros(length(Mrange),length(q_XPD),nbrOfSetups,K);

%UE Side without SIC
Rate_woSIC_MMSE=zeros(length(Mrange),length(q_XPD),nbrOfSetups,K);
Rate_woSIC_ZF=zeros(length(Mrange),length(q_XPD),nbrOfSetups,K);
Rate_woSIC_MR=zeros(length(Mrange),length(q_XPD),nbrOfSetups,K);


%Go through all setups
for n=1:nbrOfSetups
    
    
    %Compute spatial correlation matrices using the local scattering model
    
    %R_SS_normalized is the spatially-separeted DP spatial corr. matrix for Mmax x Mmax x K
    %R_SS_UE is the spatially-separeted spatial corr. matrix at UE side for 2 x 2 x K
    [R_SS_normalized,R_SS_UE,channelGaindB] = functionExampleSetup(L,K,Mmax,ASDdeg);
    
    %Compute the normalized average channel gain, where the normalization
    %is based on the noise power
    channelGainOverNoise = channelGaindB - noiseVariancedBm;
    %channelGainOverNoise= functionPowerControl( channelGaindB,noiseVariancedBm,deltadB,L);
    betas=10.^(channelGainOverNoise./10);
    
    
    
    
    for m=1:length(Mrange)
        
        M=Mrange(m);
        Mdual=M/2;
        R_uniSame=zeros(M,M,K);
        R_BSk=zeros(Mdual,Mdual,K);
        for k=1:K
            R_uniSame(:,:,k)= R_SS_normalized(1:M,1:M,k)*betas(k);
            R_BSk(:,:,k)=R_SS_normalized(1:Mdual,1:Mdual,k)*betas(k);
        end
        
        for xpd=1:length(q_XPD)
            
            
            Sigma = [sqrt(1-q_XPD(xpd)), sqrt(q_XPD(xpd)); sqrt(q_XPD(xpd)), sqrt(1-q_XPD(xpd))];
            
            [H,Rk_sqrtm] = functionChannelGeneration(Sigma,R_BSk, C_BS,C_UE,F_BS,F_UE,M,K,nbrOfRealizations);
            
            
            %Channel estimation part
            [vecHk_est,MMSEmatrixV, MMSEmatrixH,Rkv, Rkh] = functionChannelEstimation(H,Rk_sqrtm,q_XPD,M,K,nbrOfRealizations,tau_p,p1,p2);
            
            %Channel Estimation Error Matrices
            C_kV=Rkv-MMSEmatrixV;
            C_kH=Rkh-MMSEmatrixH;
            %Precoding is based on the estimated channels
            %vecHk_est 2M x nbrOfRealizations x K
            %Here are the estimated channels
            Hk_est=reshape(vecHk_est,M,2,nbrOfRealizations,K);
            %disp('Channel Estimation is DONE.')
            
            
            [W_MMSE,W_ZF,W_MR] = functionPrecoderMatrices(Hk_est,MMSEmatrixV,MMSEmatrixH,C_kV,C_kH,M,K,nbrOfRealizations,pV_UL,pH_UL,rho1,rho2);
            
            
            
            %DL Rate with SIC
            Rate_SIC_MMSE(m,xpd,n,:) = functionDownlinkMMSESICbound(H,W_MMSE,K,nbrOfRealizations,prelogFactor);
            Rate_SIC_ZF(m,xpd,n,:) = functionDownlinkMMSESICbound(H,W_ZF,K,nbrOfRealizations,prelogFactor);
            Rate_SIC_MR(m,xpd,n,:) = functionDownlinkMMSESICbound(H,W_MR,K,nbrOfRealizations,prelogFactor);
            
            %Downlink rate closed-form MR precoding
            Rate_th_MR(m,xpd,n,:) = functionDownlinkClosedForm(MMSEmatrixV,MMSEmatrixH,Rkv,Rkh,K,rho1,rho2,prelogFactor);
            
            
            
            
            %Now, Downlink without SIC part
            Rate_woSIC_MMSE(m,xpd,n,:) = functionDownlinkwithoutSICbound(H,W_MMSE,K,nbrOfRealizations,prelogFactor);
            Rate_woSIC_ZF(m,xpd,n,:) = functionDownlinkwithoutSICbound(H,W_ZF,K,nbrOfRealizations,prelogFactor);
            Rate_woSIC_MR(m,xpd,n,:) = functionDownlinkwithoutSICbound(H,W_MR,K,nbrOfRealizations,prelogFactor);
            
            
            
            %Output simulation progress
            %disp([num2str(xpd) ' XPDs out of ' num2str(length(q_XPD))]);
        end
        %         %Output simulation progress
        disp([num2str(M) ' antennas of ' num2str(Mmax)]);
        
    end
    
    %Output simulation progress
    disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
end

%UE Side SIC
Rate_SIC_MMSEsum=mean(sum(abs(Rate_SIC_MMSE),4),3);
Rate_SIC_ZFsum=mean(sum(abs(Rate_SIC_ZF),4),3);
Rate_SIC_MRsum=mean(sum(abs(Rate_SIC_MR),4),3);
Rate_th_MRsum=mean(sum(abs(Rate_th_MR),4),3);

%UE Side without SIC
Rate_woSIC_MMSEsum=mean(sum(abs(Rate_woSIC_MMSE),4),3);
Rate_woSIC_ZFsum=mean(sum(abs(Rate_woSIC_ZF),4),3);
Rate_woSIC_MRsum=mean(sum(abs(Rate_woSIC_MR),4),3);


figure;
m1=plot(Mrange, Rate_SIC_MMSEsum(:,1),'r-');
hold on
plot(Mrange, Rate_woSIC_MMSEsum(:,1),'k--')
hold on
m2=plot(Mrange, Rate_SIC_ZFsum(:,1),'r-o');
hold on
plot(Mrange, Rate_woSIC_ZFsum(:,1),'k--' );
hold on
m3=plot(Mrange, Rate_SIC_MRsum(:,1),'r->');
hold on
c=plot(Mrange, Rate_woSIC_MRsum(:,1),'k--' );
legend([m1,m2,m3 c],{'MMSE Precoding ',' ZF Precoding','MR Precoding', 'without SIC'},'Interpreter','latex','Location','NorthEast');
xlabel('Number of Antennas (M)','Interpreter','latex')
ylabel('Average DL Sum SE [bps/Hz]','Interpreter','latex');
% title(['q_{XPD}=',num2str(q_XPD(1))])
set(gca,'fontsize',18);


%save the figure
paperPos = [0.361111111111111,2.25694444444445,7.77777777777778,5];

fig = gcf;
fig.PaperPosition = paperPos;
saveas(fig,'figDownlinkXPC','fig');
saveas(fig,'figDownlinkXPC','epsc');
