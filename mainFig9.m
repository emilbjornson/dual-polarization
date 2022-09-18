% This Matlab script generates Figure 9 in the paper:
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
q_XPD=  1/(1 + db2pow(XPDdB));


%Prepare for the simulations
nbrOfRealizations=1;
nbrOfSetups=1e3;
%number of maximum antennas
Mmax=100;
Mrange=Mmax;

%Pilot transmit powers
p1=100; %mW
p2=100; %mW

%Dowlink transmit power
rhoTotal=K*100;
rhoDL=100;
% rho1=100*ones(K,1);
% rho2=100*ones(K,1);

%pilot length
tau_p=2*K;
prelogFactor= (tau_c-tau_p)/tau_c;


%Prepare to store the rates
%UE Side with SIC
% Rate_SIC_MMSE=zeros(length(Mrange),length(q_XPD),nbrOfSetups,K);
% Rate_SIC_ZF= zeros(length(Mrange),length(q_XPD),nbrOfSetups,K);
Rate_SIC_MR= zeros(length(Mrange),length(q_XPD),nbrOfSetups,K);
Rate_th_MR=zeros(length(Mrange),length(q_XPD),nbrOfSetups,K);
Rate_SIC_MRFull= zeros(length(Mrange),length(q_XPD),nbrOfSetups,K);
Rate_th_MRFull=zeros(length(Mrange),length(q_XPD),nbrOfSetups,K);


%UE Side without SIC
% Rate_woSIC_MMSE=zeros(length(Mrange),length(q_XPD),nbrOfSetups,K);
% Rate_woSIC_ZF=zeros(length(Mrange),length(q_XPD),nbrOfSetups,K);
Rate_woSIC_MR=zeros(length(Mrange),length(q_XPD),nbrOfSetups,K);
Rate_woSIC_MRFull=zeros(length(Mrange),length(q_XPD),nbrOfSetups,K);


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
            
            
            %Power Control
            ak=zeros(K,1);
            bk=zeros(K,1);
            lambdaInv=zeros(K,1);
            %First calculate the gammas
            for k=1:K
                ak(k)= abs(trace(MMSEmatrixV(:,:,k)));
                bk(k)= rhoTotal*betas(k) ;
                
                lambdaInv(k)= (1 + bk(k))/(ak(k));
                
            end
            
            %rhoTotal =\rho^\mathrm{dl}_\mathrm{tot}/2 in the paper
            rho1 = functionWaterfilling(rhoTotal,lambdaInv);
            
            rho2=rho1;
            %sum(rho1) =? rhoTotal : YES
            
            rho1Full=rhoDL*ones(K,1);
            rho2Full=rhoDL*ones(K,1);
            
            
            %Downlink rate closed-form MR precoding
            Rate_th_MR(m,xpd,n,:) = functionDownlinkClosedForm(MMSEmatrixV,MMSEmatrixH,Rkv,Rkh,K,rho1,rho2,prelogFactor);
            
            %Downlink rate closed-form MR precoding
            Rate_th_MRFull(m,xpd,n,:) = functionDownlinkClosedForm(MMSEmatrixV,MMSEmatrixH,Rkv,Rkh,K,rho1Full,rho2Full,prelogFactor);
            
            
            
            
            %Output simulation progress
            %disp([num2str(xpd) ' XPDs out of ' num2str(length(q_XPD))]);
        end
        %Output simulation progress
        %  disp([num2str(M) ' antennas of ' num2str(Mmax)]);
        
    end
    
    %Output simulation progress
    disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
end



%Plot the figure
CDFnumbers = linspace(0,1,nbrOfSetups);
MROpt=sum(reshape(Rate_th_MR(end,1,:,:),nbrOfSetups,K),2);
MRFull=sum(reshape(Rate_th_MRFull(end,1,:,:),nbrOfSetups,K),2);

figure;
hold on; box on;
plot(sort(abs(MRFull)),CDFnumbers,'k--')
plot(sort(abs(MROpt)),CDFnumbers,'r')
xlim([0 40])
legend({'Full Power','Max Sum Rate'},'Interpreter','latex','Location','NorthEast');
xlabel('Downlink Sum Rate [bps/Hz]')
ylabel('Cumulative Distribution Function');
set(gca,'fontsize',18);



%save the figure
paperPos = [0.361111111111111,2.25694444444445,7.77777777777778,5];
fig = gcf;
fig.PaperPosition = paperPos;
saveas(fig,'figDownlinkwithPC','fig');
saveas(fig,'figDownlinkwithPC','epsc');
