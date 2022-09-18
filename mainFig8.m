% This Matlab script generates Figure 8 in the paper:
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

%Uplink SE with power control
%Average Uplink Sum SE for K=10 UEs with different combining vectors.
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
nbrOfRealizations=100;
nbrOfSetups=500;
%number of maximum antennas
Mmax=100;
Mrange=Mmax;

%Pilot transmit powers
p1=100; %mW
p2=100; %mW
pUL=100; %mW
%pilot length
tau_p=2*K;
prelogFactor= (tau_c-tau_p)/tau_c;


%Prepare to store the rates
Rate_UaTFMR=zeros(length(Mrange),nbrOfSetups,K);
Rate_UaTFMRFull=zeros(length(Mrange),nbrOfSetups,K);
Rate_th_MR=zeros(length(Mrange),nbrOfSetups,K);
Rate_th_MRFull=zeros(length(Mrange),nbrOfSetups,K);

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
    
    pV_ulFull=100*ones(K,1);
    pH_ulFull=100*ones(K,1);
    %Go through different number of antennas at the BS side
    for m=1:length(Mrange)
        
        M=Mrange(m);
        Mdual=M/2;
        %Spatial Covariance Matrices
        R_BSk=zeros(Mdual,Mdual,K); %Co-located DP
        %Scale the normalized spatial correlation matrix
        for k=1:K
            R_BSk(:,:,k)=R_SS_normalized(1:Mdual,1:Mdual,k)*betas(k);
        end
        
        %Go through different channel XPDs
        for xpd=1:length(q_XPD)
            
            %2x2 Channel XPD matrix
            Sigma = [sqrt(1-q_XPD(xpd)), sqrt(q_XPD(xpd)); sqrt(q_XPD(xpd)), sqrt(1-q_XPD(xpd))];
            
            [H,Rk_sqrtm] = functionChannelGeneration(Sigma,R_BSk, C_BS,C_UE,F_BS,F_UE,M,K,nbrOfRealizations);
            
            
            
            %Channel estimation part
            [vecHk_est,MMSEmatrixV, MMSEmatrixH,Rkv, Rkh] = functionChannelEstimation(H,Rk_sqrtm,q_XPD,M,K,nbrOfRealizations,tau_p,p1,p2);
            a=zeros(K,1);
            bk=zeros(K,1);
            %First calculate the gammas
            for k=1:K
                bk(k)=abs(trace(MMSEmatrixV(:,:,k)*Rkv(:,:,k))/abs(trace(MMSEmatrixV(:,:,k))) + trace(MMSEmatrixV(:,:,k)*Rkh(:,:,k))/abs(trace(MMSEmatrixV(:,:,k))));
                
                a(k)= abs(trace(MMSEmatrixV(:,:,k)))/bk(k);
            end
            
            
            
            cvx_begin quiet
            variable x(K,1)
            variable s
            maximize log_det(eye(K)+diag(x).*diag(a))
            subject to
            
            x>=zeros(K,1);
            
            x<=s*pUL*betas;
            
            sum(x)==1-s;
            
            cvx_end
            
            
            %eta_UL_sumrate_k = x./(s*pUL.*bk);
            %pH_ul =pUL*eta_UL_sumrate_k;
            pH_ul =x./(s.*bk);
            pV_ul=pH_ul;
            
            
            
            %Channel Estimation Error Matrices
            C_kV=Rkv-MMSEmatrixV;
            C_kH=Rkh-MMSEmatrixH;
            %Precoding is based on the estimated channels
            %vecHk_est 2M x nbrOfRealizations x K
            %Here are the estimated channels
            Hk_est=reshape(vecHk_est,M,2,nbrOfRealizations,K);
            %disp('Channel Estimation is DONE.')
            
            
            
            
            
            
            
            %Uplink Spectral Efficiency Calculations
            
            
            %Second, UaTF Bound Part
            [v_kV,v_kH, vZF_kV,vZF_kH,vMR_kV,vMR_kH] = functionCombiningVectors(Hk_est,C_kV,C_kH,M,K,nbrOfRealizations,pV_ul,pH_ul);
            [v_kVFull,v_kHFull, vZF_kVFull,vZF_kHFull,vMR_kVFull,vMR_kHFull] = functionCombiningVectors(Hk_est,C_kV,C_kH,M,K,nbrOfRealizations,pV_ulFull,pH_ulFull);
            
            
            Rate_UaTFMR(m,n,:) = functionUplinkUaTFbound(H,vMR_kV,vMR_kH,K,nbrOfRealizations,pV_ul,pH_ul,prelogFactor);
            
            Rate_UaTFMRFull(m,n,:) = functionUplinkUaTFbound(H,vMR_kVFull,vMR_kHFull,K,nbrOfRealizations,pV_ulFull,pH_ulFull,prelogFactor);
            
            
            Rate_th_MR(m,n,:)= functionUplinkClosedForm(MMSEmatrixV,MMSEmatrixH,Rkv,Rkh,K,pV_ul,pH_ul,prelogFactor);
            Rate_th_MRFull(m,n,:)= functionUplinkClosedForm(MMSEmatrixV,MMSEmatrixH,Rkv,Rkh,K,pV_ulFull,pH_ulFull,prelogFactor);
            
     
            
            %Output simulation progress
            % disp([num2str(xpd) ' antennas of ' num2str(q_XPD)]);
        end
        
        %Output simulation progress
        disp([num2str(M) ' antennas of ' num2str(Mmax)]);
        
    end
    
    
    
    %Output simulation progress
    disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
end


% Plot the figure
CDFnumbers = linspace(0,1,nbrOfSetups);
MRFull=sum(reshape(Rate_UaTFMRFull(end,:,:),nbrOfSetups,K),2);
MROpt=sum(reshape(Rate_UaTFMR(end,:,:),nbrOfSetups,K),2);


figure;
hold on; box on;
plot(sort(MRFull),CDFnumbers,'k--')
plot(sort(MROpt),CDFnumbers,'r')
legend({'Full Power','Max Sum SE'},'Interpreter','latex','Location','NorthEast');
xlabel('Uplink Sum SE [bps/Hz]')
ylabel('Cumulative Distribution Function');
set(gca,'fontsize',18);

%save the figure
paperPos = [0.361111111111111,2.25694444444445,7.77777777777778,5];

fig = gcf;
fig.PaperPosition = paperPos;

saveas(fig,'figUplinkPC','fig');
saveas(fig,'figUplinkPC','epsc');

