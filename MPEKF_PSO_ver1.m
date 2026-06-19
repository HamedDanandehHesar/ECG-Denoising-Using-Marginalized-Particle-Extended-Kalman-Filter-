%%====teta and z are nonlinear
%%====w is linear
%%==== we use the proposed particle weighting strategy in 2015 paper close all


clear
 
%% Load ECG data from MAT file
[file, path] = uigetfile('*.mat','Select ECG mat file');
data = load(file);


fs = data.fs;                % Sampling frequency of ECG signal (Hz)

x = data.x;               % Extract stored signal matrix
ecg = x(1,:);             % Use first channel as ECG signal
length_sig = length(ecg); % Total number of ECG samples

ecg_bins = round(fs/2);           % Number of phase bins for ECG mean calculation

dt = 1/fs; % time step

%% adding white gaussian Noise
SNR = 6;   % you can change the Noise SNR level
x_noisy = zeros(3,length(x(1,:)));
x_noisy(1,:) = awgn(x(1,:),SNR,'measured');
%% -------- R-peak detection using Pan–Tompkins algorithm
[qrs_positions] = pantompkins_qrs(x_noisy(1,:),fs);
figure(1),plot(x_noisy(1,:),'b'),hold on,plot(qrs_positions,x_noisy(1,qrs_positions),'*r'),hold off
legend({'Noisy ECG Signal','R Peaks'})
title([file '   at SNR = ' num2str(SNR)])
axis('tight')
%% -------- Phase calculation
% Linear phase based on RR intervals
[Linearphase,~] = calculate_linear_phase_ver2(qrs_positions,length_sig,fs);
x_noisy(2,:) = Linearphase; 



[ECGsd,ECGmean,meanphase] = ECGsd_extractor_ver1(x_noisy(1,:),Linearphase,ecg_bins);

% further smoothing of ECG mean using wavelet
ECGmean = wdenoise(ECGmean,5,Wavelet="bior4.4",DenoisingMethod="BlockJS");



%% -------- ECG parameter extraction using Gaussian mixture model

MaxNumGaussian = 50;   % Maximum number of Strongest Gaussian components

% ========================building new myfun based on L Gaussians
L_num_of_Gaussian_kernels = 50;
ecg_mean_temp = 0;
ai = [];
bi = [];
tetai  = [];
for i=1:L_num_of_Gaussian_kernels
% disp(num2str(i))
ecg_mean_temp1 = ECGmean - ecg_mean_temp;
lb = [-1.5*max(ecg_mean_temp1).*ones(1,1)   0.000001*ones(1,1)   (-pi+.014)*ones(1,1)  ];
ub = [(1.5*max(ecg_mean_temp1)).*ones(1,1)  5*ones(1,1)  (pi-.014)*ones(1,1)  ];  
myfun1 = @(params)  norm(ecg_mean_temp1'-sum((repmat(params(1:1),ecg_bins,1).*exp(-(rem(repmat(meanphase,1,1)'-repmat(params(3),ecg_bins,1)+pi,2*pi)-pi) .^2 ./ (2*(repmat(params(2),ecg_bins,1)) .^ 2))),2));


% options = optimoptions('particleswarm','SwarmSize',30,'HybridFcn',@fmincon,'MaxIter',1000);
options = optimoptions('particleswarm','SwarmSize',50,'MaxIter',100,'Display','off');

OptimumParams = particleswarm(myfun1,3*1,lb,ub,options);

% L = (length(OptimumParams)/3);

ai_1 = OptimumParams(1);
bi_1 = OptimumParams(2);
tetai_1 = OptimumParams(3);
ai = [ai ai_1];
bi = [bi bi_1];
tetai  = [tetai tetai_1];
dtetai_1 = rem(meanphase - tetai_1 + pi,2*pi)-pi;
ecg_mean_temp = ecg_mean_temp + ai_1 .* exp(-dtetai_1 .^2 ./ (2*bi_1 .^ 2));
figure(41),plot(ecg_mean_temp,'b'),hold on,plot(ECGmean,'r')
legend({'Synthetic ECG','ECG Mean'}),hold off
title([num2str(i) 'th' '  Gaussian Kernel found'])
axis tight
% pause(3)
end

%% selection of the Strongest Peaks
[~,indx_strongest_peaks] = sort(abs(ai),'descend');

ai = ai(indx_strongest_peaks(1:MaxNumGaussian));
bi = bi(indx_strongest_peaks(1:MaxNumGaussian));
tetai = tetai(indx_strongest_peaks(1:MaxNumGaussian));


Alpha_i = ai;   % Gaussian amplitudes
Beta_i  = bi;   % Gaussian widths
Theta_i = tetai;   % Gaussian centers


% Sorting of parameters from based on tetai from -pi to pi

[Theta_i,idx] = sort(Theta_i,'ascend');
Alpha_i = Alpha_i(idx);
Beta_i = Beta_i(idx);
OptimumParams = [Alpha_i Beta_i Theta_i];
params = OptimumParams;
size_params = length(params);


%%  ANGULAR FREQUENCY MEASUREMENT

ind=1*qrs_positions;
ind2 = ind-[0 ind(1:end-1)];
RR = mean(ind2(1,2:end)); % mean of RR Intervals
w=2*pi*fs/RR; % angular frequency
RR_var=std(2*pi*fs./ind2(1,2:end));% standard deviation of RR Intervals







RR = mean(diff(ind(2:end-1)));

stepteta=2*pi/RR;
w_1 = fs*stepteta;
for j=ind(1,1):-1:1
    
     x_noisy(3,j) = w_1;

end   
for i=1:length(ind)-1



bins = ind(1,i+1)-ind(1,i);

 stepteta = 2*pi/(bins);
w_1 = fs*stepteta;
for j=ind(1,i)+1:ind(1,i+1)
 
    x_noisy(3,j) = w_1;

   
end

end
stepteta=2*pi/RR;
teta= 0;
w_1 = fs*stepteta;

for j=min(ind(end,end)+1,size(x_noisy,2)):size(x_noisy,2)

        x_noisy(3,j) = w_1;

 
end



%% MPEKF Initializations
 
 

RR_var = var(x_noisy(3,:));


 

Number_of_particles=400; %number of particles




%%========= synthetic ECG
ai = params(1,1:size_params/3);
bi = params(1,size_params/3+1:2*size_params/3);
tetai = params(1,2*size_params/3+1:size_params);
tetai(1:2) = tetai(1:2);
x_synth = zeros(1,size(x_noisy,2));
for i=1:length(ai)
    dtetai = rem(x_noisy(2,:)-tetai(i)+pi,2*pi)-pi;
    x_synth = x_synth+ai(i)*exp(-dtetai.^2./(2*bi(i).^2));
end


% x_noisy_copy(1,:) = x_noisy_copy_no_baseline(1,:);

% ind=3*m.ind_Rpeak;





 




y = [x_noisy(1,:);x_noisy(2,:);x_noisy(3,:)];










x_initial = [x_noisy(1,1);x_noisy(2,1);x_noisy(3,1)];
Number_of_steps = size(y,2);
%Prepare for filtering
%space for recording er(n), xe(n)
rer=zeros(3,Number_of_steps); 
MPEKF_output=zeros(3,Number_of_steps);
%------------------------------------------
x_initial_estimate = x_initial; %initial estimation
%prepare particles
%reserve space
particle_matrix=zeros(3,Number_of_particles); %particles
apriori_particle_matrix=zeros(3,Number_of_particles); %a priori particles
particle_measurement=zeros(3,Number_of_particles); %particle measurements
particle_measurement_difference=zeros(3,Number_of_particles); %meas. dif.
particle_likelihood_or_weights=zeros(1,Number_of_particles); %particle likelihoods
%particle generation
noise_for_initial_particles=randn(3,Number_of_particles); %noise (initial particles)


% state_noise_vec = [   (.1*mean(ECGsd(1,1:round(length(ECGsd)))))^2 (.01*w/fs).^2/12 (.1*RR_var)^2  (.01*ai.*ones(1,size_params/3)).^2 (.01*bi.*ones(1,size_params/3)).^2 (.000005*ones(1,size_params/3)).^2  ] ;
% measure_noise_vec = [(.2*mean(ECGsd(1,1:round(length(ECGsd))))^2) (1*w/fs).^2/12 (1*RR_var.^2)];

 

state_noise_vec = [(.1*mean(ECGsd(1,1:round(length(ECGsd)))))^2 (0.00001) (1*RR_var.^2)  (0.1*ai.*ones(1,size_params/3)).^2 (.1*bi.*ones(1,size_params/3)).^2 (.001*ones(1,size_params/3)).^2  ] ;

measure_noise_vec = [(1*mean(ECGsd(1,1:round(length(ECGsd))))^2)  (0.01) (.0001*RR_var.^2)];


nonlin_state_noise_vec = [(.1*mean(ECGsd(1,1:round(length(ECGsd)))))^2 (0.0001)  (0.1*ai.*ones(1,size_params/3)).^2 (.1*bi.*ones(1,size_params/3)).^2 (.001*ones(1,size_params/3)).^2 ];


std_Q_vec = sqrt(state_noise_vec);
std_R_vec = sqrt(measure_noise_vec);

Q = diag(state_noise_vec);
R = diag(measure_noise_vec);


for ip=1:Number_of_particles,
particle_matrix(:,ip)=x_initial+(1*std_Q_vec(1:3)'.*noise_for_initial_particles(:,ip)); %initial particles
end;

state_noise_matrix=randn(length(std_Q_vec),Number_of_steps); %process
measure_noise_matrix=randn(length(std_R_vec),Number_of_steps); %output

noise_state_for_particles=randn(3,Number_of_particles); %part of noise (process) that has function
noise_measurement_for_particles=randn(3,Number_of_particles); %noise (measurement)
particle_matrix_cell = {};
particle_matrix_weights = {};

%=======New line==============
Pp = ones(1,1,Number_of_particles); %covariance of linear part (w)
H =1;
A_linear = 1;
%====================================


for nn=1:size(y,2)
   %estimation recording
MPEKF_output(:,nn)=x_initial_estimate; %state
rer(:,nn)=x_initial-x_initial_estimate; %error
%Simulation of the system
%system
y_measure = y(:,nn);
ai_matrix = repmat(params(1,1:size_params/3),Number_of_particles,1)+repmat(std_Q_vec(4:4+size_params/3-1),Number_of_particles,1).*randn(Number_of_particles,length(ai));
bi_matrix = repmat(params(1,size_params/3+1:2*size_params/3),Number_of_particles,1)+repmat(std_Q_vec(4+size_params/3:4+2*size_params/3-1),Number_of_particles,1).*randn(Number_of_particles,length(bi));
tetai_matrix = repmat(params(1,2*size_params/3+1:size_params),Number_of_particles,1)+repmat(std_Q_vec(4+2*size_params/3:end),Number_of_particles,1).*randn(Number_of_particles,length(tetai));
dtetai_matrix = rem(repmat(particle_matrix(2,:)',1,length(bi))-tetai_matrix+pi,2*pi)-pi;
w_matrix = repmat(particle_matrix(3,:)',1,length(bi));
% 
temp_matrix = w_matrix.*(ai_matrix./(bi_matrix.^2)).*dtetai_matrix.*exp(-dtetai_matrix.^2./(2*bi_matrix.^2));

apriori_particle_matrix(1,:) = particle_matrix(1,:)-dt*(sum(temp_matrix,2))';
apriori_particle_matrix(2,:) =   particle_matrix(2,:)+particle_matrix(3,:)*dt ;
%%=======new line===============================================
% apriori_particle_matrix(3,:)= particle_matrix(3,:);  this part is linear 
apriori_particle_matrix(1:2,:) = apriori_particle_matrix(1:2,:)+(repmat(std_Q_vec(1:2)',1,Number_of_particles).*randn(2,Number_of_particles)); %additive noise
apriori_particle_matrix(2,:) = rem(apriori_particle_matrix(2,:)+pi,2*pi)-pi;

ksi = apriori_particle_matrix(1:2,:)-particle_matrix(1:2,:);

%==============================================================================

ip=1:Number_of_particles;
ai_matrix = repmat(params(1,1:size_params/3),Number_of_particles,1);
bi_matrix = repmat(params(1,size_params/3+1:2*size_params/3),Number_of_particles,1);
tetai_matrix = repmat(params(1,2*size_params/3+1:size_params),Number_of_particles,1);
dtetai_matrix = rem(repmat(particle_matrix(2,:)',1,length(bi))-tetai_matrix+pi,2*pi)-pi;
w_matrix = repmat(particle_matrix(3,:)',1,length(bi));
Al(2,1,ip) = dt;%% dF2/w

Bl(1,1,ip) = 1;%.05;1;sqrt(.5);% dF1/eta
Bl(2,2,ip) = 1;% dF2/eta1   eta1 is noise of teta
% Aw(:,:,ip) = - dt*ai_matrix.*dtetai_matrix./(bi_matrix.^2).*exp(-dtetai_matrix.^2./(2*bi_matrix.^2));
Al(1,1,ip) = sum(- dt*ai_matrix.*dtetai_matrix./(bi_matrix.^2).*exp(-dtetai_matrix.^2./(2*bi_matrix.^2)),2);%%dF1/w

Aa =-dt*w_matrix./(bi_matrix.^2).*dtetai_matrix .* exp(-dtetai_matrix.^2./(2*bi_matrix.^2));
Bl(1,3:size_params/3+2,ip) = Aa' ;%% dF1/ai
Ab = 2*dt.*ai_matrix.*w.*dtetai_matrix./bi_matrix.^3.*(1 - dtetai_matrix.^2./(2*bi_matrix.^2)).*exp(-dtetai_matrix.^2./(2*bi_matrix.^2));
Bl(1,size_params/3+3:2*size_params/3+2,ip)= Ab';%% dF1/bi
At = dt*w_matrix.*ai_matrix./(bi_matrix.^2).*exp(-dtetai_matrix.^2./(2*bi_matrix.^2)) .* (1 - dtetai_matrix.^2./bi_matrix.^2);
Bl(1,2*size_params/3+3:size_params+2,ip) = At';  %% dF1/tetai


% Al(2,1) = dt;%% dF2/w
% 
% Bl(1,1) = 1;%.05;1;sqrt(.5);% dF1/eta
% Bl(2,2) = 1;% dF2/eta1   eta1 is noise of teta
% 
% ai = params(1,1:size_params/3);%ai_matrix(ip,:);%
% bi = params(1,size_params/3+1:2*size_params/3) ;%bi_matrix(ip,:);%
% tetai = params(1,2*size_params/3+1:size_params);%tetai_matrix(ip,:);%

for ip=1:Number_of_particles

    



%%======================prediction of linear part========================

error_nonlin = ksi(:,ip);

N = Al(:,:,ip)*Pp(:,:,ip)*Al(:,:,ip)'+ Bl(:,:,ip)*diag(nonlin_state_noise_vec)*Bl(:,:,ip)'+eps ;
L = A_linear*Pp(:,:,ip)*Al(:,:,ip)'/N ;
Pp(:,:,ip) = (A_linear*Pp(:,:,ip)*A_linear' + Q(3,3) - L*N*L');
apriori_particle_matrix(3,ip) = particle_matrix(3,ip)+ ...
    L*(error_nonlin+dt*sum(particle_matrix(3,ip)*ai./(bi.^2).*dtetai_matrix(ip,:).*exp(-dtetai_matrix(ip,:).^2./(2*bi.^2))));

    
%  Kalman Filter Measurement Update

S = H*Pp(:,:,ip)*H' + 1*R(3,3);
K = Pp(:,:,ip)*H'/S;



apriori_particle_matrix(3,ip) = apriori_particle_matrix(3,ip) + K*(y_measure(3,1) - H*apriori_particle_matrix(3,ip) );


Pp(:,:,ip) = Pp(:,:,ip) - K*H*Pp(:,:,ip);




end
%Likelihood


particle_measurement(:,:)= apriori_particle_matrix(:,:);

ip=1:Number_of_particles;
% particle_likelihood_or_weights(ip)= mvnpdf(particle_measurement(:,ip)',y_measure',R); 
% particle_likelihood_or_weights(ip)= mvnpdf(particle_measurement(1,ip)',y_measure(1)',R(1,1));
particle_likelihood_or_weights(ip)= mvnpdf(particle_measurement(1,ip)',y_measure(1)',R(1,1))+mvnpdf(particle_measurement(1,ip)',x_synth(1,nn)',R(1,1));% gives best results
% particle_likelihood_or_weights(ip)= mvnpdf(particle_measurement(1,ip)',y_measure(1)',R(1,1))+mvnpdf(particle_measurement(1,ip)',x_synth(1,nn)',R(1,1));%gives better results

  

sum_of_weights=nansum(particle_likelihood_or_weights);
particle_likelihood_or_weights(isnan(particle_likelihood_or_weights))=0;




%normalization
if sum_of_weights==0
particle_likelihood_or_weights(ip) = 1/Number_of_particles;
else
particle_likelihood_or_weights(ip)=particle_likelihood_or_weights(ip)/sum_of_weights;

end

Neff = 1./(sum(particle_likelihood_or_weights.^2));

if Neff <Number_of_particles
%    disp('systematic_resampling');
%Prepare for roughening
A1=(max(apriori_particle_matrix(1,:)')-min(apriori_particle_matrix(1,:)'))';
A2=(max(apriori_particle_matrix(2,:)')-min(apriori_particle_matrix(2,:)'))';
A3=(max(apriori_particle_matrix(3,:)')-min(apriori_particle_matrix(3,:)'))';
A = [A1;A2;A3];
sig= 5*A*Number_of_particles^(-1/3);
random_numbers1= .0001*randn(1,Number_of_particles); %random numbers
random_numbers2= .0001*randn(1,Number_of_particles); %random numbers
random_numbers3= .0001*randn(1,Number_of_particles); %random numbers
random_numbers   = [random_numbers1;random_numbers2;random_numbers3];
%================================
% %Resampling (systematic)
[particle_likelihood_or_weights,id] = sort(particle_likelihood_or_weights);
apriori_particle_matrix = apriori_particle_matrix(:,id);
cumsum_weights=cumsum(particle_likelihood_or_weights);
comb=linspace(0,1-(1/Number_of_particles),Number_of_particles)+(rand(1)/Number_of_particles); %the "comb"
comb(Number_of_particles+1)=1;
ip=1; mm=1;
while(ip<=Number_of_particles) && (mm<=Number_of_particles)
if (comb(ip)<cumsum_weights(mm))
temp_particle=apriori_particle_matrix(:,mm);
temp_covariance = Pp(:,mm);
particle_likelihood_or_weights(ip) = particle_likelihood_or_weights(mm);
apriori_particle_matrix(:,ip)=temp_particle+(sig.*random_numbers(:,ip)); %roughening
Pp(:,:,ip) = temp_covariance;
if apriori_particle_matrix(2,ip)> pi
    apriori_particle_matrix(2,ip) = apriori_particle_matrix(2,ip)-2*pi;
end
ip=ip+1;
else
mm=mm+1;
end;
end;
end;
particle_matrix = apriori_particle_matrix;



%=================================
%Results
%estimated state (the particle mean)
% x_initial_estimate=sum([particle_likelihood_or_weights;particle_likelihood_or_weights;particle_likelihood_or_weights].*particle_matrix,2);
% x_initial_estimate=sum(particle_matrix,2)/Number_of_particles;
x_initial_estimate(1) = mean(particle_matrix(1,:));sum(particle_likelihood_or_weights.*particle_matrix(1,:))/sum(particle_likelihood_or_weights); 
x_initial_estimate(2) =  mean(particle_matrix(2,:));sum(particle_likelihood_or_weights.*particle_matrix(2,:))/sum(particle_likelihood_or_weights); 
x_initial_estimate(3) =  mean(particle_matrix(3,:));sum(particle_likelihood_or_weights.*particle_matrix(3,:))/sum(particle_likelihood_or_weights); 

MPEKF_output(:,nn)=x_initial_estimate; %state
  
disp(['Processing sample number:' num2str(nn) '/' num2str(size(MPEKF_output,2))])
particle_matrix_cell{1,nn}= particle_matrix;
particle_matrix_weights{1,nn} = particle_likelihood_or_weights;

end

toc


% figure(1),plot(1:size(y,2),y(1,:),'r',1:size(y,2),MPEKF_output(1,:),'b')


MPEKF_SNR = 10*log10(mean((x(1,:)-x_noisy(1,:)).^2)/mean((x(1,:)-MPEKF_output(1,:)).^2))

 



 

% figure(2),plot(1:size(y,2),y(2,:),'r',1:size(y,2),MPEKF_output(2,:),'b')
% 
% figure(3),plot(1:size(y,2),y(3,:),'r',1:size(y,2),MPEKF_output(3,:),'b')









figure(2),
subplot(3,1,1)
        plot(1:length(x),x(1,:),'k')
                legend({'Original'})
                        axis('tight')
        title([file])
        subplot(3,1,3)
        plot(1:length(x),MPEKF_output(1,:),'r')
                legend({'MPEKF'})
                        axis('tight')

        subplot(3,1,2),plot(1:length(x),x_noisy,'b')
        legend({'Noisy'})
        axis('tight')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions
function [Phase,Omega] = calculate_linear_phase_ver2(locs,length_sig,fs)

% locs       : indices of detected R-peaks
% length_sig : total number of ECG samples
% fs         : sampling frequency

ind = locs(:)';                % Convert R‑peak indices to row vector

Phase = zeros(1,length_sig);   % Phase of each ECG sample
Omega = zeros(1,length_sig);   % Instantaneous angular frequency

RR = mean(diff(ind));          % Mean RR interval (samples)

%% -------- Phase before the first R‑peak

stepTheta = 2*pi/RR;           % Average phase increment per sample
omega_val = fs*stepTheta;      % Instantaneous angular frequency

theta = 0;                     % Initialize phase

for j = ind(1)-1:-1:1          % Move backward from first R‑peak
    theta = theta - stepTheta; % Decrease phase
    theta = mod(theta+pi,2*pi)-pi; % Wrap phase into [-pi , pi]

    Phase(j) = theta;          % Store phase
    Omega(j) = omega_val;      % Store frequency
end

%% -------- Phase between consecutive R‑peaks

for k = 1:length(ind)-1

    bins = ind(k+1)-ind(k);    % Number of samples between R-peaks

    stepTheta = 2*pi/bins;     % Phase increment so phase spans one cycle
    omega_val = fs*stepTheta;  % Corresponding angular frequency

    theta = 0;
    Phase(ind(k)) = 0;         % Define phase at R‑peak as zero

    for j = ind(k)+1 : ind(k+1)-1
        theta = theta + stepTheta; % Linear phase progression
        if theta>pi
            theta = -pi;
        end
        Phase(j) = theta;
        Omega(j) = omega_val;
    end

    Phase(ind(k+1)) = 0;       % Next R‑peak also set to zero phase
end

%% -------- Phase after the last R‑peak

stepTheta = 2*pi/RR;           % Use mean RR again
omega_val = fs*stepTheta;

theta = 0;

for j = ind(end)+1:length_sig
    theta = theta + stepTheta; % Continue phase linearly
    theta = mod(theta+pi,2*pi)-pi;

    Phase(j) = theta;
    Omega(j) = omega_val;
end

end



function [ECGsd,ecg_mean,phase_mean] = ECGsd_extractor_ver1(ecg,phase,bins)

x1 = ecg;                        % ECG signal
meanPhase = zeros(1,bins);       % Mean phase per bin
ECGmean = zeros(1,bins);         % Mean ECG per bin
ECGsd = zeros(1,bins);           % ECG standard deviation per bin

% Handle wrap-around phase bin near -pi / +pi
I = find( phase >= (pi-pi/bins) | phase < (-pi+pi/bins) );

if(~isempty(I))
    meanPhase(1) = -pi;
    ECGmean(1) = mean(x1(I));
    ECGsd(1) = std(x1(I));
else
    ECGsd(1) = -1;               % Mark empty bins
end

% Loop over phase bins
for i = 1 : bins-1
    I = find( phase >= 2*pi*(i-0.5)/bins - pi & ...
              phase <  2*pi*(i+0.5)/bins - pi );

    if(~isempty(I))
        meanPhase(i+1) = mean(phase(I));
        ECGmean(i+1) = mean(x1(I));
        ECGsd(i+1) = std(x1(I));
    else
        ECGsd(i+1) = -1;
    end
end

% Interpolate missing bins
K = find(ECGsd==-1);

for i = 1:length(K)
    switch K(i)
        case 1
            meanPhase(1) = -pi;
            ECGmean(1) = ECGmean(2);
            ECGsd(1) = ECGsd(2);
        case bins
            meanPhase(bins) = pi;
            ECGmean(bins) = ECGmean(bins-1);
            ECGsd(bins) = ECGsd(bins-1);
        otherwise
            meanPhase(K(i)) = mean(meanPhase([K(i)-1 K(i)+1]));
            ECGmean(K(i))   = mean(ECGmean([K(i)-1 K(i)+1]));
            ECGsd(K(i))     = mean(ECGsd([K(i)-1 K(i)+1]));
    end
end

phase_mean = meanPhase;
ecg_mean   = ECGmean;
ECGsd      = ECGsd;

end





