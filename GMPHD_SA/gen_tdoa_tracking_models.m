function [model] = gen_tdoa_tracking_models(parameters, BayesoptResults,birthvelocity)
%gen_tdoa_tracking_models.m is a function that generates models required 
% for TDOA tracking with GMPHD-SA filter- this includes system and 
% measurement models, clutter and birth models.
%
% INPUTS:
% - parameters - a struct containing info on encounter, array, and
% parameters used for measurement extraction
% - BayesoptResults - a table containing results of Bayesian optimization 
% for certain parameters (consider this as a trained prior information).
% - birthvelocity - prior on the velocity component of newborn targets
% (learned from data)- 2 x N matrix, where first row are N bearings and
% second row are learned velocities.
%
% OUTPUTS:
% - model - a struct with 28 fields containing information required for 
% TDOA tracking with GMPHD-SA filter. 
%
%
% Pina Gruden, UH Manoa, April 2022. 


% BayesOptParams=BayesoptResults.XAtMinObjective;
BayesOptParams=BayesoptResults;

% ////////////////// ARRAY PARAMETERS /////////////////////////////////////
model.d=parameters.d; %sensor separation (in m)
model.c=parameters.c; %speed of sound (m/s)
model.m1=[0;0]; %sensor 1 position
model.m2=[model.d;0]; %sensor 2 position
model.taumax= model.d/model.c;%max possible TDOA
model.taumin= -model.d/model.c;%min possible TDOA
%//////////////////////////////////////////////////////////////////////////


%//////////////////////// SYSTEM (MOTION) MODEL ////////////////////////
%~~~~Choose CV model for tracking in the TDOA domain: ~~~~~~
%sampling period (step between windows):
model.dt= parameters.dt; 

%system matrix - determines transition from Xk_1 to Xk:
model.F = [1, model.dt; 0, 1]; 

%System noise variance matrix:
model.sig_sq_v=10^(BayesOptParams.sig_sq_v_exp);
model.Q = [model.dt^4/4,model.dt^3/2;model.dt^3/2, model.dt^2]*model.sig_sq_v;
%//////////////////////////////////////////////////////////////////////////


%//////////////////////// MEASUREMENT MODEL /////////////////////////
%WITH AMPLITUDE FEATURE:
%---------------------- TARGETS---------------------
% Expected target SNR bounds:
model.d1=10^(5/10); %SNR linear scale 10^(dB/10); lower bound on the expected SNR
model.d2=10^(20/10); %upper bound on the expected SNR
%Bounds estimated form training data

%Measurement extraction threshold (used on normalized cross-correlograms):
model.lambda= parameters.lambda; 

%Probability of target detection p_D:
% We learned it from training data since in bioacoustics p_D is a function
% of not only the SNR but also availability (probability that animals
% vocalizes- dependent on the vocalization rate)
model.p_D=BayesOptParams.p_D;
%~~~~~~~~Traditionally, if p_D only dependent on the SNR use below:~~~~~~~~
% fun_pd =@(d,A) A*1./(1+d); %Eq(18) Clark et al 2008 - this is distribution 
% %on expected SNR values; however we need to normalize it for the region
% %[d1,d2], need to get a scalar A that will then make the integral of fun_pd
% %equal to 1:
% pd_scalar= 1/integral(@(d) fun_pd(d,1),model.d1,model.d2);
% %to test wether we did ok: integral(@(d) fun_pd(d,pd_scalar),model.d1,model.d2)
% 
% %probability of detection given threshold lambda is then:
% p_D_lambda=@(d) fun_pd(d,pd_scalar).*exp(-model.lambda^2./(2.*(1+d))); %Eq(17) in Clark et al 2008
% model.p_D = integral(@(d) p_D_lambda(d),model.d1,model.d2); 
% % However this would be the case if pD only depended on the SNR - in
% % bioacoustics pD depends alos on the vocalization rate of the animals.
% %model.p_D will be trained with Bayesian optimization
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


%amplitude measeurement likelihood for targets (Clark et al.2008, Eq 19):
model.ga = @(a,d1,d2) 2.*(exp(-a.^2./(2*(1+d2)))-exp(-a.^2./(2*(1+d1))))./(a.*(log(1+d2)-log(1+d1)));

% Parameters for TDOA part of the meastuement:
model.H = [1,0];
%measurement noise variance of TDOA- estimated from training data:
model.R=10^(BayesOptParams.R_exp);  

%-------------- CLUTTER------------------------------
% Clutter amplitude likelihood (Eq(7) Clark et al 2008):
model.ca_lambda= @(a,lambda) a.*exp((lambda^2-a.^2)./2);

%mean clutter rate - learned from training data: 
model.c_lambda= BayesOptParams.c_lambda;  

%pdf of a clutter- uniform distiribution:
model.c_pdf=1/(model.taumax-model.taumin); 
%//////////////////////////////////////////////////////////////////////////


%//////////////////////////// BIRTH MODEL /////////////////////////////////

%birth intensity- learned from data:
model.birth=BayesOptParams.nbirth;

%birth velocity prior- learned from data based on bearing behaviour:
model.bpv_type=2; 
model.dottdoa_birth=birthvelocity;

% switch model.bpv_type
%     case 1
%         % Learned from one training data used in our TDOA tracking paper:
%         load('/Users/Pina/Dropbox/Pina/HAWAII/MATLAB/Code/Ground_Truth/Lasker_AC67_GT_1s_05overlap_GCCwhistles_learned_params.mat',...
%             'BestModel') %this is a GMM with 4 components learned on whistle training data
%         model.dottdoa_birth=BestModel; 
%     case 2
%         % Learned on all training data (lambda 4.5, chunked), and gives a
%         %typical velocity for each bearing
%         load('/Users/Pina/Dropbox/Pina/HAWAII/MATLAB/Code/Determine_GMPHD_params_tdoa/BirthVelocityPrior_AllTrainData.mat',...
%             'birthvelocity')
%         model.dottdoa_birth=birthvelocity;
% end
%//////////////////////////////////////////////////////////////////////////

 
%////////////////////////// GM-PHD PARAMETERS /////////////////////////////

% Probability of survival:
model.psurv = 0.99;

% Max number of Gaussian componenets used in Pruning:
model.Jmax = 100; 

% Merging threshold used in Pruning & merging:
model.U = BayesOptParams.U; 

% Pruning threshold used in Pruning & merging:
model.Tr = 0.001; 

% Weight threshold used in State estimation:
model.wth = BayesOptParams.wth; 
%//////////////////////////////////////////////////////////////////////////


end