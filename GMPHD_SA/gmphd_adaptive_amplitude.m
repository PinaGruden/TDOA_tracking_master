function [Est] = gmphd_adaptive_amplitude(model,measure)

% Pina Gruden, 2020-2021

% Based on Gruden, Nosal, Oleson (2021), JASA; Vo & Ma (2006), IEEE
% Trans.Sig.Proces.; and Ristic, Clark, Vo & Vo (2012), IEEE Trans. Aerosp. Electron. Syst.

% OUTPUT:
% Structure array Est with the following fields:
% - Est.X = estimated states (means of Gaussian components)
% - Est.w = weights associated to the estimates
% - Est.N = estimated number of targets
% - Est.Tag = identities of the estimates




% Output pre-allocation
% Est.gmphd = gmphd;
Est.X= cell(measure.T,1); %estimated states (means of Gaussian components)
Est.w= cell(measure.T,1);
Est.N = zeros(measure.T,1); %estimated number of targets
Est.Tag = cell(measure.T,1); %Identity tags of the estimates

% GM-PHD filter

% Initialization - Change this depending on the problem!
J_update=randi(3);  % number of Gaussian components (draw a random number between 1 and 3)
w_update=repmat(0.05/J_update,1,J_update); % weights of Gaussian components
m_update(1,:)=round(rand(1,J_update).* (model.taumax - model.taumin)+ model.taumin); % means of Gaussian components
m_update(2,:)=zeros(1,J_update);
P_update = repmat(0.01*model.Q, [1,1,J_update]);  % covariances of Gaussian components
Tag_update = zeros(1,J_update);
id=1;

for k=1:measure.T
%    For troubleshooting use Breakpoints option (in Editor tab)

%/////////////// Prediction ///////////////////
    %~~~~~~~~~~~~ Persistent targets ~~~~~~~~~~~~~~~ 
    w_predict = model.psurv*w_update;
    [m_predict,P_predict] = kalman_predict_multiple(model,m_update,P_update);
    Tag_predict=Tag_update;
      
    %~~~~~~~~~~~~ Newborn targets ~~~~~~~~~~~~~~~~~~
%     %GENERATE NEWBORN- draw neborn states from a GMM centered on
%     measurements
     [mb,wb,Pb,Tag_birth,id] = birth_process_amplitude(model,measure.Z{k},id);
     
%----------------------------------------------------------
%      %Choose to Update newborn and persistent jointly:
%      m_predict= cat(2,mb,m_predict);             
%      w_predict= cat(2,wb,w_predict);
%      P_predict=cat(3,Pb,P_predict);

    %/////////////// Update ///////////////////////
    % Separately for persistent and newborn (as per Ristic et al (2012))
    %~~~~~~~~~~~~~~~Persistent targets~~~~~~~~~~~~~~~~~~
    % Missed targets
    w_update=(1-model.p_D)*w_predict;
    m_update=m_predict;
    P_update=P_predict;
    Tag_update=Tag_predict;
   
    if ~isempty(measure.Z{k}) %if there are measurements to update with
        %compute means and covariance with Kalman update, and likelihood
        %function for position (tdoa) component of the measurement:
        [qz_temp,m_temp,P_temp,Tag_temp] = kalman_update_multiple(measure.Z{k}(1,:),model,m_predict,P_predict,Tag_predict);
        [qzb_temp,mb_temp,Pb_temp,Tagb_temp] = kalman_update_multiple(measure.Z{k}(1,:),model,mb,Pb,Tag_birth);
        
        for indx=1:size(measure.Z{k},2)
%            
            w_temp = model.ga(measure.Z{k}(2,indx),model.d1,model.d2)*w_predict.*qz_temp(:,indx)'; 
%             w_temp = model.p_D*w_predict(:).*qz_temp(:,indx); %Use this is w is a column vector
            
%             %Choose to Update newborn and persistent jointly:
%             w_temp= w_temp./(model.c_lambda*model.c_pdf + sum(w_temp)); 
%             w_update = cat(2,w_update,w_temp);
%             m_update = cat(2,m_update,m_temp(:,:,indx));
%             P_update = cat(3,P_update,P_temp);
            
            %Choose to Update newborn and persistent separately:
            wb_temp= model.ga(measure.Z{k}(2,indx),model.d1,model.d2)*wb.*qzb_temp(:,indx)'./...  %wb./(model.taumax-model.taumin) Previously I had wb.*qzb_temp(:,indx)'./ || repmat(model.birth/(model.taumax-model.taumin),1,length(Tagb_temp)) // Also I had wb.*qzb_temp(:,indx)'.*(model.birth/(sum(wb.*qzb_temp(:,indx)')+1e-20))./(model.taumax-model.taumin)
                (model.c_lambda*model.c_pdf*model.ca_lambda(measure.Z{k}(2,indx),model.lambda)+... clutter term with amplitude likelihood
                model.ga(measure.Z{k}(2,indx),model.d1,model.d2)*wb*qzb_temp(:,indx)+... birth term %Previously I had sum(wb.*qzb_temp(:,indx)')
                sum(w_temp)); %target term with amplitude likelihood
            w_temp= w_temp./...
                (model.c_lambda*model.c_pdf*model.ca_lambda(measure.Z{k}(2,indx),model.lambda)+... clutter term with amplitude likelihood
                model.ga(measure.Z{k}(2,indx),model.d1,model.d2)*wb*qzb_temp(:,indx)+... birth term %Previously I had sum(wb.*qzb_temp(:,indx)')
                sum(w_temp)); %target term with amplitude likelihood

            if any(isnan(wb_temp))
                wb_temp(isnan(wb_temp))=0;
            end
            
            if any(isnan(w_temp))
                w_temp(isnan(w_temp))=0;
            end

            w_update = cat(2,w_update,w_temp,wb_temp);
            m_update = cat(2,m_update,m_temp(:,:,indx),mb_temp(:,:,indx)); 
            %mb_temp(:,:,indx) gets Gaussians that were updated with measurement measure.Z{k}(:,indx)
            P_update = cat(3,P_update,P_temp,Pb_temp);
%             Tag_update = cat(2,Tag_update,Tag_temp{indx},Tagb_temp{indx});
            Tag_update = cat(2,Tag_update,Tag_temp(indx,:),Tagb_temp(indx,:));
        end
    end
    
    
        
    %////////////// Prune, Merge, and Limit Max N components ///////////////
   
    [w_update,m_update,P_update,Tag_update] = gmphd_prune(w_update,m_update,P_update,model.Tr,Tag_update);
    [w_update,m_update,P_update,Tag_update] = gmphd_merge(w_update,m_update,P_update,model.U,Tag_update);
    [w_update,m_update,P_update,Tag_update] = gmphd_limit(w_update,m_update,P_update,model.Jmax,Tag_update);
    
    %///////////// Estimate States ////////////////
    idx = find(w_update > model.wth);
    Est.X{k}=m_update(:,idx);
    Est.N(k)=size(Est.X{k},2);
    Est.w{k}=w_update(idx);
    if ~isempty(idx)
    Est.Tag{k}=Tag_update(idx);
    else
    Est.Tag{k}=[];   
    end
%     for j=1:length(idx)
%         repeat_num_targets= round(w_update(idx(j))); %this is in case there are two targets at the same location
%         Est.X{k}= [ Est.X{k} repmat(m_update(:,idx(j)),[1,repeat_num_targets]) ];
%         Est.N(k)= Est.N(k)+repeat_num_targets;
%     end
end


end

