function [Track] = tracktarget_labels_1D_smoothed(XTag, Xk, CovX, wX, model)
    IDlist = unique([XTag{:}]);
    N = numel(IDlist); % number of all tracks
    if N == 0 % there are no tracks
        Track = [];
        return
    end
    
    % Preallocate
    Track(N).est = [];
    Track(N).tdoa = [];
    Track(N).dottdoa = [];
    Track(N).time = [];
    Track(N).label = [];
    Track(N).ti = [];
    Track(N).P = []; % Covariance
    Track(N).w = []; % Weight
    
    t = 0:model.dt:size(XTag,1)*model.dt-model.dt;
    
    % Forward Pass
    for k = 1:size(XTag,1)
        if ~isempty(XTag{k})
            [ni, indxs] = histc(XTag{k}, IDlist);
            multiple = find(ni > 1);
            ID = ismember(indxs, multiple);
            
            % For IDs that are NOT repeated
            indx = indxs(ID==0);
            for m = 1:length(indx)
                M = indx(m);
                idx = find(indxs == M);
                Track(M).est = [Track(M).est, Xk{k}(:,idx)];
                Track(M).P{end+1} = CovX{k}(:,:,idx);
                Track(M).w = [Track(M).w, wX{k}(idx)];
                Track(M).time = [Track(M).time; t(k)];
                Track(M).label = [Track(M).label; IDlist(M)];
                Track(M).ti = [Track(M).ti, k];
            end
            
            % For IDs that ARE repeated
            if ~isempty(multiple)
                for m = 1:length(multiple)
                    M = multiple(m);
                    idx = find(indxs == M);
                    detections = Xk{k}(:,idx);
                    covs = CovX{k}(:,:,idx);
                    weights = wX{k}(idx);
                    
                    if isempty(Track(M).est)
                        % Initialize with the highest weighted estimate
                        [~, max_idx] = max(weights);
                        Track(M).est = [Track(M).est, detections(:,max_idx)];
                        Track(M).P{end+1} = covs(:,:,max_idx);
                        Track(M).w = [Track(M).w, weights(max_idx)];
                    else
                        % Predict, accounting for time gap
                        Nstp = k - Track(M).ti(end);
                        Fnew = model.F;
                        Fnew(1,2) = Fnew(1,2) * Nstp;
                      
                        
                        [x_pred, P_pred] = kf_predict(Track(M).est(:,end), Track(M).P{end}, Fnew, model.Q * Nstp);
                        
                        % Compute Mahalanobis distances
                        mah_dist = zeros(1, size(detections, 2));
                        for i = 1:size(detections, 2)
                            diff = detections(:,i) - x_pred;
                            S = P_pred + covs(:,:,i);
                            mah_dist(i) = sqrt(diff' / S * diff);
                        end
                        
                        % Combine Mahalanobis distance and weight
                        combined_score = mah_dist ./ weights;
                        [~, best_idx] = min(combined_score);
                        
                        % Update
                        z = detections(:,best_idx);
                        R = covs(:,:,best_idx);
                        [x_update, P_update] = kf_update(x_pred, P_pred, z, eye(2), R);
                        
                        Track(M).est = [Track(M).est, x_update];
                        Track(M).P{end+1} = P_update;
                        Track(M).w = [Track(M).w, weights(best_idx)];
                    end
                    Track(M).time = [Track(M).time; t(k)];
                    Track(M).label = [Track(M).label; IDlist(M)];
                    Track(M).ti = [Track(M).ti, k];
                end
            end
        end
    end
    
    % Backward Pass (RTS Smoothing)
    for i = 1:N
        n = size(Track(i).est, 2);
        if n > 1
            x_smooth = Track(i).est;
            P_smooth = Track(i).P;
            for k = n-1:-1:1
                [x_pred, P_pred] = kf_predict(x_smooth(:,k), P_smooth{k}, model.F, model.Q);
                G = P_smooth{k} * model.F' / P_pred;
                x_smooth(:,k) = x_smooth(:,k) + G * (x_smooth(:,k+1) - x_pred);
                P_smooth{k} = P_smooth{k} + G * (P_smooth{k+1} - P_pred) * G';
            end
            Track(i).est = x_smooth;
            Track(i).tdoa = x_smooth(1,:);
            Track(i).dottdoa = x_smooth(2,:);
            Track(i).P = P_smooth;
        else
            Track(i).tdoa = Track(i).est(1,:);
            Track(i).dottdoa = Track(i).est(2,:);
        end
    end
end

function [x_pred, P_pred] = kf_predict(x, P, F, Q)
    x_pred = F * x;
    P_pred = F * P * F' + Q;
end

function [x_update, P_update] = kf_update(x_pred, P_pred, z, H, R)
    y = z - H * x_pred;
    S = H * P_pred * H' + R;
    K = P_pred * H' / S;
    x_update = x_pred + K * y;
    P_update = (eye(size(P_pred)) - K * H) * P_pred;
end