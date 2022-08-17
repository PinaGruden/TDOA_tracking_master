function [qz_update,m_update,P_update,Tag_update] = kalman_update_multiple(z,model,m,P,Tag)

plength= size(m,2);
zlength= size(z,2);
x_dim=size(m,1);


qz_update= zeros(plength,zlength);
m_update = zeros(x_dim,plength,zlength);
P_update = zeros(x_dim,x_dim,plength);

for id=1:plength
        [qz_temp,m_temp,P_temp] = kalman_update_single(z,model.H,model.R,m(:,id),P(:,:,id));
        qz_update(id,:)   = qz_temp;
        m_update(:,id,:) = m_temp;
        P_update(:,:,id) = P_temp;        
end

Tag_update = repmat(Tag,zlength,1);

end



