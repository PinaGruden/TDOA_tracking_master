function [wn,mn,Pn,Tagn] = gmphd_limit(w,m,P,Jmax,Tag)

if numel(w) > Jmax
    [~,ranking] = sort(w,'descend'); %gives indices of ranking in descending order
    indx=ranking(1:Jmax); %indices that we want to keep 
    wn = w(indx); 
    wn=wn*(sum(w)/sum(wn));%re-weight the components
    mn=m(:,indx);
    Pn= P(:,:,indx);
    Tagn=Tag(indx);
else
    wn=w;
    mn=m;
    Pn=P;
    Tagn=Tag;
end

end

