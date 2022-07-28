function Snippets = GetSnippets(r,s,p,l)
    onset = find(s>0);
    if isempty(onset)
        Snippets = nan;
        error('Can''t predict anything on this one, no occurence.')
    end
    
    % in case r is a matrix / very general but not optimal..
    si = size(r);
    
    % find the 't' dimension
    tdim = find(si == length(s));
    if isempty(tdim)
        error('r and s dimensions don''t match...')
    end
    
    % put it first
    rper = permute(r,[tdim setdiff(1:length(si),tdim)]);
    
    S = NaN(length(onset),l,prod(si(setdiff(1:length(si),tdim))));
    for jj=1:length(onset)
        tmp = rper(max(onset(jj)-p+1,1):min(onset(jj)+l-p,size(rper,1)),:);
        if size(tmp,1) < l
            if onset(jj)-p < 1
                S(jj,end-size(tmp,1)+1:end,:) = tmp;
            elseif onset(jj)+l-p > size(rper,1)
                S(jj,1:size(tmp,1),:) = tmp;
            end
        else
            S(jj,1:size(tmp,1),:) = tmp;
        end
    end
    
    Snippets = reshape(S,[size(S,1),size(S,2),si(setdiff(1:length(si),tdim))]);
end

        
%         if p == 2;
%             % Purely predictive
%             tmp = rper(max(onset(jj)-l,1):onset(jj)-1,:);
%             S(jj,end-size(tmp,1)+1:end,:) = tmp;
%         elseif p == 1;
%             % Mix: center window on stimulus itself
%             tmp = rper(max(onset(jj)-l/2,1):min(onset(jj)+l/2,size(rper,1))-1,:);
%             if size(tmp,1) < l
%                 if onset(jj)-l/2 < 1
%                     S(jj,end-size(tmp,1)+1:end,:) = tmp;
%                 elseif onset(jj)+l/2 > size(rper,1)
%                     S(jj,1:size(tmp,1),:) = tmp;
%                 end
%             else
%                 S(jj,1:size(tmp,1),:) = tmp;
%             end
%         elseif p == 0;
%             % Not a predictive filter
%             tmp = rper(onset(jj):min(onset(jj)+l,size(rper,1))-1,:);
%             S(jj,1:size(tmp,1),:) = tmp;
%         end