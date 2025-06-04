% ---------- annotateSpecies --------------------------------------------
function hTxt = annotateSpecies(V,metric,n,ax)
% annotateSpecies(V,metric,n,ax) - label the n highest-metric points
    if nargin<4||isempty(ax)
        ax = gca;
    end
    
    if nargin<3||isempty(n)
        n = size(V,1);
    end

    [~,idx] = maxk(metric(:),n);
    hold(ax,'on');
    hTxt = text(ax,V(idx,1),V(idx,2),V(idx,3), compose('S%d',idx), ...
        'FontWeight','bold','FontSize',9,'Color','k');
end
