function newtonecolormap=colormaptonotopy(nbtones)
colormaptone = jet(nbtones);
for tid = 1:nbtones
    ToneColorMap = colormaptone(tid,:);%[linspace(1, 0, 124), (zeros(1, 132))];
    if ToneColorMap(:,1)==1
        logvalues1 = ones(1,257)';
    elseif  ToneColorMap(:,1)==0
        logvalues1 = (log([(ToneColorMap(:,1)+0.1)*nbtones:((nbtones-1)/256):nbtones])/log(nbtones))';
    else
        logvalues1 = (log([linspace(ToneColorMap(:,1)*nbtones,nbtones,257)])/log(nbtones))';
    end
    if ToneColorMap(:,2)==1
        logvalues2 = ones(1,257)';
    elseif  ToneColorMap(:,2)==0
        logvalues2 = (log([(ToneColorMap(:,2)+0.1)*nbtones:((nbtones-1)/256):nbtones])/log(nbtones))';
    else
        logvalues2 = (log([linspace(ToneColorMap(:,2)*nbtones,nbtones,257)])/log(nbtones))';
    end
    if ToneColorMap(:,3)==1
        logvalues3 = ones(1,257)';
    elseif  ToneColorMap(:,3)==0
        logvalues3 = (log([(ToneColorMap(:,3)+0.1)*nbtones:((nbtones-1)/256):nbtones])/log(nbtones))';
    else
        logvalues3 = (log([linspace(ToneColorMap(:,3)*nbtones,nbtones,257)])/log(nbtones))';
    end
    colorMap = flipud([logvalues1,logvalues2,logvalues3]);
    cmap{tid} = colorMap;
end
newtonecolormap = cell2mat(cellfun(@(x) x(end,:),cmap,'UniformOutput',0)');
end