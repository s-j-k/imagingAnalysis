function fn_getTuningPopulation(tuning,popTuning,nNeuron,nPlanes, neuronEachPlane,neuronPlane,toneLabel,toneindex,refImg,roisBound,savePath)
% fn_getTuningPopulation - Description
%
% Syntax: fn_getTuningPopulation(input)
%
% Long description

nTones = length(toneLabel);%2/20/24 removed 3 for pure tones only
nTicks = 5; tickLocation = round(linspace(1,nTones,nTicks));

figure; subplot(2,2,1)
percentResp = [sum(tuning.responsiveCellFlag), sum(~tuning.responsiveCellFlag)] ./ nNeuron;
pie([sum(tuning.responsiveCellFlag), sum(~tuning.responsiveCellFlag)],...
    {['resp ' int2str(sum(tuning.responsiveCellFlag))],['not resp ' int2str(sum(~tuning.responsiveCellFlag))]});
title([int2str(percentResp(1)*100) '% cell Responsive'])

subplot(2,2,2)
plot(1:nTones, popTuning.responsiveCount,'LineWidth',2);
set(gca, 'XTick', tickLocation); set(gca, 'XTickLabel', toneLabel(toneindex(tickLocation)));
xlim([1 nTones]); xlabel('Frequency (kHz)'); ylabel('Cell Count')
title('Significant Frequency Count of all Neurons')

subplot(2,2,3)
%errorbar(freqAxis, popTuning, popTuningStd,'LineWidth',2);
plot(1:nTones, nanmean(popTuning.median,2),'LineWidth',2); hold on; 
plot(1:nTones, nanmean(popTuning.mean,2),'LineWidth',2);
set(gca, 'XTick', tickLocation); set(gca, 'XTickLabel', toneLabel(toneindex(tickLocation)));
xlim([1 nTones]); xlabel('Frequency (kHz)'); ylabel('Cell Count')
ylabel('dF/F'); legend('Median','Mean')
title('Mean Tuning Activity')

subplot(2,2,4)
% here also we changed bfcountMedian and bfcountMean to index 1:17
plot(1:nTones, popTuning.bfCountMedian,'LineWidth',2); hold on; plot(1:nTones, popTuning.bfCountMean,'LineWidth',2); 
set(gca, 'XTick', tickLocation); set(gca, 'XTickLabel', toneLabel(toneindex(tickLocation)));
xlim([1 nTones]); xlabel('Frequency (kHz)'); ylabel('Cell Count')
title('BF Count of Responsive Neurons'); legend('Median','Mean')
saveas(gcf,[ savePath ...
        '/population/populationTuning.png']);

%---------PLOT SIGNICANT TEST RESULTS ON POPULATION LEVEL-----------

figure;
for i = 1:nPlanes

    C = colormap('jet'); colormapIndex = round(linspace(1,size(C,1),nTones));

    rotatedrefImg=imrotate(refImg{i},-180);  
    rotatedrefImg=flip(rotatedrefImg,1);
    rotatedrefImg=flip(rotatedrefImg,2);
    subplot(2,2,i);           
    imagesc(rotatedrefImg);colormap gray;hold on;
    ylim([0 size(rotatedrefImg,1)]);xlim([0 size(rotatedrefImg,2)]);
    xticks([]);yticks([]);title(['Plane' int2str(i)])
    
    subplot(2,2,2+i);             
    imagesc(rotatedrefImg);colormap gray;hold on;
    ylim([0 size(rotatedrefImg,1)]);xlim([0 size(rotatedrefImg,2)]);
    for j = 1:neuronEachPlane(i)
        cellIndex = neuronPlane(i) + j;
        if ~isempty(roisBound{i}{j})
            x = roisBound{i}{j}(:,1); %freehand rois have the outlines in x-y coordinates
            y = roisBound{i}{j}(:,2); %matlab matrices are inverted so x values are the 2nd column of mn coordinates, and y is the 1st columna
            if tuning.responsiveCellFlag(cellIndex)
                patch(x,y,C(colormapIndex(popTuning.bfMedian(cellIndex)),:),'EdgeColor','none');
                % flip horizontally & rotate 90 degrees to the right
            else
                patch(x,y,[0.8 0.8 0.8],'EdgeColor','none');
            end
        end
    end
    xticks([]);yticks([]);title(['Plane' int2str(i)])
end
saveas(gcf,[savePath...
        '/population/tuningMap.png']);
    
end