function ExtractMeanImage()

animal = 'RoofBuddy2';
datatable = readtable(['N:\Jenni\' animal '\ProcessingProgress.xlsx']);

manualroisdone = find(table2array(datatable(:,13))==1);%find(~cellfun(@isempty,cellfun(@(x) find(x=='1'),table2array(datatable(:,13)),'UniformOutput',0)));%find(table2array(datatable(:,13))==1);
daytraining = table2array(datatable(:,1));
siteind = table2array(datatable(:,5));
preprocessedrois = unique(siteind(manualroisdone));
framerate = table2array(datatable(:,8));
protocol = table2array(datatable(:,7));
path = table2array(datatable(:,9));
behaviorfile = table2array(datatable(:,10));
tonotopyfilesind = find(~cellfun(@isempty,(cellfun(@(x) strfind(x,'tonotopy')>0,behaviorfile,'UniformOutput',0))));
sweepfilesind = find(~cellfun(@isempty,(cellfun(@(x) strfind(x,'LinearSweep')>0,behaviorfile,'UniformOutput',0))));
angiefilesind = find(~cellfun(@isempty,(cellfun(@(x) strfind(x,'Angie')>0,behaviorfile,'UniformOutput',0))));
echofilesind = find(~cellfun(@isempty,(cellfun(@(x) strfind(x,'Echo')>0,behaviorfile,'UniformOutput',0))));
melfilesind =  find(~cellfun(@isempty,(cellfun(@(x) strfind(x,'Mel')>0,behaviorfile,'UniformOutput',0))));
datadrivecolumn = (table2array(datatable(:,15)));

uniquesite = unique(siteind);

for sitenb = 22%:length(uniquesite)%siteind
     DataDrive = datadrivecolumn{find(siteind==uniquesite(sitenb),1)};
      
    clear frameratesite
    frameratesite = unique(framerate(find(siteind==uniquesite(sitenb))));
    frameratesite(isnan(frameratesite)) = [];
    cd([DataDrive ':\Jenni\' animal '\' path{find(siteind==uniquesite(sitenb),1)} '\31.25\suite2p\plane0'])
    data = load('Fall.mat');
    
    for fr = 1:length(frameratesite)
        if frameratesite(fr) == 31.25
            ly = data.ops.Ly;
            lx = data.ops.Lx;
        elseif frameratesite(fr) == 62.5
            ly = data.ops.Ly/2;
            lx = data.ops.Lx;
        elseif frameratesite(fr) == 125
            ly = data.ops.Ly/4;
            lx = data.ops.Lx;
        end
        cd([DataDrive ':\Jenni\' animal '\' path{find(siteind==uniquesite(sitenb),1)} '\' num2str(frameratesite(fr)) '\suite2p\plane0'])
        fileID = [];
        fileID = fopen('data.bin','r'); % open green channel binary file
        k=0; a=1;
        nimg = 4000;%nFrames_oneplane(j+1,i); %10000 for higher frame rate
        blksize = 2000;%2000; % nb of frames loaded at a time (depend on RAM)
        to_read = min(blksize,nimg-k);
        avgA = [];
        %         avgA = nan(lx,ly);
        while to_read>0
            clear A
            A = fread(fileID,ly*lx*to_read,'*int16');
            A = reshape(A,lx,ly,[]);
            avgA(:,:,a) = mean(A,3);
            a=a+1;
            k = k+to_read;
            to_read = min(blksize,nimg-k);
        end
        % to check the bloc averages:   figure;for l=1:9,subplot(3,3,l);hold on;imagesc(avgA(:,:,l));colormap gray;title(num2str(l));end
        avg_sessions{fr} = mean(avgA,3)';
        
        if frameratesite(fr) == 31.25
        elseif frameratesite(fr) == 62.5
            avg_sessions{fr} = [zeros((ly/2),lx); mean(avgA,3)';zeros((ly/2),lx)];
        elseif frameratesite(fr) == 125
            avg_sessions{fr} = [zeros((ly/2)*3,lx); mean(avgA,3)';zeros((ly/2)*3,lx)];
        end
        
        if numel(find(isnan(mean(avgA,3))))>0
            disp('')
        end
        close all
        h = figure(fr); hold on;
        imagesc(flipud(avg_sessions{fr})); colormap gray; axis off
        axis tight
        %         imwrite(h,['N:\Jenni\RoofBuddy2\000_h5\site' num2str(uniquesite(sitenb)) '\MeanImg' num2str(frameratesite(fr)) '.tif'])
        imwrite(uint16((avg_sessions{fr})),[DataDrive ':\Jenni\' animal '\' path{find(siteind==uniquesite(sitenb),1)} '\MeanImg' num2str(frameratesite(fr)) '.tiff']);
        %         disk = 'E:';
        %         javaaddpath ([disk '\KishoreLab\Shared\Matlab\preprocessing\MIJI\mij.jar'])
%         javaaddpath ([disk '\KishoreLab\Shared\Matlab\preprocessing\MIJI\ij-1.52i.jar'])
%         MIJ.start([disk '\KishoreLab\Shared\Matlab\Fiji.app'])
%         MIJ.createImage((uint16(flipud(avg_sessions{fr}))));
    end
end
end