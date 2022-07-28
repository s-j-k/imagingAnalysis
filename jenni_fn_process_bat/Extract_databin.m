function Extract_databin(filename,channel,plane)
cd 
for i=1:nPlanes
    cd([suite2ppath 'plane' num2str(i-1)]);
    data = load('Fall.mat');
    ly = data.ops.Ly;
    lx = data.ops.Lx;
    if strcmp(channel,'green')
        fileID = fopen('data.bin','r'); % open green channel binary file
    elseif strcmp(channel,'red')
        fileID = fopen('data_chan2.bin','r'); % open green channel binary file
    end   
    for j = 1:nFiles
        k=0; a=1;
        nimg = nFrames_oneplane(j+1,i); %10000 for higher frame rate
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
        avg_sessions{j,i} = mean(avgA,3)';
        if numel(find(isnan(mean(avgA,3))))>0
          disp('')  
        end
    end   
    % Save the mean img of this plane   
    tosave = avg_sessions(:,i);
    if strcmp(channel,'green')
        save([mouse '_MeanGreenImgPerSessions_Plane' num2str(i) '.mat'],'tosave');
    elseif strcmp(channel,'red')
        save([mouse '_MeanRedImgPerSessions_Plane' num2str(i) '.mat'],'tosave');
    end
        
    fclose all;  
    % Plot mean image each sessions if asked
    if imshow
        try
            figure; 
            for j=1:4
                subplot(2,2,j);
                imagesc(avg_sessions{j,i});colormap gray; 
                title(num2str(j));
            end
        catch
            continue;
        end
    end

end