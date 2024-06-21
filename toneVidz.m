  
cd(suite2ppath);
TC = cell(nPlanes,1);
for i=1:nPlanes
    cd([suite2ppath '\plane' num2str(i-1)]);    
    data = load('Fall.mat');
    ly = data.ops.Ly;
    lx = data.ops.Lx;
    fileID = fopen(readChannel,'r'); % open binary file, 
    for j=1:nFiles        
            k=0;
            nimg = nFrames_oneplane(j+1,i);
            % here is where you specify which frames based on the tone 
            blksize = 2000;%7000; % nb of frames loaded at a time 
            to_read = min(blksize,nimg-k);  
            while to_read>0
                a = fread(fileID,ly*lx*to_read,'*int16');
                A = reshape(a,lx,ly,[]); 
                frame = im2frame(uint16(YourBinaryImage), bw);
    %             avgA = mean(A,3); % figure;subplot(2,2,1);imagesc(avgA);colormap gray;
                currFrame = getframe(fig);
                writeVideo(vidObj,currFrame);    
                fun = @(block_struct) nanmean(block_struct.data(:));
                b = blockproc(A(:,:,1),win,fun);            
                tempTC = nan(size(b,1)*size(b,2),size(A,3));
                tempTC(:,1) = b(:);
                for kk=2:size(A,3)
                    b = blockproc(A(:,:,kk),win,fun);
                    tempTC(:,kk) = b(:);
                end
                TC{i} = [TC{i} tempTC]; % figure;plot(mean(TC{1})');
                k = k+to_read;
                to_read = min(blksize,nimg-k);
            end
            disp(['File ' num2str(j) '/' num2str(nFiles) ' DONE! ' num2str(size(TC{i},1)) ' frames so far.']);
        end