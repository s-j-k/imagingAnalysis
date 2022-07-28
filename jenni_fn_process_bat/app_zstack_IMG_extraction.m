function app_zstack_IMG_extraction()
animal = 'Gray2';
datatable = readtable(['K:\Jenni\' animal '\Zstackinfopersite.xlsx']);
for st = [1,2,4,5,6]%1:size(([datatable.Site]'),2)
clear slice_green
    if strcmp(datatable.CorrespondingZstack{st},'NA')
    else
        cd([datatable.datadrive{st} ':\Jenni\' animal '\'])
        % d = dir('*.sbx');
        fn = [datatable.CorrespondingZstack{st}];%strtok(d.name,'.');
        sbxread(fn,1,1);            % read one frame to read the header of the image sequence
        global info;
        A = sbxread(fn,1,info.max_idx);
        framnb = datatable.NbFramesPerSlice(st);
        indframe = 1;
        for slicenb = 1:size(1:framnb:((size(A,4)+1)),2)-2%info.otparam(3)
            indframe = indframe+framnb;
            slice_green{slicenb} = mean(squeeze(A(1,:,:,[indframe:(indframe+framnb)-1])),3);%(A(1,:,:,[slicenb:framnb:info.max_idx])),3);
            % slice_gray{slicenb} = mean(squeeze(A(2,:,:,[indframe:(indframe+framnb)-1])),3);%(A(2,:,:,[slicenb:framnb:info.max_idx])),3);
            if slicenb ==1
                imwrite(uint16(slice_green{slicenb}),[datatable.datadrive{st} ':\Jenni\' animal '\z-stack\' fn,'_green.tif'])
                % imwrite(uint16(slice_gray{slicenb}),'test_gray.tif')
            else
                imwrite(uint16(slice_green{slicenb}),[datatable.datadrive{st} ':\Jenni\' animal '\z-stack\' fn,'_green.tif'],'WriteMode','append')
                % imwrite(uint16(slice_gray{slicenb}),'test_gray.tif','WriteMode','append')
            end
        end
        save([datatable.datadrive{st} ':\Jenni\' animal '\z-stack\' fn,'_green.mat'],'slice_green')
        % save([datatable.datadrive{st} ':\Jenni\' animal '\z-stack\' fn,'_green.tif'],'slice_green')
    end
end
end

%% save app mean image in tiff
% for nbp = 1:2
% load(['W:\Aaron\aw017\Site1\suite2p\plane' num2str(nbp-1) '\Fall.mat'])
% imwrite(uint16(ops.meanImg),['GreenSite1plane' num2str(nbp-1) '.tif'])
% imwrite(uint16(ops.meanImg_chan2),['GraySite1plane' num2str(nbp-1) '.tif'])
% end
% %
% cells = find(iscell==1)
%
% for cc = 1:length(cells)
% coord(cc,1) = (mean(stat{(cells(cc))}.xpix)/2)+214;
% coord(cc,2) = (mean(stat{(cells(cc))}.ypix)/2)+158;
% end