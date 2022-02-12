function alldff=summary_all_cells(datapath,counter,filename)
% this first part takes all of the cells matched from run_pipeline_2 output
% between one TSeries and the template TSeries 
% change directory to location of TSeries .mat files
% cd 'C:\Users\bbscott\Dropbox (BOSTON UNIVERSITY)\Kinematic_Data\Bruker\processed\mouse'
cd(datapath)
datafolders = dir(datapath);
datafolders=struct2cell(datafolders)'; newfolders=datafolders(4:13,1:3);
clear datafolders

remove = [];
r1 = 1000;r2 = 19000;
time = (0:(r2-r1))/30; % 30hz is sampling rate, 18000 is number of frames
% for uu = 1:length(filename)
%     currentname = filename;
    load(filename,'ISO_F','N','A_keep','M_green_final');
% 
%         % get microns per pixel -- this section of code is not working at
%         the moment due to errors with omeMeta 
%         cd(rawpath); rawfolders=dir('TSeries*');
%         cd(rawfolders(1).name); %go to the first folder
%         rawfiles=dir('TSeries*');
%         data = bfopen(rawfiles(1).name); %Read the header
%         omeMeta = data{1, 4};%Extract metadata
%         mpp=omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER);
%         microns_per_pixel=double(mpp);
%         clear T Y mpp data omeMeta data Y1
    % remove inactive cells
     numcells = N; % from run_pipeline_2
     cells=ones(numcells,1); alldff{1,1} = 'Date'; alldff{1,2} = 'Raw Data path'; 
     alldff{1,3}='CellNum'; alldff{1,4} = 'Cell masks';
     alldff{2,1} = filename(10:16); alldff{2,2} = [newfolders(counter,2) '\' filename];
     alldff{1,5} = 'ISO_F'; alldff{2,5} = ISO_F(:,:);
     % in this array each t series should have data in either column 4 or
     % column 5 not both
    cellnum=[]; %#ok<NASGU>
    [ISO_Fcells,cellnum] = cell_finder(A_keep,N,M_green_final,time);
    % this is where the cell masks are binarized and generated
    % cellnum is index of the cells, use this in A_keep to get cell masks
    
    % now threshold data; find index of cells that spike
     for nn = 1:length(cellnum)
         if max(ISO_Fcells(nn,r1:r2)) < mean(ISO_Fcells(nn,r1:r2))+0.75
             remove(end+1) = nn;
         end
         cells(remove)=0;
         [realcellid,ccol]=find(ISO_Fcells(:,1)); % index of real cells
     end
     clear ccol % replacing ccol with ~ sometimes does not work? fix later
     
     % then pull the entire timetrace for those cells and
     % put that into your alldff cell array
     
     for qq=1:length(realcellid)
         alldff{qq+2,3} = cellnum(realcellid(qq)); % take the cell number of cells that spike (realcellid) from cellnum
         % get masks of just the cells you want
         mask1=reshape(A_keep(:,realcellid(qq)),256,256);
         mask2=full(mask1);
         alldff{qq+2,4} = mask2;
         alldff{qq+2,5} = ISO_F(realcellid(qq),:); 
         % this is the corresponding template cell trace
         alldff{qq+2,6} = ISO_Fcells(realcellid(qq),r1:r2);
         clear mask1 mask2
     end
        clear realcellid qq cellnum
end
% %%
%  % now compare 
% %  subplot(2,1,1)
% % time=(0:18000)/30;plot(time,alldff{12,4});
% % title('Single Cell (920nm) 01/23/21 9:30am')
% % ylabel('DF/F')
% % xlabel('Seconds')
% % subplot(2,1,2)
% % time=(0:18000)/30;plot(time,alldff{12,5});
% figure(2)
% n=25;time=(0:18000)/30;
% plot(time,ISO_F(n,1000:19000));
% % plot the data
% %%
% figure;
% subplot(2,3,2)
% plot(time,ISO_F(10,1000:19000))
% ylim([-.25 7])
% xlim([0 18000/30])
% title('Cell 10 (11:30am)')
% ylabel('DF/F')
% xlabel('Seconds')
% 
% subplot(2,3,1)
% plot(time,F_dff(10,1000:19000))
% xlim([0 18000/30])
% ylim([-.25 7])
% title('Cell 10 (9:30am)')
% ylabel('DF/F')
% xlabel('Seconds')
% 
% subplot(2,3,5)
% imagesc(ISO_Fcells(realcellid,1000:19000),[0 1])
% title('Population (920nm 11:30am)')
% ylabel('Cell')
% xlabel('Seconds')
% 
% set(gca,'ytick',[1 length(realcellid)])
% % set(gca,'xtick',[1 20*30 40*30 1800])
% % set(gca,'xticklabel',{'0';'20';'40'; '60'})
% set(gca,'tickdir','out')
% 
% subplot(2,3,4)
% imagesc(F_dff(realcellid,1000:19000),[0 1])
% ylabel('Cell')
% xlabel('Seconds')
% set(gca,'ytick',[1 length(realcellid)])
% % set(gca,'xtick',[1 20*30 40*30 1800])
% % set(gca,'xticklabel',{'0';'20';'40'; '60'})
% set(gca,'tickdir','out')
% title('Population (920nm 9:30am)')
% 
% subplot(2,3,6)
% colorbar
% box off
% axis off
% %%
% ss=10;
% figure;
% subplot(2,3,2)
% plot(time,alldff{ss,5})
% ylim([-.25 7])
% xlim([0 18000/30])
% title('Cell 10 (11:30am)')
% ylabel('DF/F')
% xlabel('Seconds')
% 
% subplot(2,3,1)
% plot(time,alldff{ss,4})
% xlim([0 18000/30])
% ylim([-.25 7])
% title('Cell 10 (9:30am)')
% ylabel('DF/F')
% xlabel('Seconds')
% 
% isodff = alldff_old{2,5}; % 0123 @ 11:30am
% fdff = alldff_old2{2,5}; % 0124 @ 9:30am 
% tracesiso =isodff(realcellid,:);
% % tracesiso= cell2mat(tracesiso);
% tracesfdff =fdff(realcellid,:);
% subplot(2,3,5)
% imagesc(tracesiso(:,:),[0 1])
% title('Population (920nm 01/23 11:30am)')
% ylabel('Cell')
% xlabel('Seconds')
% 
% set(gca,'ytick',[1 length(realcellid)])
% % set(gca,'xtick',[1 20*30 40*30 1800])
% % set(gca,'xticklabel',{'0';'20';'40'; '60'})
% set(gca,'tickdir','out')
% 
% subplot(2,3,4)
% imagesc(tracesfdff(:,:),[0 1])
% ylabel('Cell')
% xlabel('Seconds')
% set(gca,'ytick',[1 length(realcellid)])
% % set(gca,'xtick',[1 20*30 40*30 1800])
% % set(gca,'xticklabel',{'0';'20';'40'; '60'})
% set(gca,'tickdir','out')
% title('Population (920nm 01/24 9:30am)')
% 
% subplot(2,3,6)
% colorbar
% box off
% axis off
% %%
% 
% 
% 
% %%
% % pull cells that fire
% ss=181;
% tracesfdff={};
% isodff={}; fdff_0124_2={};fdff_0124_4={};
% fdff_0125_2={};fdff_0125_4={};
% for ff = 3:ss
%     tracesfdff{ff-2,1} =alldff_0123_12{ff,4};
%     isodff{ff-2,1} = alldff_0123_12{ff,5}; % 0123 @ 11:30am
%     fdff_0124_2{ff-2,1} = alldff_0124_2{ff,5}; % 0124 @ 9:30am 
%     fdff_0124_4{ff-2,1} = alldff_0124_4{ff,5}; % 0124 @ 9:30am 
%     fdff_0125_2{ff-2,1} = alldff_0125_2{ff,5}; % 0124 @ 9:30am 
%     fdff_0125_4{ff-2,1} = alldff_0125_4{ff,5}; % 0124 @ 9:30am 
% % tracesiso= cell2mat(tracesiso);
% end
% % tracesfdff012312 =isodff(:,:);
% % tracesfdff01242 = fdff_0124_2(realcellid,:);
% % tracesfdff01244 = fdff_0124_4(realcellid,:);
% % tracesfdff01252 = fdff_0125_2(realcellid,:);
% % tracesfdff01254 = fdff_0125_4(realcellid,:);
% % tracesfdff012312 =isodff(realcellid,:);
% % tracesfdff01242 = fdff_0124_2(realcellid,:);
% % tracesfdff01244 = fdff_0124_4(realcellid,:);
% % tracesfdff01252 = fdff_0125_2(realcellid,:);
% % tracesfdff01254 = fdff_0125_4(realcellid,:);
% %%
% tracesfdff012312 =isodff(:,:);
% tracesfdff01242 = fdff_0124_2;
% tracesfdff01244 = fdff_0124_4;
% tracesfdff01252 = fdff_0125_2;
% tracesfdff01254 = fdff_0125_4;
% % tracesfdff012312 =isodff;
% tracesfdff01242 = fdff_0124_2;
% tracesfdff01244 = fdff_0124_4;
% tracesfdff01252 = fdff_0125_2;
% tracesfdff01254 = fdff_0125_4;
% %%
% tracesfdff012312 =isodff(:,:);
% tracesfdff012312 =tracesfdff012312';
% tracesfdff01242 = cell2mat(fdff_0124_2);
% tracesfdff01244 = cell2mat(fdff_0124_4);
% tracesfdff01252 = cell2mat(fdff_0125_2);
% tracesfdff01254 = cell2mat(fdff_0125_4);
% tracesfdff012312 = cell2mat(isodff);
% tracesfdff01242 = cell2mat(fdff_0124_2);
% tracesfdff01244 = cell2mat(fdff_0124_4);
% tracesfdff01252 = cell2mat(fdff_0125_2);
% tracesfdff01254 = cell2mat(fdff_0125_4);
% %%
% figure;
% subplot(2,4,1)
% cc=179;
% % plot(time,alldff{ss,5})
% % ylim([-.25 7])
% % xlim([0 18000/30])
% % title('Cell 10 (11:30am)')
% % ylabel('DF/F')
% % xlabel('Seconds')
% imagesc(tracesfdff012312(:,1:179),[0 1])
% title('Population (920nm 01/23 11:30am)')
% ylabel('Cell')
% xlabel('Seconds')
% set(gca,'ytick',[1 cc])
% % set(gca,'xtick',[1 20*30 40*30 1800])
% % set(gca,'xticklabel',{'0';'20';'40'; '60'})
% set(gca,'tickdir','out')
% 
% 
% 
% subplot(2,4,2)
% imagesc(tracesfdff01244(:,:),[0 1])
% title('Population (920nm 01/24 11:30am)')
% ylabel('Cell')
% xlabel('Seconds')
% set(gca,'ytick',[1 cc])
% % set(gca,'xtick',[1 20*30 40*30 1800])
% % set(gca,'xticklabel',{'0';'20';'40'; '60'})
% set(gca,'tickdir','out')
% 
% subplot(2,4,3)
% imagesc(tracesfdff01254(:,:),[0 1])
% title('Population (920nm 01/25 11:30am)')
% ylabel('Cell')
% xlabel('Seconds')
% set(gca,'ytick',[1 cc])
% % set(gca,'xtick',[1 20*30 40*30 1800])
% % set(gca,'xticklabel',{'0';'20';'40'; '60'})
% set(gca,'tickdir','out')
% 
% 
% subplot(2,4,5)
% imagesc(F_dff(realcellid,1000:19000),[0 1])
% ylabel('Cell')
% xlabel('Seconds')
% set(gca,'ytick',[1 cc])
% % set(gca,'xtick',[1 20*30 40*30 1800])
% % set(gca,'xticklabel',{'0';'20';'40'; '60'})
% set(gca,'tickdir','out')
% title('Population (920nm 01/23 9:30am)')
% 
% subplot(2,4,7)
% imagesc(tracesfdff01252(:,:),[0 1])
% ylabel('Cell')
% xlabel('Seconds')
% set(gca,'ytick',[1 cc])
% % set(gca,'xtick',[1 20*30 40*30 1800])
% % set(gca,'xticklabel',{'0';'20';'40'; '60'})
% set(gca,'tickdir','out')
% title('Population (920nm 01/25 9:30am)')
% 
% subplot(2,4,6)
% imagesc(tracesfdff01242(:,:),[0 1])
% ylabel('Cell')
% xlabel('Seconds')
% set(gca,'ytick',[1 cc])
% % set(gca,'xtick',[1 20*30 40*30 1800])
% % set(gca,'xticklabel',{'0';'20';'40'; '60'})
% set(gca,'tickdir','out')
% title('Population (920nm 01/24 9:30am)')
% 
% % subplot(2,4,4)
% % imagesc(tracesfdff01252(:,:),[0 1])
% % ylabel('Cell')
% % xlabel('Seconds')
% % set(gca,'ytick',[1 length(realcellid)])
% % % set(gca,'xtick',[1 20*30 40*30 1800])
% % % set(gca,'xticklabel',{'0';'20';'40'; '60'})
% % set(gca,'tickdir','out')
% % title('Population (920nm 01/24 9:30am)')
% 
% subplot(2,4,8)
% colorbar
% box off
% axis off

