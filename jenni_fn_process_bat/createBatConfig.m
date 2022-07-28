function createBatConfig

mouse = 'Gray2';
Rootaddress = ['K:\Jenni\' mouse '\004_h5\site5\31.25']
roiName =  ['roi.zip'];%[Rootaddress, '\roi.zip'];
tcName =  ['TC\' mouse '_TC_plane0.mat TC\'];
spikeFile = ['TC\' mouse '_spike_plane0.mat TC\'];

T = func_createMouseConfig(mouse,...
    'savepath',' E:\KishoreLab\Jenni\process2pRoiTracking\config\mouse',...
    'h5path',[Rootaddress],'sbxpath',[Rootaddress(1:end-19)],... % 20 for RoofBuddy1
    'suite2ppath',[Rootaddress '\suite2p'],...
    'datapath',Rootaddress,...
    'behavpath',[Rootaddress(1:end-13) '\behavior'],...
    'roiFile',roiName,'tcFile',tcName,'spikeFile',spikeFile,'nPlanes',1);
end%