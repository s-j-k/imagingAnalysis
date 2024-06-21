sbxpath='Y:\sjk\202207\z.sk73\20220718\000_004';
targetPath='Y:\sjk\202207\z.sk73\20220718\000_004';
cd(sbxpath)
files=dir([sbxpath '\*.sbx']);

for jj = 1:length(files)
    fname=files(jj).name;
    disp(['processing' fname]); tic;
    sbxname = strsplit(fname,'.');
    sbxname=cell2mat(sbxname(1));
    sbx2tif(sbxname);
    disp('done');
end
