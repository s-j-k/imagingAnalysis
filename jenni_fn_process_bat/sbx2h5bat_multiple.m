cd L:\su\DATA\imagingData\passive\sk023\day3\%K:\Jenni\Gray2\
h5files = dir('*.sbx');
fnames = cellfun(@(x) x(1:end-4),{h5files.name},'UniformOutput',0);%{'pink7_002_012'};%};%,'JL009_007_004'};%,'fz006_009_002'};

nfiles = length(fnames);

for i=1:nfiles
    if exist([fnames{i} '.h5'])>0
    else
    sbx2h5bat(fnames{i});
    end
end
