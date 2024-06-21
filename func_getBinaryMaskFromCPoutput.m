% This function converts binary masks outputed by cellprofiler to a matlab
% mat file
% Input:
%   tiffpath - the path which the Cell Profiler tiff masks are stored at
%   xr_padding - the number of pixels to be cropped on the right side
% Output:
%   A maskdat.mat is generated under the tiffpath containing a cell array
%   with size 1 x num_masks where each element is a logical array with size
%   mask_height x (mask_width + xr_padding)
% Example function call:
%   func_getBinaryMaskFromCPoutput('tiffpath', [path], 'xr_padding', [int])
%   masksdat = func_getBinaryMaskFromCPoutput('tiffpath', [path], 'xr_padding', [int])

function masks = func_getBinaryMaskFromCPoutput(varargin)
% Parse tiffpath and xr_padding input
p = inputParser;
p.KeepUnmatched = true;
p.addParameter('tiffpath', pwd)
p.addParameter('xr_padding',0)
% p.addParameter('yd_padding',0)
p.parse(varargin{:});

% outputPath is an dir object of the names of the tiffmasks
outputPath = dir([p.Results.tiffpath '/*.tiff']);
xr_exten = p.Results.xr_padding;
% yd_exten = p.Results.yd_padding;
disp('arguments parsed')
% cd to tiffpath and get all tiff mask names
cd(p.Results.tiffpath)
filenames = {outputPath.name};
masks = cell(1,length(filenames));

% loop through all tiff mask names; let tiffname be the current tiff mask
% name
for i = 1:length(filenames); tiffname = filenames{i};
    % read tiff mask and convert it to logical array
    mask = imread(tiffname);
    mask = logical(mask);
    % add false padding to the right
    xr_padding = false(size(mask,1), xr_exten);
    mask = [mask xr_padding];
    % store it in masks array
    masks{i} = mask;
end
% save the masks array to maskdat.mat
save('maskdat.mat', 'masks','-v7.3');
disp(['maskdat saved at' p.Results.tiffpath]);
end