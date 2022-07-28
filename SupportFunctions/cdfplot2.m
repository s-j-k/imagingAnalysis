function [xCDF,yCDF,handleCDF,stats] = cdfplot2(x)

% Get sample cdf, display error message if any
[yy,xx,~,~,eid] = cdfcalc(x);
if isequal(eid,'VectorRequired')
    error(message('stats:cdfplot:VectorRequired'));
elseif isequal(eid,'NotEnoughData')
    error(message('stats:cdfplot:NotEnoughData'));
end

% Create vectors for plotting
k = length(xx);
n = reshape(repmat(1:k, 2, 1), 2*k, 1);
xCDF    = [-Inf; xx(n); Inf]/17;
yCDF    = [0; 0; yy(1+n)];

%
% Now plot the sample (empirical) CDF staircase.
%

hCDF = plot(xCDF , yCDF);
if (nargout>0), handleCDF=hCDF; end
% grid  ('on')
xlabel('Area under the tuning curve (normalized)')
ylabel('% of cells')
% title (getString(message('stats:cdfplot:Title')))

%
% Compute summary statistics if requested.
%

if nargout > 1
   stats.min    =  min(x);
   stats.max    =  max(x);
   stats.mean   =  mean(x);
   stats.median =  median(x);
   stats.std    =  std(x);
   
end