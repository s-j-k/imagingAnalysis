function b = redwhite(m)
%REDBLUE   Bilinear red/blue color map.
%
%   REDBLUE(M) returns an M-by-3 matrix containing the colormap.
%   REDBLUE, by itself, is the same length as the current colormap.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redblue)
%
%
%   See also BROWNBLUE, YARG, YARGPRINT, GRAYPRINT.
%   See also HSV, HOT, COOL, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.

if nargin < 1, m = size(get(gcf,'colormap'),1); end
top = floor(m/2);
bot = m-top;

btop = [1.0*ones(top,1),([1:top]'/top).^2,ones(top,1)];   % 1.0 hue good rebot = [1*ones(bot,1),ones(bot,1),ones(bot,1)]; % 0.7 hue good blue
b=hsv2rgb([bbot;btop]);