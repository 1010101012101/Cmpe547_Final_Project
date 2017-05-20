function marginal_plot(x, t, v, r1, cx)
% MARGINAL_PLOT		
%
%  [] = marginal_plot(x, t, v, r, cx)
%
% Inputs :
%	x : array of size M_by_N
%          t,v : arrays of size M_by_I and I_by_N
%          r : Relative size of the main panel 0<r<1
%          cx : color axis
%
% Outputs :
%	:
%
% Usage Example : [] = marginal_plot();
%
%
% Note	:
% See also

% Uses :

% Change History :
% Date		Time		Prog	Note
% 08-Nov-2007	12:39 PM	ATC	Created under MATLAB 6.5.0 (R13)

% ATC = Ali Taylan Cemgil,
% SPCL - Signal Processing and Communications Lab., University of Cambridge, Department of Engineering
% e-mail : atc27@cam.ac.uk

%clf
pos = get(gca, 'pos');
delete(gca)

if nargin<4,
    r1 = 0.9;
end;
r2 = 1 - r1;
r3 = 0.8*r2;

if ~exist('cx', 'var'), cx = []; end;


axd = axes('position',[pos(1) pos(2)+r2*pos(4) pos(3)*r3 pos(4)*r1]);
imagesc(t);
set(axd,'ydir','normal','xtick', [], 'ytick', []);
if ~isempty(cx),
    caxis(cx);
end;

axt = axes('position',[pos(1)+r2*pos(3) pos(2) pos(3)*r1 pos(4)*r3]);
%set(axt,'vis','off','ydir','normal');
imagesc(v)
set(axt,'ydir','normal','xtick', [], 'ytick', []);
if ~isempty(cx),
    caxis(cx);
end;

axtd = axes('pos',[pos(1)+r2*pos(3) pos(2)+r2*pos(4) pos(3)*r1 pos(4)*r1] );
imagesc(x);
set(axtd,'ydir','normal','xtick', [], 'ytick', []);

if ~isempty(cx),
    caxis(cx);
end;