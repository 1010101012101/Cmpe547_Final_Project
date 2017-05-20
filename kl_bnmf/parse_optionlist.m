function [varargout] = parse_optionlist(defaults, varargin)
% PARSE_OPTIONLIST		Parses Parameter-Value Pairs
%
%  [VARARGOUT] = PARSE_OPTIONLIST(DEFAULTS, <Parameter>, <Value>, ...)
%
% Inputs :
%	:
%
% Outputs :
%	:
%
% Usage Example : [] =
% parse_optionlist(parse_optionlist({'a',1;'b','2';'c',3}, varargin{:});
%
%
% Note	:
% See also

% Uses :

% Change History :
% Date		Time		Prog	Note
% 10-May-2000	 1:05 PM	ATC	Created under MATLAB 5.3.1.29215a (R11.1)

% ATC = Ali Taylan Cemgil,
% SNN - University of Nijmegen, Department of Medical Physics and Biophysics
% e-mail : cemgil@mbfys.kun.nl 

names =  defaults(:,1);
varargout = defaults(:,2);
  
for i=1:2:length(varargin),
  idx = strmatch(varargin{i},names,'exact');
  if isempty(idx), error(['Unknown Option : ' varargin{i}]); end;
  varargout{idx} =  varargin{i+1};
end;