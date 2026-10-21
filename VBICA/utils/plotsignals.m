function plotsignals(data,columns,SCALE);
% plotsignals(data,columns,SCALE);
%
% Plots the L rows of data in L subfigures organised
% in COLUMNS columns.
%
%
% -----------
% Input
% -----------
%
% Necessary parameters
%
% data  LxT dataset
%
%
% Optional parameters
%
% columns Number of subfigure columns
% SCALE   1: Plot all subfigures to the same scale (Default)
%         0:   Plot all subfigures to their own scale
%
%
% --------------------------------------------------------------
%
% Original code by Rizwan Choudrey 
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)
  
  
[sources points] = size(data);


if nargin<2
    columns = 2;
end
if nargin<3
  SCALE=1;
end

rows = ceil(sources/columns);

minx = 0;
maxx = points;
miny = floor(min(min(data)));
maxy = ceil(max(max(data)));

for i=1:sources
  subplot(rows,columns,i);
  plot(data(i,:));
  if(SCALE)
  axis([minx maxx miny maxy]);
  end
end
