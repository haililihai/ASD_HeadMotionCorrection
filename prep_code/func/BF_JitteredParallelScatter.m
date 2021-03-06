function [ff,xx] = BF_JitteredParallelScatter(dataCell,addMeans,doveTail,makeFigure,extraParams)
% Plots a scatter of a set of distributions with data offset randomly in x
% input is a cell with each element containing a collection of data.
%
%---INPUTS:
% dataCell, a cell where each element is a vector of numbers specifying a distribution
% addMeans, a flag for whether to add a strip for the mean of each group
% doveTail, whether to show a kernel smoothed distributions
% makeFigure, a flag to specify whether to generate a new figure to plot into
% extraParams, a structure of additional custom parameters
%
%---EXAMPLE USAGE:
% dataCell = {randn(1000,1),randn(500,1)*2+1};
% BF_JitteredParallelScatter(dataCell,true,true,true);

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
% ------------------------------------------------------------------------------

if nargin < 2
    addMeans = true;
    % Add strip for mean of each group by default
end

if nargin < 3
    doveTail = true;
    % Add kernel distribution
end

if nargin < 4
    makeFigure = true;
end

if nargin < 5
    extraParams = struct;
end

% ------------------------------------------------------------------------------
% Set extra parameters

numGroups = length(dataCell);

% Custom marker for points within the distribution:
if ~isfield(extraParams,'customSpot')
    customSpot = '.';
else
    customSpot = extraParams.customSpot;
end

% Custom x-offset
if ~isfield(extraParams,'customOffset')
    customOffset = 0;
else
    customOffset = extraParams.customOffset;
end

% Custom offset range; width of each jitter:
% (proportion of each consecutive unit interval to use for the plot)
if ~isfield(extraParams,'offsetRange')
    offsetRange = 0.5;
else
    offsetRange = extraParams.offsetRange;
end

% Custom colormap
if ~isfield(extraParams,'theColors')
    if numGroups <= 3
        theColors = BF_getcmap('set1',numGroups,1);
    else
        theColors = BF_getcmap('spectral',numGroups,1);
    end
    if length(theColors) < numGroups
        theColors = arrayfun(@(x)zeros(3,1),1:numGroups,'UniformOutput',0);
    end
else
    theColors = extraParams.theColors;
    % fprintf(1,'Using custom colors\n');
end

% add zero line
if ~isfield(extraParams,'add0Line')
    add0Line = false;
else
    add0Line = extraParams.add0Line;
end

% ------------------------------------------------------------------------------

% Reset random number generator for reproducibility:
rng('default');

if makeFigure
    figure('color','w');
end
hold on; box('on');

lwidth = 2;

% ------------------------------------------------------------------------------
% Add kernel distribution
% ------------------------------------------------------------------------------
if doveTail
    ff = cell(numGroups,1);
    xx = cell(numGroups,1);
    for i = 1:numGroups
        kk=ceil(i/2);
        if isempty(dataCell{i})
            continue
        end
        if any(isnan(dataCell{i}))
            warning('NaNs in dataCell')
        end
        [f, x] = ksdensity(dataCell{i},'npoints',500);
        f = f/max(f);
        % Only keep range where data exists:
        minKeep = max([1,find(x>=min(dataCell{i}),1,'first')]);
        maxKeep = min([length(x),find(x>=max(dataCell{i}),1,'first')]);
        keepR = minKeep:maxKeep;
        % keepR = (x>=min(dataCell{i}) & x<=max(dataCell{i}));
        x = x(keepR);
        f = f(keepR);
        plot(customOffset+i+f*offsetRange/2,x,'-','color',theColors{kk},'LineWidth',lwidth)
        plot(customOffset+i-f*offsetRange/2,x,'-','color',theColors{kk},'LineWidth',lwidth)

        % Plot top and bottom
        plot(customOffset+i+[-f(1),+f(1)]*offsetRange/2,min(x)*ones(2,1),'-','color',theColors{kk},'LineWidth',0.1)
        plot(customOffset+i+[-f(end),f(end)]*offsetRange/2,max(x)*ones(2,1),'-','color',theColors{kk},'LineWidth',0.1)

        % Keep for dovetailing the jittered scatter points:
        xx{i} = x; ff{i} = f;
    end
end

% ------------------------------------------------------------------------------
% Plot jittered scatter for each group
% ------------------------------------------------------------------------------
if ~isempty(customSpot)
    for i = 1:numGroups
        kk=ceil(i/2);
        if doveTail
            xRand = zeros(length(dataCell{i}),1);
            for j = 1:length(dataCell{i})
                try
                    xRand(j) = (rand(1)*offsetRange-offsetRange/2)*ff{i}(find(xx{i} >= dataCell{i}(j),1));
                end
            end
            %i + rand([length(dataCell{i}),1])*offsetRange-offsetRange/2;
        else
            xRand = rand([length(dataCell{i}),1])*offsetRange-offsetRange/2;
        end
        plot(customOffset + i + xRand,dataCell{i},customSpot,'color',theColors{kk})
    end
end

% ------------------------------------------------------------------------------
% Add strips for means and stds:
% ------------------------------------------------------------------------------
brightenAmount = -0.3;
if any(cellfun(@(x)any(isnan(x)),dataCell)), warning('NaNs in data'); end
for i = 1:numGroups
    kk=ceil(i/2);
    try
        brightColor = brighten(theColors{kk},brightenAmount);
    catch
        brightColor = theColors{kk};
    end

    % plot([customOffset + i - offsetRange/2,customOffset + i + offsetRange/2],nanmean(dataCell{i})*ones(2,1),'-',...
                            % 'color',brightColor,'LineWidth',lwidth)
    % plot([customOffset + i - offsetRange/2,customOffset + i + offsetRange/2],(nanmean(dataCell{i})-nanstd(dataCell{i}))*ones(2,1),'--',...
                            % 'color',brightColor,'LineWidth',lwidth)
    % plot([customOffset + i - offsetRange/2,customOffset + i + offsetRange/2],(nanmean(dataCell{i})+nanstd(dataCell{i}))*ones(2,1),'--',...
                            % 'color',brightColor,'LineWidth',lwidth)
    % plot([customOffset + i - offsetRange/2,customOffset + i + offsetRange/2],nanmedian(abs(dataCell{i}))*ones(2,1),'-',...
                            % 'color',brightColor,'LineWidth',lwidth)
    plot([customOffset + i - offsetRange/2,customOffset + i + offsetRange/2],nanmedian(dataCell{i})*ones(2,1),'-',...
                            'color',brightColor,'LineWidth',lwidth)

end

% ------------------------------------------------------------------------------
% Add dotted line at 0
% ------------------------------------------------------------------------------
if add0Line
    plot([0:numGroups+1],zeros(1,numGroups+2),':','Color','k')
end


end
