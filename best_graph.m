%% %% %% --- --- --- --- --- --- BEST GRAPH --- --- --- --- --- --- %% %% %%
%
%                            building 2D plots
% 
% 
% 
%
% type               = 1 'sub' subplot graphs
%                      2 'one' graphs on one figure
%                      3 'dif' graphs on different figures
% 
% 
% type_format        = 1
%                      2 type correction
%                      3
% 
% graph_type         = 1 'usl' usual plot
%                      2 'cdf' cdf of data
%                      3 'log' log
% 
% line_format        = Color
%                      LineStyle
%                      LineWidth
%                      Marker
%                      MarkerIndices
%                      MarkerEdgeColor
%                      MarkerFaceColor
%                      MarkerSize
% 
% % %                      ' ' or [] - use previous(!) or default if no previous
% % %                      n = none
% % % 
% % %                      COLOR (default black) (two symbols)
% % %                      white              w1
% % %                      light grey         w2
% % %                      dark grey          w3  W
% % %                      black              w4
% % % 
% % %                      light red          r1
% % %                      red2               r2
% % %                      red3               r3  R
% % %                      dark red           r4
% % %                      
% % %                      light blue         b1
% % %                      blue2              b2
% % %                      blue3              b3  B
% % %                      dark blue          b4                     
% % % 
% % %                      light green        g1
% % %                      green2             g2
% % %                      green3             g3  G
% % %                      dark green         g4    
% % % 
% % %                      light yellow       y1
% % %                      yellow2            y2
% % %                      yellow3            y3  Y
% % %                      dark yellow        y4
% % % 
% % %                      light purple       p1
% % %                      purple2            p2
% % %                      purple3            p3  P
% % %                      dark purple        p4   
% % % 
% % %                      
% % %                      LINESTYLE (default -) (one symbol)
% % %                      -      -
% % %                      --     =
% % %                      -.     ,
% % %                      :      :
% % %                      none   n
% % % 
% % % 
% % %                      LINEWIDTH (default 2) (one symbol)
% % %                      positive value 1 = 1/72 inch                    !!!max 9!!!
% % % 
% % % 
% % %                      MARKER (default none) (one symbol)
% % %                      o + * . x _ | s d ^ v > < p h n [none]
% % % 
% % % 
% % %                      MARKERINDICES (default []) (one symbol)
% % %                      positive int values of how many markers to show !!!max 9!!!
% % % 
% % % 
% % %                      MARKERENGECOLOR (default red2) (two symbols)
% % %                      same as colors
% % % 
% % % 
% % %                      MARKERFACECOLOR (defult light red) (two symbols)
% % %                      same as colors
% % % 
% % %                      MARKERSIZE (default 6) (one symbol)
% % %                      positive values 1 = 1/72 inch                   !!!max 9!!!  
% 
% 
% TO-DO
% different size of graphs
% figures of subplots
% establish default colors and line_types
% in one cell many graphs
% graph_type - unconvinient

function [] = best_graph(type, x, y, graph_type, type_format, line_format, name, legen, xlab, ylab)

%% ---------------------------------- TECHNICAL ----------------------------------

% -- default line formats -- 
def_cl = '#000000';
def_ls = '-';
def_lw = 2;
def_mc = 'none';
def_mi = 'default';
def_me = 'auto';
def_mf = 'auto';
def_ms = 6;
NO_PAR = 8;
MX_PAR = 11;

FormatBank = cell(128, 5);

FormatBank{double('W')} = '#696969';
FormatBank{double('w'), 1 + 1} = '#ffffff';
FormatBank{double('w'), 1 + 2} = '#bdbdbd';
FormatBank{double('w'), 1 + 3} = FormatBank(double('W'));
FormatBank{double('w'), 1 + 4} = '#000000';
FormatBank{double('G')} = '#00ff88';
FormatBank{double('g'), 1 + 1} = '#a3ff75';
FormatBank{double('g'), 1 + 2} = '#6fb378';
FormatBank{double('g'), 1 + 3} = FormatBank(double('G'));
FormatBank{double('g'), 1 + 4} = '#007d3a';
FormatBank{double('R')} = '#ff0000';
FormatBank{double('r'), 1 + 1} = '#ff7575';
FormatBank{double('r'), 1 + 2} = '#915656';
FormatBank{double('r'), 1 + 3} = FormatBank(double('R'));
FormatBank{double('r'), 1 + 4} = '#7d0000';
FormatBank{double('B')} = '#0011ff';
FormatBank{double('b'), 1 + 1} = '#75a3ff';
FormatBank{double('b'), 1 + 2} = '#566a91';
FormatBank{double('b'), 1 + 3} = FormatBank(double('B'));
FormatBank{double('b'), 1 + 4} = '#002a7d';
FormatBank{double('Y')} = '#ffff00';
FormatBank{double('y'), 1 + 1} = '#ffff75';
FormatBank{double('y'), 1 + 2} = '#bdb060';
FormatBank{double('y'), 1 + 3} = FormatBank(double('Y'));
FormatBank{double('y'), 1 + 4} = '#877500';
FormatBank{double('O')} = '#ff8800';
FormatBank{double('o'), 1 + 1} = '#ffbf75';
FormatBank{double('o'), 1 + 2} = '#b39c6f';
FormatBank{double('o'), 1 + 3} = FormatBank(double('O'));
FormatBank{double('o'), 1 + 4} = '#7d4b00';
FormatBank{double('P')} = '#d400ff';
FormatBank{double('p'), 1 + 1} = '#e875ff';
FormatBank{double('p'), 1 + 2} = '#aa72b5';
FormatBank{double('p'), 1 + 3} = FormatBank(double('P'));
FormatBank{double('p'), 1 + 4} = '#8c00a8';
FormatBank{double('-')} = '-';
FormatBank{double('=')} = '--';
FormatBank{double(',')} = '-.';
FormatBank{double(':')} = ':';
FormatBank{double('o')} = 'none';
FormatBank{double('+')} = 'o';
FormatBank{double('*')} = '+';
FormatBank{double('.')} = '*';
FormatBank{double('x')} = '.';
FormatBank{double('_')} = 'x';
FormatBank{double('|')} = '_';
FormatBank{double('s')} = '|';
FormatBank{double('d')} = 's';
FormatBank{double('^')} = 'd';
FormatBank{double('v')} = '^';
FormatBank{double('>')} = 'v';
FormatBank{double('<')} = '>';
FormatBank{double('p')} = '<';
FormatBank{double('h')} = 'p';
FormatBank{double('n')} = 'h';
% ----------------------------------------------------------------------

%% --------------------------------- DEFAULTS --------------------------------
if nargin < 2
    error('??? You must specify type and x !!!');
end

if ~exist('x', 'var')               || isempty(x)
    x               = {[]};
end

if ~exist('y', 'var')               || isempty(y)
    y               = {[]};
end

if ~exist('graph_type', 'var')      || isempty(graph_type)
    graph_type      = {1};
end

if ~exist('type_format', 'var')     || isempty(type_format)
    type_format     = {1};
end

if ~exist('line_format', 'var')     || isempty(line_format)
    line_format     = {'R-'};
end

if ~exist('name', 'var')            || isempty(name)
    name            = {'Title'};
end

if ~exist('legen', 'var')           || isempty(legen)
    legen           = {'data1'};
end

if ~exist('xlab', 'var')            || isempty(xlab)
    xlab            = {'x axis'};
end

if ~exist('ylab', 'var')            || isempty(ylab)
    ylab            = {'y axis'};
end

%% --------------------------- START -------------------------------------
if isa(x, 'cell') || isa(y, 'cell')

    if ~isempty(x)
        N_incell = length(x);
    elseif ~isempty(y)
        N_incell = length(y);
    end

    for icell = 1 : N_incell

        % -- extract data --
        supp        = x{icell};
        supp_type   = graph_type{icell};
        values      = y{icell};

        if isempty(supp) && ~isempty(values)
            if size(values, 1) >= size(values, 2)
                values = values.';
            end
            supp = repelem(1, size(values, 2), size(values ,1), 1);
        end
        if size(supp, 1) >= size(supp, 2)
            supp = supp.';
        end

        PlotNumber = size(supp, 1);
        PlotLength = size(supp, 2);

        if supp_type == 2
            supp = sort(supp, 2);
        elseif supp_type == 3
            error('!!! NOT COMPLETED !!!');
        end

        if isempty(values) && supp_type == 2
            for n = 1 : PlotNumber
                values(n, :) = 0 : 1/PlotLength : 1-1/PlotLength;
            end
        elseif isempty(values) && supp_type == 1
            for n = 1 : PlotNumber
                values(n, :) = 1 : PlotLength;
            end
        end

        % -- extract line format -- 

        LineFormat = cell(PlotNumber, NO_PAR);
        prevs      = {def_cl, def_ls, def_lw, def_mc, def_mi, def_me, def_mf, def_ms};
        InitLength = length(line_format);

        for iPlot = 1 : PlotNumber

            if InitLength >= iPlot && isa(line_format{iPlot}, 'char')

                format  = line_format{iPlot};
                last    = length(format);
                par     = 1;
                i       = 0;

                for j = 1 : MX_PAR % -- max size of format --
                    
                    i = i + 1;
                    if i > last
                        for other = par : NO_PAR % -- all others specifications are from previous -- 
                            LineFormat{iPlot, other} = prevs{other};
                        end
                        break;

                    elseif format(i) == ' ' % -- prev or default --

                        LineFormat{iPlot, par} = prevs{par};
                        par = par + 1;
                        continue;

                    elseif par == 1 || par == 6 || par == 7 % -- colors --

                        val = double(format(i));
                        if val < double('Z') && val > double('A')
                            LineFormat{iPlot, par} = FormatBank{val};
                        elseif val < double('z') && val > double('a')
                            LineFormat{iPlot, par} = FormatBank{val, str2double(format(1 + i)) + 1};
                            i = i + 1;
                        end
                        par = par + 1;
                        continue;

                    elseif par == 3 || par == 8 % -- number -- 

                        LineFormat{iPlot, par} = str2double(format(i));
                        par = par + 1;
                        continue;

                    elseif par == 2 || par == 4 % -- symbol --

                        LineFormat{iPlot, par} = FormatBank{double(format(i))};
                        par = par + 1;
                        continue;

                    elseif par == 5 % -- markers amount --

                        val = str2double(format(i)) + 1;
                        dif = PlotLength / val;
                        LineFormat{iPlot, par} = int32(dif : dif : PlotLength - dif);
                        par = par + 1;
                        continue;

                    end
                end
            
            else
                for other = 1 : NO_PAR % -- all of specifications are form previous -- 
                    LineFormat{iPlot, other} = prevs{other};
                end
            end

            prevs = {LineFormat{iPlot, 1}, LineFormat{iPlot, 2}, LineFormat{iPlot, 3}, LineFormat{iPlot, 4},...
                     LineFormat{iPlot, 5}, LineFormat{iPlot, 6}, LineFormat{iPlot, 7}, LineFormat{iPlot, 8}};
        end


        
        %% -------------------------------- building graph --------------------------        
          
        if type == 1
            error('!!! NOT COMPLETED !!!');
        elseif type == 2

            figure;
            plot(supp(1, :), values(1, :), 'Color', LineFormat{1, 1}, 'LineStyle', LineFormat{1, 2}, ...
                'LineWidth', LineFormat{1, 3}, 'Marker', LineFormat{1, 4}, 'MarkerIndices', LineFormat{1, 5}, ...
                'MarkerEdgeColor', LineFormat{1, 6}, 'MarkerFaceColor', LineFormat{1, 7}, 'MarkerSize', LineFormat{1, 8});
            hold on;

            for i = 2 : PlotNumber
                plot(supp(i, :), values(i, :), 'Color', LineFormat{i, 1}, 'LineStyle', LineFormat{i, 2}, ...
                'LineWidth', LineFormat{i, 3}, 'Marker', LineFormat{i, 4}, 'MarkerIndices', LineFormat{i, 5}, ...
                'MarkerEdgeColor', LineFormat{i, 6}, 'MarkerFaceColor', LineFormat{i, 7}, 'MarkerSize', LineFormat{i, 8});
            end

            legend(legen, 'Location', 'southoutside', 'FontSize', 12, 'NumColumns', 2);
            set(gca, 'FontSize', 14);
            set(gca, 'Box', 'on');
            set(gcf, 'Color', [1, 1, 1]);
            grid on;

            title(name{icell}, 'FontSize', 16);
            xlabel(xlab{icell}, 'FontSize', 16);
            ylabel(ylab{icell}, 'FontSize', 16);
            hold off;

        elseif type == 3
            error(' !!! NOT COMPLETED !!!');
        end
    end
end

end



                   
































