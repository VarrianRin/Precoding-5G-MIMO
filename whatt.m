function [CALC, GRAPH, LOAD_R, GRAPH_ALL, BUILD, CAPACITY, LLOYD] = whatt(todo)

% 1 - load saved correlations matrix and build precoder and plot
% 4 - load saved correlations matrix and build precoder
% 6 - no calculation, plot
% 7 - no calculation, plot every saved cdf with V
% 8 - generate correlation matrix and build precoder and plot

LLOYD = 1;

switch todo

    case 1
    CALC        = 1;
    GRAPH       = 1;
    LOAD_R      = 1;
    GRAPH_ALL   = 0;
    BUILD       = 1;
    CAPACITY    = 1;

    case 4
    CALC        = 1;
    GRAPH       = 0;
    LOAD_R      = 1;
    GRAPH_ALL   = 0;
    BUILD       = 1;
    CAPACITY    = 1;

    case 6
    CALC        = 0;
    GRAPH       = 1;
    LOAD_R      = 0;
    GRAPH_ALL   = 0;
    BUILD       = 0;
    CAPACITY    = 0;

    case 7
    CALC        = 0;
    GRAPH       = 1;
    LOAD_R      = 0;
    GRAPH_ALL   = 1;
    BUILD       = 0;
    CAPACITY    = 0;

    case 8
    CALC        = 1;
    GRAPH       = 1;
    LOAD_R      = 0;
    GRAPH_ALL   = 0;
    BUILD       = 1;
    CAPACITY    = 1;    

end
