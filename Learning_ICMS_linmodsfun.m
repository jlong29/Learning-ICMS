function [dif, expect] = Learning_ICMS_linmodsfun(params, SelectParams, p, Info4FitFun)
% This code specifies a series of linear models applied to the figures of
% the 'Learning ICMS' paper. 
%%%%%%%%%%%%%%
%%% Inputs %%%
%%%%%%%%%%%%%%
% params:       the values of the current iterations free parameters
% SelectParams: indices of the free parameters
% p:            the complete vector of parameters for a model
% Info4FitFun:  a structure with fields relevant to current model

% Info4FitFun.fields:
    % .data   = data for fitting the model
    % .SEtype = sets whether SE are known or not
    % .SEdat  = the standard errors for the data
    % .nIter  = indexes the model for this run
    
%%%%%%%%%%%%%%%
%%% Outputs %%%
%%%%%%%%%%%%%%%
% dif:         scaled residual between model and data
% expect:      the bestfit model values

% Check data size
data             = Info4FitFun.data;
[Nrows,Ncolumns] = size(Info4FitFun.data);

% Subsequent models are specified in terms of a single parameter vector p.
% This iterations floating parameters are set as per 'SelectParams' while
% all fixed parameters are left at their preset values.
p(SelectParams)=params;  %This iterations selected params are taken from params and put in p

% Determine model type
modelInd = Info4FitFun.nIter;

if modelInd == 1
    % a simple single pararmeter for a DC shift    
    expect = p(1)*ones(1,Ncolumns);
elseif modelInd == 2
    % two parameters: second and third columns have same parameter
    expect = p(1)*[0 1 1] + p(2)*[1 0 0];
elseif modelInd == 3
    % two parameters: first and third columns have same parameter
    expect = p(1)*[1 0 1] + p(2)*[0 1 0];
elseif modelInd == 4
    % two parameters: first and second columns have same parameter
    expect = p(1)*[1 1 0] + p(2)*[0 0 1];
end

% Compare model to data
if Info4FitFun.SEtype == 1
    SE   = Info4FitFun.SEdat;
else
    SE   = 1;
end

% Calculate Residuals betwene model and data
dif  = (expect-data)./SE;

