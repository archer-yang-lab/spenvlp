%% choose_spenv
% Select the dimension of the envelope subspace using Bayesian information
% criterion.

%% Syntax
%         [u,results] = choose_spenv(X, Y, range)
%         [u,results] = choose_spenv(X, Y, range, Opts)
%
%% Input
%
% *X*: Predictors. An n by p matrix, p is the number of predictors and n 
% is the number of observations. The predictors can be univariate or 
% multivariate, discrete or continuous.
% 
% *Y*: Multivariate responses. An n by r matrix, r is the number of
% responses. The responses must be continuous variables.
% 
% *Opts*: A list containing the optional input parameters, to control the
% iterations in sg_min. If one or several (even all) fields are not
% defined, the default settings are used.
% 
% * Opts.maxIter: Maximum number of iterations.  Default value: 300.
% * Opts.ftol: Tolerance parameter for F.  Default value: 1e-10. 
% * Opts.gradtol: Tolerance parameter for dF.  Default value: 1e-7.
% * Opts.verbose: Flag for print out dimension selection process, 
% logical 0 or 1. Default value: 0.
%
%% Output
%
% *u*: Dimension of the envelope. An integer between 0 and r.

%% Description
% This function implements the Bayesian information criteria (BIC) to select
% the dimension of the envelope subspace.  

%% Example
%         load wheatprotein.txt
%         X = wheatprotein(:, 8);
%         Y = wheatprotein(:, 1:6);
%         u = choose_spenv(X, Y,0.1,0.05)

function [u_lambda,results] = choose_spenv(X, Y, range, Opts)
[n r] = size(Y);
if nargin < 2
    error('Inputs: X, Y should be specified!');
elseif nargin == 2
	range=[0,r-1];
    Opts = [];
elseif nargin == 3
    Opts = [];
end
if range(2)>r-1
    range(2)=r-1;
end

Opts = make_opts(Opts);
printFlag = Opts.verbose;
Opts.verbose = 0;
lambda_choose=[1e-2,5e-2,0.1,0.2,0.5];

dims=range(1):range(2);
bic_seq=zeros(length(dims)+1,1);
aic_seq=zeros(length(dims)+1,1);
l_seq=zeros(length(dims)+1,1);
pvalue_seq=ones(length(dims)+1,1);
df_seq=zeros(length(dims)+1,1);
lambda_seq=zeros(length(dims)+1,1);
%if range(1)==0 
%    c=length(dims)-1;
%else
%    c=length(dims);
%end
%lambda_seq=zeros(c,1);


ModelOutput0 = spenv(X, Y, r, lambda_choose(1), Opts);
bic_seq(length(dims)+1) = - 2 * ModelOutput0.l + log(n) * ModelOutput0.paramNum;
aic_seq(length(dims)+1)= -2 * ModelOutput0.l + 2 * ModelOutput0.paramNum;
l_seq(length(dims)+1) = - 2 * ModelOutput0.l ;

%%
for i = 1 : length(dims) 
    i
    %if printFlag == 1
%         fprintf(['================Current dimension ' int2str(dims(i)) '\n']);
    %end
    % Chooose lambda
    if dims(i)>0 && length(lambda_choose)>1
%         fprintf(['======Current lambda ' num2str(lambda_choose(1)) '\n']);
		m0=env(X,Y,dims(i));
        Opts.init=m0.Gamma;
        ModelOutput=spenv(X,Y,dims(i),lambda_choose(1), Opts);
        lambda=lambda_choose(1);
        for k=2:length(lambda_choose)
            ModelOutput1 = spenv(X,Y,dims(i),lambda_choose(k), Opts);
%             fprintf(['=====Current lambda ' num2str(lambda_choose(k)) '\n']);
            if ModelOutput1.BIC<ModelOutput.BIC
                lambda=lambda_choose(k);
                ModelOutput=ModelOutput1;
            end
        end
        lambda_seq(i)=lambda;
    else
        ModelOutput = spenv(X, Y, dims(i),lambda_choose(1), Opts);
    end
    bic_seq(i) = - 2 * ModelOutput.l + log(n) * ModelOutput.paramNum;
	aic_seq(i) = -2 * ModelOutput.l + 2 * ModelOutput.paramNum;
	l_seq(i) = - 2 * ModelOutput.l ;
	chisq = l_seq(i) - (- 2 * ModelOutput0.l) ;
	df = ModelOutput0.paramNum - ModelOutput.paramNum;
    df_seq(i)=df;
	pvalue_seq(i)=1-chi2cdf(chisq, df); 
end


if length(dims) == 1
    
    u = [dims dims dims dims];
    lambda = [lambda_seq(1) lambda_seq(1) lambda_seq(1) lambda_seq(1)]
    u_lambda=[u;lambda];
    results=[dims bic_seq(1) aic_seq(1) l_seq(1) pvalue_seq(1) df_seq(1) lambda_seq(1)];
    
else
%%
    dims2=[dims r];
    [t1,index1]=min(bic_seq);
    u1=dims2(index1);
    [t2,index2]=min(aic_seq);
    u2=dims2(index2);
    index3=find(pvalue_seq>0.05);
    index3=index3(1);
    u3=dims2(index3); 
    index4=find(pvalue_seq>0.01);
    index4=index4(1);
    u4=dims2(index4); 
    u=[u1 u2 u3 u4];
    lambda=[lambda_seq(index1) lambda_seq(index2) lambda_seq(index3) lambda_seq(index4)];
    u_lambda=[u;lambda];
    results=[dims2' bic_seq aic_seq l_seq pvalue_seq df_seq lambda_seq];

end

