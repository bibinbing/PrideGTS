function [lambda_hat,H,x,source,mn] = init_vbica(data,source_type,states,isonoise,priors,net_init,source_init,eta)
% [lambda_hat,H,x,source,mn] = init_vbica(data,source_type,states,isonoise,priors,net_init,source_init,eta)
%
% Initialisation of vbICA model parameters. Called from 'initialise_ica'.
%
%
% -----------
% Input
% -----------
%
% Necessary paramteres
%
% data        No. rows = No. sensors, No. columns = No. data points
% source_type Mixture of 'g'aussians     (Default)
%             Mixture of 'e'xponentials
%             Mixture of 't'runcated gaussians
%             'hg','he','ht' - HMM versions
% states      No. states per source model eg [2 4 3 5] for 4 sources
% isonoise    1: learn isotropic noise
%             0: learn non-isotropic diagonal noise
% priors      Data structure of priors
% net_init    'SVD', 'rand'om
% source_init 'kmeans', 'co_mean', 'rand'om
%
%
% Optional parameters
%
% eta         Class responsibilities if being trained as part of a
%             Mixture of ICAs, otherwise = 1
%
%
%
% -----------
% Output
% -----------
%
% lambda_hat  Initialised expectation of noise precision
% H           Initialised expectation of mixmatrix
% x           Initialised expectation of hidden source signals
% source      Cell array of initialised source models
% mn          Initialised network mean vector
%
%
% ------------------------------------------------------
%
% Original code by Rizwan Choudrey
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)


global CHECK_PROGRESS;


if nargin<8
    eta = 1;
end


sources = length(states);
[sensors,points] = size(data);

SQUARE = 0;
if sensors == sources
    SQUARE = 1;
end

if sources>sensors
    net_init = 'rand';
    fprintf('More sources than sensors - using %s...\n',net_init);
end

if length(source_type)==2
    src_type = source_type(2);
else
    src_type = source_type;
end


MIXTURES=0;
if length(eta)==1
    switch src_type
        case 'g'
            mn = mean(data')';
            init_data = data - repmat(mn,1,points);
        otherwise
            mn = zeros(sensors,1);
            init_data = data;
    end
else
    MIXTURES = 1;
    switch src_type
        case 'g'
            big_eta = repmat(eta,sensors,1);
            mn = sum(big_eta.*data,2)./sum(eta);
            data = data - repmat(mn,1,points);
        otherwise
            mn = zeros(sensors,1);
    end
    %[a,index] = max(eta);
    index = find(eta==1);
    init_data = data(:,index);
end
[sensors,points] = size(init_data);

fprintf('Initialising observation network parameters using %s...\n',net_init);

switch src_type

    case 'g'



        switch net_init

            case 'SVD'
                %set mixing matrix, x and noise using SVD solution


                [U,lambda,V] =svd(init_data',0);
                clear U

                %set mixing matrix to appropriate vectors in SVD
                m_H_hat=V(:,1:sources);
                H = m_H_hat;
                clear V;


                if SQUARE
                    % For square mixing, it is vital to have a very low initial
                    % noise precision otherwise we'll never infer that
                    % there is added noise
                    var_noise = var(init_data')./5;
                    isovalue = mean(1./var_noise);
                else
                    %set noise to residual singular values
                    values = diag(lambda);
                    explained = values(1:sources);
                    unexplained = values(sources+1:sensors);
                    clear lambda;

                    %Need to hack average noise precision to allow 'softer' convergence
                    mean_crap = (mean(unexplained))^2/points;
                    isovalue = 1/(mean_crap);
                    clear explained unexplained mean_crap;
                end


                noise = ones(sensors,1)*isovalue;
                noise_prec_vect = noise;
                lambda_hat = noise_prec_vect;
                clear isovalue noise;


                %set x to reconstruct from SVD
                x=inv(H'*H)*H'*init_data;
                mu_x_hat=x;
            case 'rand'
                H = randn(sensors,sources);
                x=inv(H'*H)*H'*init_data;
                m_H_hat = H;
                mu_x_hat=x;
                var_noise = var(init_data')./5;
                lambda_hat = mean(1./var_noise)*ones(sensors,1);


        end


        %initialise sources
        fprintf('Initialising source parameters using %s..',source_init);
        for i=1:sources
            fprintf('.%i',i);
            if CHECK_PROGRESS
                hDebugFig = findobj('Type', 'figure', 'Name', 'VBICA PDFS');
                if isempty(hDebugFig)
                    hDebugFig = figure('Name', 'VBICA PDFS', 'NumberTitle', 'off');
                end
                set(0, 'CurrentFigure', hDebugFig(1));
                subplot(1,sources,i)
            end
            source{i} = initialise_mix1d(x(i,:),states(i),source_type,source_init,priors.source);
            if MIXTURES
                gams = repmat(source{i}.pi',1,length(eta));
                gams(:,index) = source{i}.gammas;
                source{i}.gammas = gams;
            end

        end

        if MIXTURES
            allX = repmat(mean(x')',1,length(eta));
            allX(:,index)=x;
            x=allX;
            mu_x_hat = x;
        end



    otherwise
        m_H_hat = randn(sensors,sources);
        H = rect_expect(m_H_hat,1000*ones(size(m_H_hat)));
        mu_x_hat=randn(sources,points);
        %x = randn(sources,points);
        x = rect_expect(mu_x_hat,1000*ones(size(mu_x_hat)));
        %x = init_data(1:sources,:)/10;
        x = mu_x_hat.^2;
        var_noise = var(init_data')./100;
        var_noise = 0.00001;
        lambda_hat = mean(1./var_noise)*ones(sensors,1);

        %initialise Positive Distbns
        fprintf('Initialising positive-source parameters using %s and %s..',source_type,source_init);

        if MIXTURES
            allX = repmat(mean(x')',1,length(eta));
            allX(:,index)=x;
            x=allX;
            mu_x_hat = x;
        end


        for i=1:sources
            fprintf('.%i',i);
            if CHECK_PROGRESS
                hDebugFig = findobj('Type', 'figure', 'Name', 'VBICA PDFS');
                if isempty(hDebugFig)
                    hDebugFig = figure('Name', 'VBICA PDFS', 'NumberTitle', 'off');
                end
                set(0, 'CurrentFigure', hDebugFig(1));
                subplot(1,sources,i)
            end
            source{i} = initialise_mix1d(x(i,:),states(i),source_type,source_init,priors.source);
        end

end
fprintf('\n\n\n');
