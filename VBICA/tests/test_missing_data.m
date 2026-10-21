function test_missing_data
% Test the unchanged RunvbICA interface and masked vbICA2 training only.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(root));
oldRng = rng;
cleanup = onCleanup(@() rng(oldRng));
rng(41);
N = 240;
s = [sin((1:N)/9); 2*(rand(1,N)>.65)-.7];
data = [1 .4; -.5 1.1; .8 -.7]*s + [.2;-.4;.7] + .03*randn(3,N);
data(rand(size(data))<.25) = NaN;
data(2,70:100) = NaN;
data(:,120:122) = NaN;

% Independent one-source posterior calculation.
y = [2 NaN 4; NaN 3 5];
H = [1;2]; H2 = [1.2;4.3]; lambda = [2;3];
src = struct('type','g','m',1,'gammas',ones(1,3), ...
    'precs',4,'centres',.5);
[~,~,~,~,mu,beta] = recon_source2(y,zeros(1,3),H,H2,lambda,{src},1e-8,1,0);
assert(max(abs(beta{1}-[6.4 16.9 19.3]))<1e-10);
assert(max(abs(mu{1}-[6 20 40]./[6.4 16.9 19.3]))<1e-10);
fprintf('PASS: source posterior uses observed entries only.\n');

for iso=0:1
    for ard=0:1
        rng(17);
        evalc('net = vbica2(data,[2 2],''g'',0,1e-4,8,iso,ard);');
        assert(all(isfinite([net.recon(:);net.mixmatrix(:);net.noise(:);net.energy_path(:)])));
        if iso
            expectedCount = nnz(~isnan(data));
        else
            expectedCount = sum(~isnan(data),2);
        end
        assert(all(abs(net.posts.noise.c-(.001+expectedCount/2))<1e-10));
    end
end
fprintf('PASS: random gaps, continuous gaps, missing epochs, noise and ARD settings.\n');

% Exactly the five original output fields; no completion output or new option.
rng(17);
evalc('result = RunvbICA(data'',1,2,''g'',0,1e-4,8,1);');
expectedFields = {'Contribution';'FilteredData';'OriginTimeSeries'; ...
    'SpatialResponsesMatrix';'Components'};
assert(isequal(fieldnames(result),expectedFields));
assert(isequaln(result.OriginTimeSeries,data'));
assert(isequal(isnan(result.FilteredData),isnan(data')));
assert(isequal(size(result.Components),[N 3]));
assert(all(isfinite(result.Components(:))));
assert(all(isfinite(result.SpatialResponsesMatrix(:))));
assert(all(isfinite(result.Contribution)));
fprintf('PASS: original arguments and output fields; finite components from NaN input.\n');

% No fully observed epoch: complete-case row deletion cannot pass this case.
rng(52);
incomplete = randn(4,160);
for t=1:size(incomplete,2)
    incomplete(mod(t-1,4)+1,t) = NaN;
end
incomplete(2,40:80) = NaN;
incomplete(1,100) = 0; % An observed zero is still a valid observation.
assert(all(any(isnan(incomplete),1)));
evalc('result = RunvbICA(incomplete'',1,2,''g'',0,1e-4,8,1);');
assert(all(isfinite(result.Components(:))));
assert(all(isfinite(result.SpatialResponsesMatrix(:))));
assert(isequaln(result.OriginTimeSeries,incomplete'));
assert(isequal(fieldnames(result),expectedFields));
assert(isequal(size(result.Components),[160 4]));
fprintf('PASS: decomposition without any complete epoch; observed zeros retained.\n');
end
