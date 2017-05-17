% BayesianSTC.m

fileName = sprintf('5Min_PinkNoiseMovie%d.mat',1);
load(fileName);
load('5Min_UnbiasedPinkNoiseMovie1.mat');
xaxis = linspace(-maxPix/2,maxPix/2,DIM(2));
yaxis = linspace(-minPix/4,3*minPix/4,DIM(1));

%S = normrnd(0,1,[N,N,numStimuli]);

totalCentisecs = 5*60*100;

stimTimes = round((0:1/60:5*60-1/60).*100);
pointProcessStimTimes = zeros(totalCentisecs,1);
for ii=1:numStimuli-1
    pointProcessStimTimes(stimTimes(ii)+1:stimTimes(ii+1)) = ii;
end
pointProcessStimTimes(stimTimes(numStimuli):end) = numStimuli;

% DEFINE SPATIOTEMPORAL RECEPTIVE FIELD (STRF) AND PLOT EXAMPLE CONVOLUTION
%   OF STRF WITH SEQUENCE STIMULUS
%  the temporal component is a zero-gain, time-causal filter, which
%  endows the "neuron" with its preferred temporal frequency, its transient
%  responses to visual images, and its latency
%   if the v (velocity) is zero, then this is a space-time separable
%   filter, but if v ~= 0, it's space-time inseparable
gaborFun = @(x,y,t,k,n,v,A,xc,yc,sigmax,sigmay,spatFreq,theta,phi) ...
    exp(-((x-xc).*cos(A)-(y-yc).*sin(A)).^2./(2*sigmax*sigmax)-...
    ((x-xc).*sin(A)+(y-yc).*cos(A)).^2/(2*sigmay*sigmay))...
    .*sin((2*pi.*spatFreq).*(cos(theta-pi/2).*(x-xc)+sin(theta-pi/2).*(y-yc)+v.*t)-phi)...
    .*(k.*t).^n.*exp(-k.*t).*(1/gamma(n+1)-(k.*t).^2./(gamma(n+3)));

t = 200:-1:0;
sigmax = 100; % size (arbitrary units) of the spatial Gabor
sigmay = 80;
k = 0.3; % k and n define the temporal filter
n = 20;
A = pi/4;
spatFreq = 0.1;
theta = 60*pi/180; % neuron's preferred orientation
phi = 0; % phase of the sinusoid
v = 0; % velocity (shift in the phase over time)


[X,Y,T] = meshgrid(xaxis,yaxis,t);
filter1 = gaborFun(X,Y,T,k,n,v,A,0,0,sigmax,sigmay,spatFreq,theta,phi);

%filter1 = filter1./max(abs(filter1(:)));
%filter1 = filter1./sum(filter1(:));

filter1 = filter1(1:end-30,35:end-35,1:10:end);
filterLen = size(filter1,3);
onScreenStims = zeros(totalCentisecs,filterLen);
convResult = zeros(totalCentisecs,1);
y = zeros(totalCentisecs,1); % neural response

sigmoid = @(x,A,B,C) A.*exp((x-B).*C)./(1+exp((x-B).*C));
historyB = [-2,-0.3,0,0.1,0.4,0.25,0.15,0.1,0.075,0.05]';
baseRate = 5/100;

historyLen = length(historyB);
y(1:filterLen) = poissrnd(baseRate,[filterLen,1]);
currentHistory = fliplr(y(filterLen-historyLen+1:filterLen)');
A = 3;B = 12;C = 0.5;

startTime = filterLen+1;
for ii=startTime:totalCentisecs
%     if mod(ii,10) == 0
%         onScreenStims(ii,:) = pointProcessStimTimes(ii-filterLen+1:ii)';
%         temp = double(S(:,:,onScreenStims(ii,:)));
% 
%         temp1 = temp.*filter1;
%         convResult(ii) = sum(temp1(:));
%     else
%        onScreenStims(ii,:) = onScreenStims(ii-1,:);
%        convResult(ii) = convResult(ii-1); 
%     end
    onScreenStims(ii,:) = pointProcessStimTimes(ii-filterLen+1:ii)';
    temp = double(S(1:end-30,35:end-35,onScreenStims(ii,:)));
      
    temp1 = temp.*filter1;
    convResult(ii) = sum(temp1(:));
    lambda = baseRate*exp(currentHistory*historyB)*exp(sigmoid(convResult(ii),A,B,C));
    y(ii) = poissrnd(lambda);
    currentHistory = [y(ii),currentHistory(1:end-1)];
end
%y = logical(y);
historyParams = historyLen;

historyDesign = zeros(totalCentisecs,historyParams);

for kk=1:historyParams
    temp = y;shift = zeros(kk,1);
    history = [shift;temp];
    historyDesign(:,kk) = history(1:(end-kk));
end

y = y(filterLen+1:end);
historyDesign = historyDesign(filterLen+1:end,:);
onScreenStims = onScreenStims(filterLen+1:end,:);

clear history pointProcessStimTimes stimTimes X Y T;

STA = 0;

spikeTimes = find(y>0);
numSpikes = length(spikeTimes);
for ii=1:numSpikes
   tempInds = onScreenStims(spikeTimes(ii),:);
   movie = unbiasedS(1:end-30,35:end-35,tempInds);
   STA = STA+y(ii).*movie(:);
end
STA = STA./sum(y);

clearvars -except numSpikes STA spikeTimes unbiasedS onScreenStims y;
STC = zeros(length(STA),length(STA),'single');

for ii=1:numSpikes
    tempInds = onScreenStims(spikeTimes(ii),1:2:end);
    movie = unbiasedS(1:end-30,35:end-35,tempInds);
    STC = STC+(movie(:)-STA)*(movie(:)-STA)';
end
STC = STC./numSpikes;
[V,D] = eig(STC);

figure();plot(diag(D));
figure();histogram(STA);

% for a quadratic LNP model of spiking activity
%   f(x) = exp(0.5*x'*C*x+b'*x+a);
% then the ML estimates are
%  C = inv(phi)-inv(STC);
%  b = inv(STC)*STA
%  a = log(numSpikes/totalSeconds*sqrt(det(phi*inv(STC)))-0.5*STA'*inv(phi)*inv(STC)*STA;

%   the model provides a direct estimate of the spiking non-linearity