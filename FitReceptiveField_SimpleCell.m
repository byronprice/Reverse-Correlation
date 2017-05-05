function [PosteriorSamples,PosteriorMean,PosteriorInterval,Likelihood] =...
    FitReceptiveField_SimpleCell(AnimalName,Date,NoiseType)
%FitReceptiveField_SimpleCell.m
%  Use data from a receptive-field mapping experiment and fit a Gabor model
%   to obtain the spatiotemporal receptive field for a simple cell in L4 of V1.
%  MCMC

% declare global variables
global numStimuli totalMillisecs pointProcessStimTimes h ...
    historyDesign unbiasedS numParameters historyParams ...
    movement stimulusTime X Y T gaborFun nonLinFun spikeTrain ...
    trainMillisecs;
% read in the .plx file
beta = 0;

EphysFileName = strcat('NoiseMovieData',NoiseType,num2str(Date),'_',num2str(AnimalName),'-sort.mat');

if exist(EphysFileName,'file') ~= 2
    readall(strcat(EphysFileName(1:end-4),'.plx'));pause(1);
end

StimulusFileName = strcat('NoiseMovieStim',NoiseType,num2str(Date),'_',num2str(AnimalName),'.mat');
load(EphysFileName,'nunits1','allts','adfreq','allad','svStrobed','tsevs')
load(StimulusFileName)

gaborFun = @(x,y,t,k,n,v,A,xc,yc,sigmax,sigmay,spatFreq,theta,phi) ...
    exp(-((x-xc).*cos(A)-(y-yc).*sin(A)).^2./(2*sigmax*sigmax)-...
    ((x-xc).*sin(A)+(y-yc).*cos(A)).^2/(2*sigmay*sigmay))...
    .*sin((2*pi.*spatFreq).*(cos(theta-pi/2).*(x-xc)+sin(theta-pi/2).*(y-yc)+v.*t)-phi)...
    .*(k.*t).^n.*exp(-k.*t).*(1/gamma(n+1)-(k.*t).^2./(gamma(n+3)));

%nonLinFun = @(x,baseX,scale) exp((x-baseX)/scale);
nonLinFun = @(x,base,slope,rise,drop) rise.*exp((x-base)./slope)./(1+exp((x-base)./slope))-drop;
%nonLinFun = @(x,slope,intercept) max(0,slope.*x+intercept);

% horzDegrees = atan((screenPix_to_effPix*DIM(1)*conv_factor/10)/DistToScreen);
% vertDegrees = atan((screenPix_to_effPix*DIM(2)*conv_factor/10)/DistToScreen);
xaxis = linspace(-screenPix_to_effPix*DIM(2)/2,...
    screenPix_to_effPix*DIM(2)/2,DIM(2));
yaxis = linspace(-screenPix_to_effPix*DIM(1)/4,...
    3*screenPix_to_effPix*DIM(1)/4,DIM(1));
taxis = linspace(200,0,21);

stimulusTime = 10;

[X,Y,T] = meshgrid(xaxis,yaxis,taxis);

% CREATE UNBIASED VERSION OF MOVIE BY DIVIDING OUT POWER SPECTRUM
%  USED TO GENERATE THE MOVIE
DIM = [effectivePixels(1),effectivePixels(2),numStimuli];

u = [(0:floor(DIM(1)/2)) -(ceil(DIM(1)/2)-1:-1:1)]'/DIM(1);
v = [(0:floor(DIM(2)/2)) -(ceil(DIM(2)/2)-1:-1:1)]'/DIM(2);
tau = [(0:floor(DIM(3)/2)) -(ceil(DIM(3)/2)-1:-1:1)]'/(DIM(3));
[V,U,TAU] = meshgrid(v,u,tau);
S_f = (U.^spaceExp+V.^spaceExp+TAU.^timeExp).^(beta/2);

S_f(S_f==inf) = 1;

unbiasedS = real(ifftn(fftn(double(S))./S_f));
bigTemp = zeros(DIM(1),DIM(2),DIM(3)+1);
temp = ones(DIM(1),DIM(2)).*mean(unbiasedS(:));
bigTemp(:,:,1) = temp;
bigTemp(:,:,2:end) = unbiasedS;

unbiasedS = bigTemp;

clear S S_f U V TAU u v tau temp bigTemp;

% REORGANIZE SPIKING DATA
temp = ~cellfun(@isempty,allts);
Chans = find(sum(temp,1));numChans = length(Chans);
totalUnits = sum(sum(temp));

temp = cell(totalUnits,1);
count = 1;
for ii=1:numChans
   for jj=1:nunits1
       if isempty(allts{jj,Chans(ii)}) == 0
           temp{count} = allts{jj,Chans(ii)};
           count = count+1;
       end
   end
end
allts = temp;

strobeStart = 33;
strobeData = tsevs{1,strobeStart};


% GATHER LFP AND MOVEMENT DATA
nonEmptyAD = ~cellfun(@isempty,allad);
inds = find(nonEmptyAD==1);
LFP = cell(length(inds),1);
for ii=1:length(inds)
   LFP{ii} = allad{inds};
end

totalTime = length(LFP{1})./adfreq;

if isempty(allad{49}) == 0
    movement = allad{49};

    difference = length(LFP{1})-length(movement);
    
    if mod(difference,2) == 0
        addOn = difference/2;
        movement = [zeros(addOn-1,1);movement;zeros(addOn+1,1)];
    else
        addOn = floor(difference/2);
        movement = [zeros(addOn,1);movement;zeros(addOn+1,1)];
    end
    tempMov = conv(abs(movement),ones(adfreq/2,1),'same');
    tempMov = tempMov-mean(tempMov);
    stdEst = 1.4826*mad(tempMov,1);
    movement = single(tempMov>(3*stdEst));
    clear tempMov stdEst;
else
    movement = zeros(length(LFP{1}),1);
end

% COLLECT DATA IN THE PRESENCE OF VISUAL STIMULI
totalMillisecs = round(totalTime*1000);

trainMillisecs = round(totalMillisecs*0.7);
% timeVec = 202:stimulusTime:trainMillisecs;
% timeVecLen = length(timeVec);


stimTimes = round(strobeData.*1000);
pointProcessStimTimes = ones(totalMillisecs,1);

for kk=1:numStimuli-1
   pointProcessStimTimes(stimTimes(kk):stimTimes(kk+1)-1) = kk+1;
end
pointProcessStimTimes(stimTimes(numStimuli):(stimTimes(numStimuli)+1000/movie_FrameRate)) = numStimuli+1;

pointProcessSpikes = zeros(totalMillisecs,totalUnits);

for ii=1:totalUnits
   spikeTimes = round(allts{ii}.*1000);
   for jj=1:length(spikeTimes)
      pointProcessSpikes(spikeTimes(jj),ii) = 1;
   end
end


historyParams = 25;
gaborParams = 11;
nonLinParams = 4;
numParameters = historyParams+gaborParams+nonLinParams+1;

twopi = 2*pi;
Bounds = zeros(numParameters,2);
Bounds(1:historyParams,1) = -200;
Bounds(1:historyParams,2) = 200;
Bounds(historyParams+1,:) = [1e-2,50];% k time-filter parameter
Bounds(historyParams+2,:) = [1e-1,100]; % n time-filter parameter
Bounds(historyParams+3,:) = [-100,100]; % v time-filter parameter
Bounds(historyParams+4,:) = [0,twopi]; % A orientation of Gabor
Bounds(historyParams+5,:) = [min(xaxis)-5,max(xaxis)+5]; % x center
Bounds(historyParams+6,:) = [min(yaxis)-5,max(yaxis)+5]; % y center
Bounds(historyParams+7,:) = [1,500]; % standard deviation x
Bounds(historyParams+8,:) = [1,500]; % standard deviation y
Bounds(historyParams+9,:) = [0,100]; % spatial frequency
Bounds(historyParams+10,:) = [0,twopi]; %  orientation theta
Bounds(historyParams+11,:) = [0,twopi]; % phase shift phi
Bounds(historyParams+12,:) = [-200,200]; % sigmoid base
Bounds(historyParams+13,:) = [1e-4,200]; % sigmoid slope
Bounds(historyParams+14,:) = [0,200]; % sigmoid rise
Bounds(historyParams+15,:) = [-200,200]; % sigmoid drop
Bounds(end,:) = [-200,200]; % parameter for movement modulation

numIter = 5e5;burnIn = 1e5;skipRate = 50;

PosteriorSamples = zeros(totalUnits,numParameters,length(burnIn:skipRate:numIter));
PosteriorMean = zeros(totalUnits,numParameters);
PosteriorInterval = zeros(totalUnits,numParameters,2);
Likelihood = zeros(totalUnits,length(burnIn:skipRate:numIter));
for zz=1:totalUnits
        % ORGANIZE DATA
        
        spikeTrain = pointProcessSpikes(:,zz);
        historyDesign = zeros(totalMillisecs,historyParams);
            
        % get history dependence (past 25ms)
        historyDesign(:,1) = ones(totalMillisecs,1,'single');
        for kk=2:historyParams
            temp = spikeTrain;shift = zeros(kk-1,1);
            history = [shift;temp];
            historyDesign(:,kk) = history(1:(end-(kk-1)));
        end
        
        [b,dev,stats] = glmfit([historyDesign,movement],spikeTrain,'poisson','constant','off');
        
        
        % RUN MCMC
        h = ones(numParameters,1)./1000;

        numIter = 5e5;burnIn = 1e5;
        parameterVec = zeros(numParameters,numIter);
        
        parameterVec(1:historyParams,1) = b(1:historyParams);
        parameterVec(end,1) = b(end);
        
        parameterVec(historyParams+1,1) = normrnd(0.5,0.1);
        parameterVec(historyParams+2,1) = normrnd(15,5);
        parameterVec(historyParams+3,1) = normrnd(0,0.1);
        parameterVec(historyParams+4,1) = normrnd(pi,pi/4);
        parameterVec(historyParams+5,1) = normrnd(0,300);
        parameterVec(historyParams+6,1) = normrnd(0,200);
        parameterVec(historyParams+7,1) = normrnd(150,50);
        parameterVec(historyParams+8,1) = normrnd(150,50);
        parameterVec(historyParams+9,1) = normrnd(0.4,0.1);
        parameterVec(historyParams+10,1) = normrnd(pi,pi/4);
        parameterVec(historyParams+11,1) = normrnd(pi,pi/4);
        parameterVec(historyParams+12,1) = normrnd(2,1);
        parameterVec(historyParams+13,1) = normrnd(1,0.5);
        parameterVec(historyParams+14,1) = normrnd(1,0.5);
        parameterVec(historyParams+15,1) = normrnd(0,1);
        
        parameterVec(:,1) = max(Bounds(:,1),min(parameterVec(:,1),Bounds(:,2)));
        
        likelihoodXprev = GetLikelihood(parameterVec(:,1));
        tempLikelihood = zeros(numIter,1);
        tempLikelihood(1) = likelihoodXprev;
        
        scatter(1,tempLikelihood(1));hold on;pause(0.1);
        
        sigma = zeros(numParameters,1);
        sigma(1:historyParams) = abs(0.5.*b(1:historyParams));
        sigma(historyParams+1) = 0.2;
        sigma(historyParams+2) = 5;
        sigma(historyParams+3) = 0.5;
        sigma(historyParams+4) = twopi/4;
        sigma(historyParams+5) = 50;
        sigma(historyParams+6) = 50;
        sigma(historyParams+7) = 20;
        sigma(historyParams+8) = 20;
        sigma(historyParams+9) = 0.5;
        sigma(historyParams+10) = twopi/4;
        sigma(historyParams+11) = twopi/4;
        sigma(historyParams+12) = 1;
        sigma(historyParams+13) = 1;
        sigma(historyParams+14) = 1;
        sigma(historyParams+15) = 1;
        sigma(end) = abs(0.5*b(end));
        
        rejectRate = 0;
        for iter = 2:numIter
            xStar = parameterVec(:,iter-1)+normrnd(0,sigma,[numParameters,1]);
            xStar([historyParams+4,historyParams+10,historyParams+11]) = ...
                mod(xStar([historyParams+4,historyParams+10,historyParams+11]),twopi);
            temp = xStar-Bounds(:,1);
            
            if sum(temp<0) > 0 % outside bounds
                parameterVec(:,iter) = parameterVec(:,iter-1);
                rejectRate = rejectRate+1;
            else
                likelihoodXstar = GetLikelihood(xStar);
                scatter(iter,likelihoodXstar);hold on;pause(0.1);
                %         A = min(1,(likelihoodXstar*priorXstar)./(likelihoodXprev*priorXprev));
                %             A = min(1,likelihoodXprev/likelihoodXstar);
                if log(rand) < (likelihoodXstar-likelihoodXprev)
                    parameterVec(:,iter) = xStar;
                    likelihoodXprev = likelihoodXstar;
                    display('improve');
                else
                    parameterVec(:,iter) = parameterVec(:,iter-1);
                    rejectRate = rejectRate+1;
                end
            end
            tempLikelihood(iter) = likelihoodXprev;
        end
        
        fprintf('Rejection Rate: %3.2f',rejectRate/numIter);
        
        PosteriorSamples(zz,:,:) = parameterVec(:,burnIn:skipRate:numIter);
        Likelihood(zz,:) = tempLikelihood(:,burnIn:skipRate:numIter);
        
        alpha = 0.05;
        
        numColumns = 4;
        numRows = floor(numParameters/numColumns)+mod(numParameters,numColumns);
        for ii=1:numParameters
           subplot(numRows,numColumns,ii);
           histogram(squeeeze(PosteriorSamples(zz,ii,:)));
           title('Marginal Posterior, Parameter %d',ii);
           PosteriorMean(zz,ii) = mean(PosteriorSamples(ii,:));
           PosteriorInterval(zz,ii,:) = quantile(squeeze(PosteriorSamples(zz,ii,:)),...
               [alpha/2,1-alpha/2]);
        end
end


end

function [loglikelihood] = GetLikelihood(parameterVec)
global unbiasedS pointProcessStimTimes historyParams historyDesign ...
    stimulusTime X Y T gaborFun nonLinFun spikeTrain ...
    movement trainMillisecs;

baseMu = historyDesign(1:trainMillisecs,:)*parameterVec(1:historyParams);
moveMu = parameterVec(end).*movement(1:trainMillisecs);

gabor = gaborFun(X,Y,T,parameterVec(historyParams+1),parameterVec(historyParams+2),...
    parameterVec(historyParams+3),parameterVec(historyParams+4),parameterVec(historyParams+5),...
    parameterVec(historyParams+6),parameterVec(historyParams+7),...
    parameterVec(historyParams+8),parameterVec(historyParams+9),...
    parameterVec(historyParams+10),parameterVec(historyParams+11));
%filterEnergy = sum(gabor(:).*gabor(:));

filterOutput = zeros(trainMillisecs,1);
%stimLen = 29000:stimulusTime:trainMillisecs;

for kk=30000:stimulusTime:trainMillisecs
%     onScreenStim = squeeze(newS(pointProcessStimTimes(zz(kk)),:,:));
%     if mod(kk,stimulusTime) == 1
%         filterOutput = sum(newS(pointProcessStimTimes(zz(kk),:))'.*gaborFilter)./numPixels;
%     end
    filterOutput(kk-stimulusTime/2:kk+stimulusTime/2-1) = ...
        sum(sum(sum(unbiasedS(:,:,pointProcessStimTimes(kk-201:stimulusTime:kk))...
        .*gabor)));
end

filterOutput(1:30000) = sum(sum(sum(unbiasedS(:,:,pointProcessStimTimes(1:stimulusTime:201))...
        .*gabor)));

% for kk=1:numStimuli
%    filterOutput = sum(newS(kk,:)'.*gaborFilter)./numPixels;
%    loglikelihood(pointProcessStimTimes==kk) = Y(pointProcessStimTimes==kk).*log(baseMu*exp(filterOutput))-...
%        baseMu*exp(filterOutput);
% end

%figure();plot(filterOutput);

filterOutput = nonLinFun(filterOutput,parameterVec(historyParams+12),...
    parameterVec(historyParams+13),parameterVec(historyParams+14),...
    parameterVec(historyParams+15));


loglikelihood = sum(spikeTrain(1:trainMillisecs).*(baseMu+moveMu+filterOutput)-...
    exp(baseMu+moveMu+filterOutput));
end

