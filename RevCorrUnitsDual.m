function [] = RevCorrUnitsDual(AnimalName,Date,NoiseType)
%RevCorrUnitsDual.m
%   Same as RevCorrUnits, but will try to rethink how the linear receptive
%    field is calculated by having two receptive fields, which represent
%    light and dark being presented
%   %   Analysis of single unit recording data in response to a series of white
%   or pink noise stimuli (see Noise_RevCorr.m for Psychtoolbox
%   stimulus)
%    See Smyth et al. 2003 Receptive Field Organization ... We'll solve a
%    matrix equation using the regularized pseudo-inverse technique.
%
% Briefly, if the receptive field for a given neuron is described by a
%  p-by-1 vector (p being the number of pixels), f, and the stimuli described
%  by a s-by-p (s being the number of stimuli presented) matrix, S, and the 
%  response of that neuron described by an s-by-1 vector, r, then we can write
%  the responses in terms of the stimuli and the receptive field as 
%  r = S * f .  We want to infer the spatial receptive field structure,
%  f, so we need S^-1 * r = S^-1 * S * f => f = S^-1 * r ... So, the
%  following code will get the spatial receptive fields f for a series of 
%  neurons whose responses were measured with single-unit recordings during
%  the display of the white/pink noise stimuli (in S).
%
%INPUT: AnimalName - unique identifier for the animal as a number, e.g.
%            12345
%       Date - date of the experiment, e.g. 20160525
%    
%      Input through saved file 'NoiseStimData_AnimalName.mat'
%      S - number of stimuli-by-number of pixels matrix of white noise
%        stimuli used by the function Noise_RevCorr.m
%      sorted spike data from Plexon, as a cell array, called allts, time
%      stamp data, also from Plexon, as a cell array, called tsevs
%OUTPUT: a saved file with the name similar to the above and a suffix, -Results
%       F - number of neurons-by-number of pixels matrix representing the
%         receptive field for each neuron recorded from
%            The receptive field for a given neuron can be 
%             visualized by:
%            figure();imagesc(reshape(F(1,:),[DIM(1),DIM(2)]);      
%    
% Created: 2016/03/04, 24 Cummington, Boston
%  Byron Price
% Updated: 2017/07/17
% By: Byron Price

% assume we sample at some frequency, Fs, and have a point process
%  representation of the spiking behavior of a given neuron, we 
%  have the set of stimuli used in the input matrix S and the 
%  timeStamps for when each was displayed
%  If neuron i fired after the onset of stimulus j, then we'll take 
%  the number of spikes in a certain window after the onset of that
%  stimulus as the response
%
%   Note: we could record from these neurons for a period with no visual
%   stimulus present and get a background firing rate for each and then
%   normalize the response by that background firing rate, though this may
%   not add anything


% read in the .plx file
beta = 0;

EphysFileName = strcat('NoiseData',NoiseType,num2str(Date),'_',num2str(AnimalName),'-sort.mat');

if exist(EphysFileName,'file') ~= 2
    readall(strcat(EphysFileName(1:end-4),'.plx'));pause(1);
end

StimulusFileName = strcat('NoiseStim',NoiseType,num2str(Date),'_',num2str(AnimalName),'.mat');
load(EphysFileName,'nunits1','allts','adfreq','allad','svStrobed','tsevs')
load(StimulusFileName)

% gaborFun = @(x,y,B,A,xc,yc,sigmax,sigmay,spatFreq,theta,phi) ...
%     B.*exp(-((x-xc).*cos(A)-(y-yc).*sin(A)).^2./(2*sigmax*sigmax)-...
%     ((x-xc).*sin(A)+(y-yc).*cos(A)).^2/(2*sigmay*sigmay))...
%     .*sin((2*pi.*spatFreq).*(cos(theta-pi/2).*(x-xc)-sin(theta-pi/2).*(y-yc))-phi);
% 
% %nonLinFun = @(x,baseX,scale) exp((x-baseX)/scale);
% nonLinFun = @(x,base,slope,rise) rise./(1+exp(-(x-base).*slope));
%nonLinFun = @(x,slope,intercept) max(0,slope.*x+intercept);

% horzDegrees = atan((screenPix_to_effPix*DIM(1)*conv_factor/10)/DistToScreen);
% vertDegrees = atan((screenPix_to_effPix*DIM(2)*conv_factor/10)/DistToScreen);

DIM = zeros(2,1);
if size(effectivePixels,2) == 1
    DIM(1) = sqrt(effectivePixels);
    DIM(2) = sqrt(effectivePixels);
elseif size(effectivePixels,2) == 2
    DIM(1) = effectivePixels(2);
    DIM(2) = effectivePixels(1);
end
    
xaxis = linspace(-round(screenPix_to_effPix*DIM(2)/2)+1,...
    round(screenPix_to_effPix*DIM(2)/2),DIM(2));
yaxis = linspace(round(3*screenPix_to_effPix*DIM(1)/4),...
    -round(screenPix_to_effPix*DIM(1)/4)+1,DIM(1));

[X,Y] = meshgrid(xaxis,yaxis);

% CREATE UNBIASED VERSION OF MOVIE BY DIVIDING OUT POWER SPECTRUM
%  USED TO GENERATE THE MOVIE

u = [(0:floor(DIM(1)/2)) -(ceil(DIM(1)/2)-1:-1:1)]'/DIM(1);
v = [(0:floor(DIM(2)/2)) -(ceil(DIM(2)/2)-1:-1:1)]'/DIM(2);
[V,U] = meshgrid(v,u);
S_f = (U.^2+V.^2).^(beta/2);

S_f(S_f==inf) = 0;
S_f = 1./S_f;
S_f(S_f==inf) = 0;
S_f = sqrt(S_f);
a = 0;b = 255;

[~,numPixels] = size(S);
unbiasedS = zeros(numStimuli,numPixels*2);
revisedS = zeros(numStimuli,numPixels*2);
for ii=1:numStimuli
    temp = reshape(real(ifftn(fftn(double(reshape(S(ii,:),[DIM(1),DIM(2)]))).*S_f)),[DIM(1)*DIM(2),1])';
    currentMin = min(temp);currentMax = max(temp);
    temp = ((b-a).*(temp-currentMin))/(currentMax-currentMin)+a;
    
    light = max(temp,127)-127; % how much brightness
    dark = max(255-temp,127)-127; % how much darkness
    unbiasedS(ii,1:numPixels) = light;
    unbiasedS(ii,numPixels+1:end) = dark;
    
    light = max(S(ii,:),127)-127;
    dark = max(255-S(ii,:),127)-127;
    revisedS(ii,1:numPixels) = light;
    revisedS(ii,numPixels+1:end) = dark;
end

if strcmp(NoiseType,'pinkHC') == 1
   unbiasedS(unbiasedS<60) = 0;unbiasedS(unbiasedS>=60 & unbiasedS<196) = 127;
   unbiasedS(unbiasedS>=196) = 255;
end

clear U V u v S;

% REORGANIZE SPIKING DATA
temp = ~cellfun(@isempty,allts);
Chans = find(sum(temp,1));numChans = length(Chans);
totalUnits = sum(sum(temp))-numChans;

temp = cell(totalUnits,1);
count = 1;
for ii=1:numChans
   for jj=2:nunits1
       if isempty(allts{jj,Chans(ii)}) == 0
           temp{count} = allts{jj,Chans(ii)};
           count = count+1;
       end
   end
end
allts = temp;

strobeStart = 33;
strobeData = tsevs{1,strobeStart};

if length(strobeData) ~= numStimuli && length(unique(svStrobed)) == 1
    strobeData = strobeData(1:2:end);
elseif length(unique(svStrobed)) == 3
    strobeData = strobeData(svStrobed==1);
end

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
timeMultiplier = 1000;
totalMillisecs = round(totalTime*timeMultiplier);

stimTimes = round(strobeData.*timeMultiplier);

if exist('flipIntervals','var')==1
    flipInterval = mean(flipIntervals);
end

pointProcessSpikes = zeros(totalMillisecs,totalUnits);

for ii=1:totalUnits
   spikeTimes = max(1,round(allts{ii}.*timeMultiplier));
   for jj=1:length(spikeTimes)
      pointProcessSpikes(spikeTimes(jj),ii) = 1;
   end
end

reducedMov = zeros(numStimuli,1*timeMultiplier);
reducedSpikeCount = zeros(totalUnits,numStimuli,1*timeMultiplier);
for ii=1:totalUnits
    for jj=1:numStimuli
        timeInds = stimTimes(jj):stimTimes(jj)+1*timeMultiplier-1;
        reducedSpikeCount(ii,jj,:) = pointProcessSpikes(timeInds,ii)';
        reducedMov(jj,:) = movement(timeInds);
    end
end

clearvars -except EphysFileName totalUnits numStimuli ...
    reducedSpikeCount DIM unbiasedS allts strobeData xaxis yaxis ...
    X Y totalMillisecs reducedMov movement svStrobed revisedS S_f;

% GLM with Gaussian basis functions
% fullSize = DIM(1)*DIM(2);
% basisStdDevs = [250,500,750,1000,1250,1500,1750,2000,2500];
% numStdDevs = length(basisStdDevs);
% 
% finalResultsPoissonB = cell(totalUnits,numStdDevs);
% finalResultsPoissonDev = cell(totalUnits,numStdDevs);
% finalResultsPoissonSE = cell(totalUnits,numStdDevs);
% temp = randperm(numStimuli);
% divide = round(0.7*numStimuli);
% train = temp(1:divide);test = temp(divide+1:end);clear temp divide;
% for zz=1:totalUnits
%    spikeTrain = squeeze(reducedSpikeCount(zz,:,:));
%    spikeTrain = sum(spikeTrain(:,50:300),2);
%    movDesign = sum(reducedMov(:,50:300),2);
%    
%    %    r = spikeTrain;
%    %    fhat = unbiasedS\r;
%    gaussFun = @(x,y,xc,yc,std) exp(-((x-xc).*(x-xc))./(2*std*std)-...
%        ((y-yc).*(y-yc))./(2*std*std));
%    
%    center1 = xaxis(1:3:end);
%    center2 = yaxis(1:3:end);
%    numBasis1 = length(center1);
%    numBasis2 = length(center2);
%    totalParams = numBasis1*numBasis2;
%    
%    basisFuns = zeros(fullSize,totalParams);
%    for stddev = 1:numStdDevs
%        count = 1;
%        for ii=1:numBasis1
%            for jj=1:numBasis2
%                temp = gaussFun(X,Y,center1(ii),center2(jj),basisStdDevs(stddev));
%                basisFuns(:,count) = temp(:)./max(temp(:));
%                count = count+1;
%            end
%        end
%        design = [ones(numStimuli,1),unbiasedS*basisFuns];
%        [bPoiss,~,statsPoiss] = glmfit(design(train,:),spikeTrain(train),'poisson','constant','off');
%        finalResultsPoissonB{zz,stddev} = bPoiss;
%        finalResultsPoissonSE{zz,stddev} = statsPoiss.se;
%        
%        mu = exp(design(test,:)*bPoiss);
%        deviance = spikeTrain(test).*log(spikeTrain(test)./mu)-(spikeTrain(test)-mu);
%        deviance(isnan(deviance)) = mu(isnan(deviance));
%        devPoiss = 2.*sum(deviance);
%        finalResultsPoissonDev{zz,stddev} = devPoiss;
%        clear bPoiss devPoiss statsPoiss count design;
%    end
%    clear spikeTrain movDesign basisFuns center1 center2 numBasis1 numBasis2 totalParams;
% end
% 
% fileName = strcat(EphysFileName(1:end-9),'-GLMResults.mat');
% save(fileName,'finalResultsPoissonB',...
%     'finalResultsPoissonDev',...
%     'finalResultsPoissonSE','allts','totalUnits',...
%     'reducedSpikeCount','DIM','unbiasedS','movement',...
%     'xaxis','yaxis','basisStdDevs','reducedMov','strobeData','svStrobed','totalMillisecs');

% REGULARIZED PSEUDO-INVERSE SOLUTION, CREATE CONVOLUTION MATRIX L
S = double(S);
fullSize = DIM(1)*DIM(2);
L = sparse(fullSize,fullSize);

%operator = [0,-1,0;-1,4,-1;0,-1,0];
bigCount = 1;
for jj=1:DIM(2)
    for ii=1:DIM(1)
        tempMat = zeros(DIM(1),DIM(2));
        tempMat(ii,jj) = 4;
        if ii > 1
            tempMat(ii-1,jj) = -1;
        end
        if ii < DIM(1)
            tempMat(ii+1,jj) = -1;
        end
        if jj > 1
            tempMat(ii,jj-1) = -1;
        end
        if jj < DIM(2)
            tempMat(ii,jj+1) = -1;
        end
        L(bigCount,:) = tempMat(:)';bigCount = bigCount+1;
    end
end

numLambda = 20;
loglambda = logspace(3,7,numLambda);
F = zeros(totalUnits,fullSize);
bestLambda = zeros(totalUnits,1);
heldOutDeviance = zeros(totalUnits,5);
heldOutExplainedVariance = zeros(totalUnits,1);
sigmoidNonlin = zeros(totalUnits,4);
for zz=1:totalUnits
   fprintf('Running unit %d ...\n',zz);
   spikeTrain = squeeze(reducedSpikeCount(zz,:,:));
   
   % divide data into training and test
   temp = randperm(numStimuli);
   divide = round(0.7*numStimuli);
   train = temp(1:divide);test = temp(divide+1:end);clear temp divide;
   
   baseRate = sum(sum(spikeTrain(:,51:400)))./(numStimuli*0.35);
   spikeTrain = sum(spikeTrain(:,50:300),2);
   
   r = [spikeTrain(train);sparse(fullSize,1)];
   
   tempF = zeros(numLambda,fullSize);
   tempSigmoid = zeros(numLambda,4);
   tempDev = zeros(numLambda,4);
   tempExplainVar = zeros(numLambda,1);
%    allOnesTrain = ones(length(train),1);
%    allOnesTest = ones(length(test),1);
   for jj=1:numLambda
       % calculate regularized pseudoinverse solution
      constraints = [unbiasedS(train,:);loglambda(jj).*L];
      fhat = constraints\r;fhat = full(fhat);
      
      % remove bias created by using pink noise
      temp = reshape(fhat,[DIM(1),DIM(2)]);
      temp = real(ifft2(fft2(full(temp)).*S_f));
      fhat = temp(:);
      
      % fit sigmoid nonlinearity
      result = S(train,:)*fhat;
      myFun = @(x) x(1)./(1+exp(-(result-x(2)).*x(3)))+x(4)-spikeTrain(train);
      x0 = [5,0.5,10,baseRate];
      lb = [0,0,0,0];ub = [5e2,Inf,Inf,Inf];
%       options.Algorithm = 'levenberg-marquardt';
      sigmoidParams = lsqnonlin(myFun,x0,lb,ub);
      
      % save parameters
      tempF(jj,:) = fhat;
      tempSigmoid(jj,:) = sigmoidParams;
      
      % calculate held-out Poisson deviance
      temp = S(test,:)*fhat;
      temp = sigmoidParams(1)./(1+exp(-(temp-sigmoidParams(2)).*sigmoidParams(3)))+sigmoidParams(4);
      
      initialDev = spikeTrain(test).*log(spikeTrain(test)./temp)-(spikeTrain(test)-temp);
      initialDev(isnan(initialDev) | isinf(initialDev)) = temp(isnan(initialDev) | isinf(initialDev));
      tempDev(jj,1) = 2*sum(initialDev);
      
      % calculate held-out explained variance
      rr = corrcoef(temp,spikeTrain(test));
      tempExplainVar(jj,1) = rr(2,1)^2;
      
      % rotate the receptive field and recalculate the held-out deviance
      rf = reshape(fhat,[DIM(1),DIM(2)]);
      for kk=2:4
        rf = imrotate(rf,90);
        tempfhat = rf(:);
        temp = S(test,:)*tempfhat;
        temp = sigmoidParams(1)./(1+exp(-(temp-sigmoidParams(2)).*sigmoidParams(3)))+sigmoidParams(4);
        
        initialDev = spikeTrain(test).*log(spikeTrain(test)./temp)-(spikeTrain(test)-temp);
        initialDev(isnan(initialDev) | isinf(initialDev)) = temp(isnan(initialDev) | isinf(initialDev));
        tempDev(jj,kk) = 2*sum(initialDev);
      end
   end
  
   [~,bestMap] = min(tempDev(:,1));
   F(zz,:) = tempF(bestMap,:);
   bestLambda(zz) = loglambda(bestMap);
   heldOutDeviance(zz,1:4) = tempDev(bestMap,:);
   sigmoidNonlin(zz,:) = tempSigmoid(bestMap,:);
   heldOutExplainedVariance(zz,1) = tempExplainVar(bestMap);
   
   f = ones(length(train),1)\spikeTrain(train);
   heldOutDeviance(zz,5) = sum((spikeTrain(test)-allOnesTest*f).^2);
end

% save the results
fileName = strcat(EphysFileName(1:end-9),'-PseudoInvResults.mat');
save(fileName,'F','totalUnits','bestLambda',...
    'reducedSpikeCount','DIM','unbiasedS','movement',...
    'xaxis','yaxis','reducedMov','allts','strobeData','totalMillisecs',...
    'svStrobed','heldOutDeviance','numStimuli','S_f','sigmoidNonlin',...
    'heldOutExplainedVariance');

% REVERSE CORRELATION SOLUTION
% for ii=1:numChans
%     numNeurons = Neurons(ii,1);
%     figure();plotRows = ceil(numNeurons/2);
%     for jj=1:numNeurons
%         r = squeeze(Response(ii,jj,:));
% % %         r = newS*f
%         fhat = newS\r;
%         F(ii,jj,:) = fhat;
%         subplot(plotRows,2,jj);imagesc(reshape(fhat,[N,N]));
%         title(sprintf('RevCorr: Chan %d, Unit %d, Animal %d',ii,jj,AnimalName));
%     end
% end
% Img = reshape(S(tt,:),[minPix/screenPix_to_effPix,minPix/screenPix_to_effPix]);
% Img = kron(double(Img),ones(screenPix_to_effPix));

end