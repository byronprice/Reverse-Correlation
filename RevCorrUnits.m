function [] = RevCorrUnits(AnimalName,Date,NoiseType)
%RevCorrUnits.m
%   %   Analysis of single unit recording data in response to a series of white
%   or pink noise stimuli (see Noise_RevCorr.m for Psychtoolbox
%   stimulus)
%    See Smyth et al. 2003 Receptive Field Organization ... We'll solve a
%    matrix equation using a regularized pseudo-inverse technique.
%
% Briefly, if the receptive field for a given neuron is described by a
%  p-by-1 vector (p being the number of pixels), f, and the stimuli described
%  by a s-by-p (s being the number of stimuli presented) matrix, S, and the 
%  response of that neuron described by an s-by-1 vector, r, then we can write
%  the responses in terms of the stimuli and the receptive field as 
%  r = S * f .  We want to infer the spatial receptive field structure
%  f, so we need S^-1 * r = S^-1 * S * f => f = S^-1 * r ... So, the
%  following code will get the spatial receptive fields f for a series of 
%  neurons whose responses were measured with single-unit recordings during
%  the display of the white/pink noise stimuli in S.
%
%INPUT: AnimalName - unique identifier for the animal as a number, e.g.
%            12345
%       Date - date of the experiment, e.g. 20160525
%       Chans - channel numbers, input as [6,8], defaults to 6 and 8
%    
%      Input through saved file 'NoiseStimData_AnimalName.mat'
%      S - number of stimuli-by-number of pixels matrix of white noise
%        stimuli used by the function Noise_RevCorr.m
%       R - matrix of size number of frames-by-number of neurons that
%        contains neuronal data which has already been processed into a
%        point process and which will be converted to a new matrix of responses, 
%        number of stimuli-by-number of neurons matrix of neuronal
%        numbers obtained from calcium imaging (count threshold
%        crossings as spikes ... we will want to know the number of spikes
%        in response to the onset and offset of each stimuli.

%OUTPUT: F - number of pixels-by-number of neurons matrix representing the
%         receptive field for each neuron recorded from in R
%            if N = sqrt(size(F,1));
%            then the receptive field for a given neuron can be 
%            visualized by:
%            figure();imagesc(reshape(F(:,ii),[N,N]);      
%    
% Created: 2016/03/04, 24 Cummington, Boston
%  Byron Price
% Updated: 2017/07/17
% By: Byron Price

% CREATE GIANT CONVOLUTION MATRIX TO TAKE DISCRETE
%  LAPLACIAN OF THE SPATIAL RECEPTIVE FIELD f for each
%  neuron

% assume we sample at some frequency, Fs, and have a point process
%  representation of the spiking behavior of a given neuron, we 
%  have the set of stimuli used in the input matrix S and the 
%  timeStamps for when each was displayed
%  If neuron i fired after the onset of stimulus j, then we'll take 
%  the number of spikes in a certain window after the onset of that
%  stimulus as the response, with one caveat: we'll take the spikes 
%  at the offset of the stimulus and count those as negative responses.
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

S_f(S_f==inf) = 1;
a = 0;b = 255;
unbiasedS = zeros(size(S));
for ii=1:numStimuli
    temp = reshape(real(ifftn(fftn(double(reshape(S(ii,:),[DIM(1),DIM(2)])))./S_f)),[DIM(1)*DIM(2),1])';
    currentMin = min(temp);currentMax = max(temp);
    temp = ((b-a).*(temp-currentMin))/(currentMax-currentMin)+a;
    unbiasedS(ii,:) = temp;
end

if strcmp(NoiseType,'pinkHC') == 1
   unbiasedS(unbiasedS<60) = 0;unbiasedS(unbiasedS>=60 & unbiasedS<196) = 127;
   unbiasedS(unbiasedS>=196) = 255;
end

clear S S_f U V u v;

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
    X Y totalMillisecs reducedMov movement svStrobed;

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

% REGULARIZED PSEUDO-INVERSE SOLUTION
fullSize = DIM(1)*DIM(2);
L = zeros(fullSize,fullSize,'single');

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

temp = randperm(numStimuli);
divide = round(0.7*numStimuli);
train = temp(1:divide);test = temp(divide+1:end);clear temp divide;

numLambda = 10;
loglambda = logspace(3,6,numLambda);
F = zeros(totalUnits,fullSize);
bestLambda = zeros(totalUnits,1);
for zz=1:totalUnits
   fprintf('Running unit %d ...\n',zz);
   spikeTrain = squeeze(reducedSpikeCount(zz,:,:));
   spikeTrain = sum(spikeTrain(:,50:300),2);
   r = [spikeTrain(train);zeros(fullSize,1)];
   
   tempF = zeros(numLambda,fullSize);
   RMS = zeros(numLambda,1);
   for jj=1:numLambda
      constraints = [unbiasedS(train,:);loglambda(jj).*L];
      fhat = constraints\r;tempF(jj,:) = fhat;
      RMS(jj) = norm(spikeTrain(test)-unbiasedS(test,:)*fhat);
   end
   [~,bestMap] = min(RMS);
   F(zz,:) = tempF(bestMap,:);
   bestLambda(zz) = loglambda(bestMap);
end

fileName = strcat(EphysFileName(1:end-9),'-PseudoInvResults.mat');
save(fileName,'F','totalUnits','bestLambda',...
    'reducedSpikeCount','DIM','unbiasedS','movement',...
    'xaxis','yaxis','reducedMov','allts','strobeData','totalMillisecs',...
    'svStrobed');

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

% SIMULATE DATA and test algorithm
% numStimuli = 5000;
% N = 48;
% newS = randn([numStimuli,N*N]);
% gaborFilter = @(x,y) exp(-x.^2./(2*3*3)-y.^2./(2*3*3)).*sin(2*pi*0.05.*(x.*cos(pi/4)+y.*sin(pi/4)));
% gaussFilter = @(x,y) exp(-x.^2./(2*3*3)-y.^2./(2*3*3));
% x = linspace(-20,20,N);y = linspace(-20,20,N);
% [X,Y] = meshgrid(x,y);
% gabor = gaborFilter(X,Y);
% gauss = gaussFilter(X,Y);
% r = zeros(numStimuli,1);
% for ii=1:numStimuli
%     temp = newS(ii,:);
%     temp = temp-min(temp);temp = (temp./max(temp)).*255;
%     temp = temp-(mean(temp)-127);
%     newS(ii,:) = temp-127;
%     tempIm = reshape(temp,[N,N]);
%     gaussOutput = 1;%sum(sum(conv2(tempIm,gauss)));
%     gaborOutput = sum(newS(ii,:)'.*gabor(:));
%     lambda = exp((gaborOutput/gaussOutput)./(N*N));
%     r(ii) = poissrnd(lambda);
% end
% 
% L = zeros(N*N,N*N,'single');
% 
% %operator = [0,-1,0;-1,4,-1;0,-1,0];
% bigCount = 1;
% for jj=1:N
%     for ii=1:N
%         tempMat = zeros(N,N);
%         tempMat(ii,jj) = 4;
%         if ii > 1
%             tempMat(ii-1,jj) = -1;
%         end
%         if ii < N
%             tempMat(ii+1,jj) = -1;
%         end
%         if jj > 1
%             tempMat(ii,jj-1) = -1;
%         end
%         if jj < N
%             tempMat(ii,jj+1) = -1;
%         end
%         L(bigCount,:) = tempMat(:)';bigCount = bigCount+1;
%     end
% end
% bigLambda = [0,5e1,1e2,1e3,1e4,5e4,1e5,1e6,1e7,1e8];
% RMS = zeros(length(bigLambda),1);
% parfor ii=1:length(bigLambda)
%     A = [newS;bigLambda(ii).*L];
%     constraints = [r;zeros(N*N,1)];
%     fhat = pinv(A)*constraints;
%     RMS(ii) = (N*N)^(-0.5)*norm(r-newS*fhat);
% end
% rmsDiff = diff(RMS);lambdaDiff = diff(bigLambda)';
% deltaRMSdeltaLambda = rmsDiff./lambdaDiff;
% [maxVal,index] = max(abs(deltaRMSdeltaLambda));
% onepercent = 0.01*maxVal;
% firstBelow = find(abs(deltaRMSdeltaLambda(index:end))<onepercent,1);
% bestMap = index+firstBelow-1;
% 
% A = [newS;bigLambda(bestMap).*L];
% fhat = pinv(A)*constraints;
% figure();subplot(2,1,1);imagesc(gabor);
% subplot(2,1,2);imagesc(reshape(fhat,[N,N]));
% title(sprintf('Lambda: %3.0e',bigLambda(bestMap)));
