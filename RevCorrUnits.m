function [] = RevCorrUnits(AnimalName,Date,NoiseType)
%RevCorrUnits.m
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
a = 0;b = 255;
unbiasedS = zeros(size(S));
for ii=1:numStimuli
    temp = reshape(real(ifft2(fft2(double(reshape(S(ii,:),[DIM(1),DIM(2)]))).*S_f)),[DIM(1)*DIM(2),1])';
    currentMin = min(temp);currentMax = max(temp);
    temp = ((b-a).*(temp-currentMin))/(currentMax-currentMin)+a;
    unbiasedS(ii,:) = temp;
end

if strcmp(NoiseType,'pinkHC') == 1
   unbiasedS(unbiasedS<60) = 0;unbiasedS(unbiasedS>=60 & unbiasedS<196) = 127;
   unbiasedS(unbiasedS>=196) = 255;
end

clear U V u v;

% DIVIDE DATA INTO TEST AND TRAIN
if exist('test','var') == 1 && exist('train','var') == 1
    fprintf('Train and test predetermined.\n');
else
    % divide data into training and test
   temp = randperm(numStimuli);
   divide = round(0.75*numStimuli);
   train = temp(1:divide);test = temp(divide+1:end);clear temp divide;
end

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
    X Y totalMillisecs reducedMov movement svStrobed S S_f test train beta;

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
        
        if ii==1 && jj==1
            tempMat(ii,jj) = 4;
            tempMat(ii+1,jj) = -1;
            tempMat(ii,jj+1) = -1;
        elseif ii==DIM(1) && jj==1
            tempMat(ii,jj) = 4;
            tempMat(ii-1,jj) = -1;
            tempMat(ii,jj+1) = -1;
        elseif ii==1 && jj==DIM(2)
            tempMat(ii,jj) = 4;
            tempMat(ii,jj-1) = -1;
            tempMat(ii+1,jj) = -1;
        elseif ii == DIM(1) && jj == DIM(2)
            tempMat(ii,jj) = 4;
            tempMat(ii-1,jj) = -1;
            tempMat(ii,jj-1) = -1;
        elseif ii==1
            tempMat(ii,jj) = 4;
            tempMat(ii,jj-1) = -1;
            tempMat(ii+1,jj) = -1;
            tempMat(ii,jj+1) = -1;
        elseif jj==1
            tempMat(ii,jj) = 4;
            tempMat(ii-1,jj) = -1;
            tempMat(ii,jj+1) = -1;
            tempMat(ii+1,jj) = -1;
        elseif ii==DIM(1)
            tempMat(ii,jj) = 4;
            tempMat(ii-1,jj) = -1;
            tempMat(ii,jj-1) = -1;
            tempMat(ii,jj+1) = -1;
        elseif jj==DIM(2)
            tempMat(ii,jj) = 4;
            tempMat(ii-1,jj) = -1;
            tempMat(ii+1,jj) = -1;
            tempMat(ii,jj-1) = -1;
        else
            tempMat(ii,jj) = 4;
            tempMat(ii-1,jj) = -1;
            tempMat(ii+1,jj) = -1;
            tempMat(ii,jj+1) = -1;
            tempMat(ii,jj-1) = -1;
        end
        L(bigCount,:) = tempMat(:)';bigCount = bigCount+1;
%         imagesc(tempMat);caxis([-1 4]);pause(0.1);
    end
end

warning('off','all');

inverseSf = 1./S_f;
inverseSf(inverseSf==inf) = 0;

matrix = zeros(DIM(1),DIM(2));
matrix(1,:) = 1;matrix(end,:) = 1;matrix(:,1) = 1;matrix(:,end) = 1;
edgeInds = matrix(:)==1;

fullSizeRevised = sum(~edgeInds);

numLambda = 25;
loglambda = logspace(2,6,numLambda);
F = zeros(totalUnits,fullSizeRevised);
STA = zeros(totalUnits,fullSize);

if fullSize<=10000
    STC = zeros(totalUnits,fullSize,fullSize);
else
    STC = 0;
end

bestLambda = zeros(totalUnits,1);
heldOutDeviance = zeros(totalUnits,5);
heldOutExplainedVariance = zeros(totalUnits,1);
sigmoidNonlin = zeros(totalUnits,4);
for zz=1:totalUnits
   fprintf('Running unit %d ...\n',zz);
   spikeTrain = squeeze(reducedSpikeCount(zz,:,:));
   
   baseRate = sum(sum(spikeTrain(:,51:400)))./(numStimuli*0.35);
   spikeTrain = sum(spikeTrain(:,50:200),2);
   
   % PCA solution
%    [V,D] = eig(cov(S));
%    eigenvals = diag(D);
%    start = find(eigenvals>10,1,'first');
%    q = size(eigenvals(start:end),1);
%    meanEig = mean(eigenvals(1:start-1));
%    W = V(:,start:end)*sqrtm(D(start:end,start:end)-meanEig.*eye(q));
%    W = fliplr(W);
%    Winv = pinv(W);
%    x = zeros(q,numStimuli);
%    for ii=1:numStimuli
%       x(:,ii) = Winv*(S(ii,:)'-127); 
%    end
%    x = x';
%    fhat = x(train,:)\spikeTrain(train); % then proceed as usual, with no smoothing
   
   
   r = [spikeTrain(train);sparse(fullSize,1)];
   
   tempF = zeros(numLambda,fullSizeRevised);
   tempSigmoid = zeros(numLambda,4);
   tempDev = zeros(numLambda,4);
%    allOnesTrain = ones(length(train),1);
%    allOnesTest = ones(length(test),1);
   for jj=1:numLambda
       % calculate regularized pseudoinverse solution
      constraints = [unbiasedS(train,:);loglambda(jj).*L];
      fhat = constraints\r;fhat = full(fhat);
      
      % remove bias created by using pink noise
      temp = reshape(fhat,[DIM(1),DIM(2)]);
      temp = real(ifft2(fft2(temp).*S_f));
      fhat = temp(:);fhat = fhat(~edgeInds);
      
      % fit sigmoid nonlinearity
      result = S(train,~edgeInds)*fhat;
      myFun = @(x) x(1)./(1+exp(-(result-x(2)).*x(3)))+x(4)-spikeTrain(train);
      x0 = [5,0.5,10,baseRate];
      lb = [0,0,0,0];ub = [5e2,Inf,Inf,Inf];
%       options.Algorithm = 'levenberg-marquardt';
      %options.MaxFunctionEvaluations = 1000;
      options = optimoptions(@lsqnonlin,'MaxIterations',1000);
      options.Display = 'off';
      sigmoidParams = lsqnonlin(myFun,x0,lb,ub,options);
      
      temp = sigmoidParams(1)./(1+exp(-(result-sigmoidParams(2)).*sigmoidParams(3)))+sigmoidParams(4);
      initialDev = spikeTrain(train).*log(spikeTrain(train)./temp)-(spikeTrain(train)-temp);
      initialDev(isnan(initialDev) | isinf(initialDev)) = temp(isnan(initialDev) | isinf(initialDev));
      sigmoidParams = FitSigmoid(sigmoidParams,2*sum(initialDev),spikeTrain(train),result);
      
      % save parameters
      tempF(jj,:) = fhat;
      tempSigmoid(jj,:) = sigmoidParams;
      
      % calculate held-out Poisson deviance
      temp = S(test,~edgeInds)*fhat;
      temp = sigmoidParams(1)./(1+exp(-(temp-sigmoidParams(2)).*sigmoidParams(3)))+sigmoidParams(4);
      
      initialDev = spikeTrain(test).*log(spikeTrain(test)./temp)-(spikeTrain(test)-temp);
      initialDev(isnan(initialDev) | isinf(initialDev)) = temp(isnan(initialDev) | isinf(initialDev));
      tempDev(jj,1) = 2*sum(initialDev);
      
      % rotate the receptive field and recalculate the held-out deviance
      rf = reshape(fhat,[DIM(1)-2,DIM(2)-2]);
      for kk=2:4
        rf = imrotate(rf,90);
        tempfhat = rf(:);
        temp = S(test,~edgeInds)*tempfhat;
        temp = sigmoidParams(1)./(1+exp(-(temp-sigmoidParams(2)).*sigmoidParams(3)))+sigmoidParams(4);
        
        initialDev = spikeTrain(test).*log(spikeTrain(test)./temp)-(spikeTrain(test)-temp);
        initialDev(isnan(initialDev) | isinf(initialDev)) = temp(isnan(initialDev) | isinf(initialDev));
        tempDev(jj,kk) = 2*sum(initialDev);
      end
   end
  
   [~,bestMap] = min(tempDev(:,1));
   fprintf('Best Map: %d\n',bestMap);
   F(zz,:) = tempF(bestMap,:);
   bestLambda(zz) = loglambda(bestMap);
   heldOutDeviance(zz,1:4) = tempDev(bestMap,:);
   sigmoidNonlin(zz,:) = tempSigmoid(bestMap,:);
   
   % get held-out deviance of model with single regressor
   [~,dev,~] = glmfit(ones(length(test),1),spikeTrain(test),'poisson','constant','off');
   heldOutDeviance(zz,5) = dev;
   
   % get held-out explained variance
   %  see R-squared measures for Count Data Regression Models
   %   A. Colin Cameron & Frank A.G. Windmeijer April 1995
   %   Journal of Business and Economic Statistics
   heldOutExplainedVariance(zz,1) = 1-tempDev(bestMap,1)/dev;
   fprintf('Fraction of Explained Variance: %3.2f\n\n',1-tempDev(bestMap,1)/dev);
   
   totalSpikes = sum(spikeTrain);
   for ii=1:numStimuli
      if spikeTrain(ii) > 0
         STA(zz,:) = STA(zz,:)+S(ii,:).*(spikeTrain(ii)/totalSpikes);
      end
   end
   
   if fullSize<=10000
       for ii=1:numStimuli
           if spikeTrain(ii) > 0
               difference = S(ii,:)-STA(zz,:);
               STC(zz,:,:) = STC(zz,:,:)+(1/(totalSpikes-1))*...
                   (spikeTrain(ii)).*(difference*difference');
           end
       end
   end
end

% save the results
fileName = strcat(EphysFileName(1:end-9),'-PseudoInvResults.mat');
save(fileName,'F','totalUnits','bestLambda',...
    'reducedSpikeCount','DIM','unbiasedS','movement',...
    'xaxis','yaxis','reducedMov','allts','strobeData','totalMillisecs',...
    'svStrobed','heldOutDeviance','numStimuli','S_f','sigmoidNonlin',...
    'heldOutExplainedVariance','STA','STC','beta');

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

function [sigmoidParams] = FitSigmoid(sigmoidParams,initialDev,spikeTrain,result)
N = 5e4;
numParams = length(sigmoidParams);
paramVec = zeros(numParams,N);
devVec = zeros(N,1);

variance = 0.01;

paramVec(:,1) = sigmoidParams';
devVec(1) = initialDev;

W = diag(ones(numParams,1));
loglambda = log(2.38^2)*ones(numParams,1);
optimalAccept = 0.44;
for ii=2:N
    index = random('Discrete Uniform',numParams);
    
    pStar = paramVec(:,ii-1)+W(:,index).*normrnd(0,sqrt(exp(loglambda(index))*variance));
    
    if sum(pStar<=0) == 0
        temp = pStar(1)./(1+exp(-(result-pStar(2)).*pStar(3)))+pStar(4);
        
        pDev = spikeTrain.*log(spikeTrain./temp)-(spikeTrain-temp);
        pDev(isnan(pDev) | isinf(pDev)) = temp(isnan(pDev) | isinf(pDev));
        
        logA = devVec(ii-1) - 2*sum(pDev);
        if log(rand) < logA
           paramVec(:,ii) = pStar;
           devVec(ii) = 2*sum(pDev);
        else
            paramVec(:,ii) = paramVec(:,ii-1);
            devVec(ii) = devVec(ii-1);
            
        end
        loglambda(index) = loglambda(index)+0.01*(exp(min(0,logA))-optimalAccept);
    else
       loglambda(index) = loglambda(index)+0.01*(-optimalAccept);
       paramVec(:,ii) = paramVec(:,ii-1);
       devVec(ii) = devVec(ii-1);      
    end
end

[~,ind] = min(devVec);
sigmoidParams = paramVec(:,ind);

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
