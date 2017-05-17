% STRF_SimulateResponses.m

% INITIALIZE VARIABLES AND MAKE SEQUENCE VIDEO, binned at 1ms
horzPix = 2560;vertPix = 1440;
mmPerPixel = 0.2363;
spatFreq = 0.05; % cycles per degree
DistToScreen = 250; % mm
orientations = [60,150,90,120].*pi./180;
sequenceTime = 150; %ms

temp = (tan((1/spatFreq)*pi/180)*(DistToScreen))/mmPerPixel;
spatFreq = 1/temp;

xaxis = 1:10:2560;yaxis = 1:10:1440;
xaxis = (xaxis-median(xaxis));
yaxis = (yaxis-median(yaxis));

xaxis = xaxis(128-35:128+35);yaxis = yaxis(72-35:72+35);

xLen = length(xaxis);
yLen = length(yaxis);
sequenceVid = zeros(yLen,xLen,sequenceTime*4+1200);
totalTime = size(sequenceVid,3);
[XPOS,YPOS] = meshgrid(xaxis,yaxis);
count = 601;
for ii=1:4
    orientVec = cos(orientations(ii)-pi/2).*XPOS+sin(orientations(ii)-pi/2).*YPOS;
    Z = sin(2*pi*spatFreq.*orientVec);
    for jj=1:150
        sequenceVid(:,:,count) = Z;
        count = count+1;
    end
end

% DISPLAY THE SEQUENCE STIMULUS AS A VIDEO
%   on my computer, displays slowly (not in real-time) as the computer
%   obviously cannot do 1000 Hz
% figure();pause(1);
% for ii=1:totalTime
%    imagesc(sequenceVid(:,:,ii));caxis([-1 1]);set(gca,'YDir','normal');pause(1/1000); 
% end
% 

% gaborFun = @(x,y,t,xc,yc,tc,sigma,tau,gamma,spatFreq,v,theta,phi) ...
%     (gamma./(sqrt(2*pi*sigma*sigma))).*exp(-((x-xc).^2+gamma.*gamma.*(y-yc).^2)/(2*sigma*sigma))...
%     .*sin((2*pi.*spatFreq).*(cos(theta-pi/2).*(x-xc)+sin(theta-pi/2).*(y-yc)+v.*t)+phi)...
%     .*(1/(sqrt(2*pi*tau))).*exp(-(t-tc).^2/(2*tau*tau));

% gaborFun = @(x,y,t,xc,yc,tc,sigmax,sigmay,sigmat,spatFreq,v,theta,phi) ...
%     exp(-(x-xc).^2./(2.*sigmax.^2)-(y-yc).^2./(2.*sigmay.^2)-(t-tc).^2./(2.*sigmat.^2))...
%     .*sin((2*pi.*spatFreq).*(cos(theta-pi/2).*(x-xc)+sin(theta-pi/2).*(y-yc)+v.*t)+phi);

% DEFINE SPATIOTEMPORAL RECEPTIVE FIELD (STRF) AND PLOT EXAMPLE CONVOLUTION
%   OF STRF WITH SEQUENCE STIMULUS
%  the temporal component is a zero-gain, time-causal filter, which
%  endows the "neuron" with its preferred temporal frequency, its transient
%  responses to visual images, and its latency
%   if the v (velocity) is zero, then this is a space-time separable
%   filter, but if v ~= 0, it's space-time inseparable
gaborFun = @(x,y,t,xc,yc,sigma,k,n,gama,spatFreq,v,theta,phi) ...
    exp(-((x-xc).^2+gama.*gama.*(y-yc).^2)/(2*sigma*sigma))...
    .*sin((2*pi.*spatFreq).*(cos(theta-pi/2).*(x-xc)+sin(theta-pi/2).*(y-yc)+v.*t)+phi)...
    .*(k.*t).^n.*exp(-k.*t).*(1/gamma(n+1)-(k.*t).^2./(gamma(n+3)));

t = 300:-1:0;
sigma = 100; % size (arbitrary units) of the spatial Gabor
k = 0.4; % k and n define the temporal filter
n = 20;
gama = 1; % leave at 1, would stretch the extent of the Gabor in y dimension
theta = 60*pi/180; % neuron's preferred orientation
phi = 0; % phase of the sinusoid
v = 0; % velocity (shift in the phase over time)


[X,Y,T] = meshgrid(xaxis,yaxis,t);
filter1 = gaborFun(X,Y,T,0,0,sigma,k,n,gama,spatFreq,v,theta,phi);

% DISPLAY MODEL STRF AS A VIDEO
% minVal = min(filter1(:));
% maxVal = max(filter1(:));
%filter1(abs(filter1)<1e-4) = 0;
% figure();
% for ii=400:-1:150
%     imagesc(filter1(:,:,ii));set(gca,'YDir','normal');
%     colormap('bone');caxis([minVal maxVal]);pause(1/1000);
% end

% displays the temporal component of the filter reversed from 
%  its normal orientation in time ... the fact
%  that the filter peaks early and then decays signifies that this causal
%  filter weights the recent past more strongly than it weights the distant
%  past

% figure();
% tempT = 0:150;
% plot(tempT,(k.*tempT).^n.*exp(-k.*tempT).*(1/factorial(n)-(k.*tempT).^2./factorial(n+2)));

filter1Energy = sum(sum(sum(filter1.*filter1)));


convResult = zeros(900,1);
count = 1;
for ii=500:500+900-1
    temp = sequenceVid(:,:,ii-300:ii);

    temp1 = temp.*filter1./filter1Energy;
    convResult(count,1) = sum(temp1(:));

    count = count+1;
end

figure();
time = -100:(800-1);
plot(time',convResult(:,1));title('Example Convolution of STRF with Sequence Stimulus');
set(gca,'XTickLabel',[-100,0,100,200,300,400,500,600,700]);
xlabel('Time from stimulus onset (ms)');


% CALCULATE CONVOLUTION FOR A SET OF STRFs WITH RANDOM (BUT REASONABLE)
%  PROPERTIES ... convolution is somewhat slow, try fewer "neurons" maybe
numCells = 20;
paramVec = zeros(numCells,6);
allFilters = cell(numCells,1);

for ii=1:numCells
    % still haven't figured out all the intricacies of the temporal
    %  filter function, so need to randomly select values but make sure
    %  they are reasonable (proper latency, temporal frequency not too
    %  high)
   latency = 0;timing = 0;maxVal = 500;
   while latency<40 || latency>90 || timing < 50 || maxVal > 0.3
       paramVec(ii,1) = 0.5+normrnd(0,0.2); % k
       paramVec(ii,2) = 20+normrnd(0,10); % n
       tempT = 0:250;
       temp = (paramVec(ii,1).*tempT).^paramVec(ii,2)...
           .*exp(-paramVec(ii,1).*tempT).*(1/gamma(paramVec(ii,2)+1)-(paramVec(ii,1).*tempT).^2./gamma(paramVec(ii,2)+3));
       tempT = diff(sign(temp));
       latency = single(find(tempT==-2));
       if isempty(latency)==1
           latency = 1;
       end
       [maxVal,maxInd] = max(temp);
       [~,minInd] = min(temp);
       timing = minInd-maxInd;
       paramVec(ii,6) = timing;
   end
   
   paramVec(ii,3) = normrnd(0,0.1); % velocity
   paramVec(ii,4) = rand*2*pi; % orientation theta
   paramVec(ii,5) = rand*2*pi-pi; % phase phi

   filterOutput = gaborFun(X,Y,T,0,0,sigma,paramVec(ii,1),paramVec(ii,2),gama,spatFreq,paramVec(ii,3),paramVec(ii,4),paramVec(ii,5));
   
   allFilters{ii,1} = filterOutput./sum(sum(sum(filterOutput.*filterOutput)));
end

% perform convolution for all cells ... slide the filter across
%  the video and multiply at each time point
convResult = zeros(900,numCells);

count = 1;
for ii=500:500+900-1
    temp = sequenceVid(:,:,ii-300:ii);
    singleCellOutput = zeros(1,numCells);
    for jj=1:numCells
        temp1 = allFilters{jj,1}.*temp;
        singleCellOutput(jj) = sum(temp1(:));
    end
    convResult(count,:) = singleCellOutput;
    count = count+1;
end



figure();
plot(time,mean(real(convResult),2));
minVal = min(mean(convResult,2));
maxVal = max(mean(convResult,2));
axis([time(1) time(end) minVal-10 maxVal+10]);
hold on;plot(0.*ones(10,1),linspace(minVal-10,maxVal+10,10),'k--','LineWidth',2);
plot(150.*ones(10,1),linspace(minVal-10,maxVal+10,10),'k--','LineWidth',2);
plot(300.*ones(10,1),linspace(minVal-10,maxVal+10,10),'k--','LineWidth',2);
plot(450.*ones(10,1),linspace(minVal-10,maxVal+10,10),'k--','LineWidth',2);
title('Averaged Response to Sequence Stimulus');

% ADD NORMALIZATION STEP, if desired, neurons that fire strongly to a given
%  element would inhibit other neurons due to normalization
% sigma_norm = 10; % this controls how much different neurons affect each
%  other ... a value of 0 means a pure normalization, where at time each
%  point the total "activity" across the group of neurons must sum to 1
%  higher values of sigma_norm approach each neuron being independent
% newConvResult = zeros(size(convResult));
% for ii=1:900
%     temp = convResult(ii,:);
%     newConvResult(ii,:) = temp./sqrt(sigma_norm^2+sum(temp)^2);
% end
% time = -100:800-1;
% figure();
% for ii=1:25
%     subplot(5,5,ii);plot(time,newConvResult(:,ii));
%     minVal = min(newConvResult(:,ii));
%     maxVal = max(newConvResult(:,ii));
%     hold on;plot(0.*ones(10,1),linspace(minVal,maxVal,10),'k--','LineWidth',2);
% end

% VISUALIZE SIMULATED SPIKING ACTIVITY OF THE "NEURONS" WITH THE 
%   STRONGEST RESPONSES TO THE STIMULUS (top 25%)
% maxVals = max(convResult,[],1);
% threshold = quantile(maxVals(~isnan(maxVals)),0.75);
% 
% inds = find(maxVals>threshold);

% [maximumModulation,~] = max(abs(maxVals));
% inds = random('Discrete Uniform',numCells,[20,1]);
     
gaussKernel = @(t,delta) (1./(sqrt(2*pi)*delta)).*exp(-t.^2./(2*delta*delta));
delta = 2:1:20;
convolveKernel = @(t,delta) (1./(sqrt(pi)*2*delta)).*exp(-t.^2./(4*delta*delta));
     
sigmoidParam = 10;
for neuronNum = 1:10
    % limit total modulation to be 20 ... see below, exp(4) would be
    %  multiplied by the base firing rate (~20 fold increase) ... sigmoid
    %  function saturates highest possible firing rate
    targetNeuronModulation = 10.*(exp(convResult(:,neuronNum)./sigmoidParam)...
        ./(1+exp(convResult(:,neuronNum)./sigmoidParam)));
    %targetNeuronModulation = conv(targetNeuronModulation,gausswin(51,5),'same');
    
    % SIMULATE SPIKE TRAIN BASED ON A HISTORY-DEPENDENT GENERAL POINT
    %   PROCESS MODEL
    %  neurons spike at a base rate of 1 Hz (1/1000 given the 1ms time
    %  bins) and their activity is modulated by past spiking and also by
    %  the convolution of their STRFs with the sequence stimulus
    
    % add history dependence to the spike trains, slightly more realistic,
    %  this particular history dependence yields a refractory period and 
    %  a tendency to burst
    
    % -1 the first bin, i.e. 1ms ago
    historyB = [-1,-0.5,-0.3,-0.1,0,0.1,0.15,0.3,0.3,0.25,0.2,0.15,0.1,0.05,0.025]';
    
    trials = 200;
    time = -100:800-1;
    spikeTrain = zeros(trials,length(time));
    baseRate = 3/1000;
    
    for jj=1:trials
        currentHistory = zeros(1,length(historyB));
        for ii=1:length(time)
            % lambda defines the rate parameter for the generation of a
            %  poisson-distributed random variable
            lambda = baseRate*exp(currentHistory*historyB)*targetNeuronModulation(ii);
            spikeTrain(jj,ii) = poissrnd(lambda);
            currentHistory = [spikeTrain(jj,ii),currentHistory(1:end-1)];
        end
    end
    figure();subplot(2,1,1);imagesc(time,1:trials,~spikeTrain);colormap('bone');caxis([0 1]);
    ylabel('Trial');xlabel('Time from Stimulus Onset (ms)');
    title(sprintf('%3.0f %3.0f %3.0f %3.0f Neuron Orientation: %3.0f',...
        orientations(1)*180/pi,orientations(2)*180/pi,orientations(3)*180/pi,orientations(4)*180/pi,paramVec(neuronNum,4)*180/pi));
    hold on;plot(0.*ones(10,1),linspace(0,trials,10),'k--');
    plot(150.*ones(10,1),linspace(0,trials,10),'k--');
    plot(300.*ones(10,1),linspace(0,trials,10),'k--');
    plot(450.*ones(10,1),linspace(0,trials,10),'k--');
    
    % calculate optimal bandwidth delta of the Gaussian kernel density
    %  estimator
    % optimal estimate found by minimizing the mean integrated squared
    % error (MISE) ... from Shigeru Shinomoto "Estimating the Firing Rate"
    %  Chapter 2 from the book Analysis of Parallel Spike Trains
    spikeTimes = find(spikeTrain==1);
    [~,J] = ind2sub(size(spikeTrain),spikeTimes);
    cost = zeros(length(delta),1);
    
    for ii=1:length(delta)
        tempCost1 = 0;
        tempCost2 = 0;
        for jj=1:length(J)
            inds = find(abs(J(jj)-J)<(5*delta(ii)));
            for kk=1:length(inds)
                if jj ~= inds(kk)
                    tempCost2 = tempCost2+gaussKernel(J(jj)-J(inds(kk)),delta(ii));
                end
                tempCost1 = tempCost1+convolveKernel(J(jj)-J(inds(kk)),delta(ii));
            end
        end
        cost(ii) = (1/trials^2)*tempCost1-(2/trials^2)*tempCost2;
    end
    [~,ind] = min(cost);
    finalDelta = delta(ind);
    
    binned = sum(spikeTrain,1);
    binned = conv(binned,gaussKernel(round(-5*finalDelta):1:round(5*finalDelta),finalDelta),'same');
    binned = 1000.*binned./trials;
    subplot(2,1,2);plot(time,binned,'LineWidth',2);title('Post-stimulus Kernel Density Estimate');
    xlabel('Time');ylabel('Firing Rate (Hz)');
    axis([-100 800-1 0 max(binned)+5]);
end
