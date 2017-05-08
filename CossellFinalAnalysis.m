% CossellFinalAnalysis.m

fileName = {'13_05_22','13_05_20','13_05_19','13_05_14',...
    '13_05_13','13_05_10','13_05_07','13_05_06','13_05_01','13_04_30',...
    '13_04_29','13_04_26','13_04_21','13_04_17','13_04_16','13_04_14',...
    '13_04_05'};
numExps = length(fileName);

numImPlanes = 8;

load('NaturalImageStack.mat');
numIms = size(tfImageStack,3);

for ii=1:numExps
    for jj=1:numImPlanes
       tempFile = strcat(fileName{ii},'-',num2str(jj),'.mat');
       load(tempFile);
       numNeurons = size(F,1);
       evokedSlot = 61;
       numLambdas = size(F,2);
       
       correlationCoeffs = zeros(numNeurons,numLambdas);
       bestFilters = zeros(numNeurons,size(F,3));
       bestCoeffs = zeros(numNeurons,2);
       train = round(numIms*0.7);
       for kk=1:numNeurons
           r = reducedSpikeData(kk,evokedSlot:evokedSlot+numIms-1)';
           for ll=1:numLambdas
               tempF = reshape(squeeze(F(kk,ll,:)),[DIM(1),DIM(2)]);
               filterOutput = zeros(numIms,1);
               for mm=1:numIms
                   filterOutput(mm) = sum(sum(tempF.*tfImageStack(:,:,mm)));
               end
               y = r(1:train);
               d = filterOutput(1:train);
               myFun = @(x) x(1).*exp((d-x(2)).*x(3))./(1+exp((d-x(2)).*x(3)))-y;
               x0 = [max(r),mean(d),1e3];
               
               x = lsqnonlin(myFun,x0);
               
               tempFun = @(x,d) x(1).*exp((d-x(2)).*x(3))./(1+exp((d-x(2)).*x(3)));
               d = filterOutput(train+1:end);
               actualResponse = r(train+1:end);
               predictedResponse = tempFun(x,d);
               tempr = corrcoef(actualResponse,predictedResponse);
               correlationCoeffs(kk,ll) = tempr(1,2);
           end
           [bestCoeffs(kk,2),bestCoeffs(kk,1)] = max(correlationCoeffs(kk,:));
           bestFilters(kk,:) = squeeze(F(kk,bestCoeffs(kk,1),:));
       end
       saveFile = strcat(fileName{ii},'-',num2str(jj),'final.mat');
       save(saveFile,'correlationCoeffs','bestFilters','bestCoeffs');
    end
end