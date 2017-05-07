% RevCorr_Cossell.m

fileName = {'13_05_22','13_05_20','13_05_19','13_05_14',...
    '13_05_13','13_05_10','13_05_07','13_05_06','13_05_01','13_04_30',...
    '13_04_29','13_04_26','13_04_21','13_04_17','13_04_16','13_04_14',...
    '13_04_05'};
nExps = length(fileName);
[csCaData, sParameters, tfImageStack] = ...
   Load_Cossell_etal_2015_NaturalImageResponses(fileName);

load('NaturalImageStack.mat');
convertedStack = zeros(size(tfImageStack,3),size(tfImageStack,1)*size(tfImageStack,2));
for ii=1:size(tfImageStack,3)
   tempIm = tfImageStack(:,:,ii);
   convertedStack(ii,:) = tempIm(:)';
end

unbiased_convertedStack = zeros(size(tfImageStack,3),size(tfImageStack,1)*size(tfImageStack,2));
for ii=1:size(unbiasedImageStack,3)
   tempIm = unbiasedImageStack(:,:,ii);
   convertedStack(ii,:) = tempIm(:)';
end

for ii=1:nExps
    CaData = csCaData{ii};
    nImgPlanes = size(CaData,2);
    for jj=1:nImgPlanes
        spikeData = CaData(jj).spikes;
        nROIs = size(spikeData,1);
        nFrames = size(spikeData,3);
        nImageSlot = size(spikeData,2); % 1:60 and 1861:1920 are spontaneous activity
                        % 61:1860 are evoked activity
        spontaneousSlots = sParameters.vnSpontaneousFrames;
        evokedSlots = setdiff(1:max(spontaneousSlots),spontaneousSlots);
        imStackIndex = sParameters.vnPresentedImage;
        
        numNeurons = size(spikeData,1);
        numIms = length(evokedSlots);
        reducedSpikeData = sum(spikeData(:,:,3:5),3)./3;
        
        DIM = [size(tfImageStack,1),size(tfImageStack,2)];
        L = zeros(DIM(1)*DIM(2),DIM(1)*DIM(2));
        bigCount = 1;
        for kk=1:DIM(2)
            for ll=1:DIM(1)
                tempMat = zeros(DIM(1),DIM(2));
                tempMat(ll,kk) = 4;
                if kk > 1 && kk < DIM(2)
                    tempMat(ll,kk-1) = -1;
                    tempMat(ll,kk+1) = -1;
                end
                if ll > 1 && ll < DIM(1)
                    tempMat(ll-1,kk) = -1;
                    tempMat(ll+1,kk) = -1;
                end
                L(bigCount,:) = tempMat(:)';bigCount = bigCount+1;
            end
        end
        
        
        bigLambda = logspace(3,7,15);
        RMS = zeros(numNeurons,length(bigLambda));
        F = zeros(numNeurons,length(bigLambda),DIM(1)*DIM(2));
        
        train = round(numIms*0.75);
        for kk=1:numNeurons
            r = reducedSpikeData(kk,evokedSlots(1):evokedSlots(1)+train-1)';
            constraints = [r;zeros(DIM(1)*DIM(2),1)];
            tempA = convertedStack(1:train,:);
            
            for ll=1:length(bigLambda)
                A = [tempA;bigLambda(ll).*L];
                F(kk,ll,:) = A\constraints;clear A;
                rTest = reducedSpikeData(kk,evokedSlots(1)+train:evokedSlots(1)+numIms-1)';
                A = convertedStack(train+1:numIms,:);
                RMS(kk,ll) = norm(rTest-A*squeeze(F(kk,ll,:)))./sqrt(DIM(1)*DIM(2));
                clear A;
            end
             clear tempA A;
        end
        
        RMS_Unbiased = zeros(numNeurons,length(bigLambda));
        F_Unbiased = zeros(numNeurons,length(bigLambda),DIM(1)*DIM(2));
        for kk=1:numNeurons
            r = reducedSpikeData(kk,evokedSlots(1):evokedSlots(1)+train-1)';
            constraints = [r;zeros(DIM(1)*DIM(2),1)];
            tempA = unbiased_convertedStack(1:train,:);
            
            for ll=1:length(bigLambda)
                A = [tempA;bigLambda(ll).*L];
                F(kk,ll,:) = A\constraints;clear A;
                rTest = reducedSpikeData(kk,evokedSlots(1)+train:evokedSlots(1)+numIms-1)';
                A = unbiased_convertedStack(train+1:numIms,:);
                RMS(kk,ll) = norm(rTest-A*squeeze(F(kk,ll,:)))./sqrt(DIM(1)*DIM(2));
                clear A;
            end
             clear tempA A;
        end
        
        newFileName = strcat(fileName{ii},'-',num2str(jj));
        save(newFileName,'F','bigLambda','DIM','RMS','reducedSpikeData',...
            'RMS_Unbiased','F_Unbiased');
    end
end