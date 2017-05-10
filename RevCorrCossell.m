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
 
%  unbiased_convertedStack = zeros(size(tfImageStack,3),size(tfImageStack,1)*size(tfImageStack,2));
%  for ii=1:size(unbiasedImageStack,3)
%      tempIm = unbiasedImageStack(:,:,ii);
%      unbiased_convertedStack(ii,:) = tempIm(:)';
%  end
 
 unbiasedPink_convertedStack = zeros(size(tfImageStack,3),size(tfImageStack,1)*size(tfImageStack,2));
 for ii=1:size(unbiasedPinkStack,3)
     tempIm = unbiasedPinkStack(:,:,ii);
     unbiasedPink_convertedStack(ii,:) = tempIm(:)';
 end
 
 clear tfImageStack unbiasedPinkStack unbiasedImageStack;
 
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
        
        bigLambda = logspace(3,6,10);
        RMS = zeros(numNeurons,length(bigLambda));
        F = zeros(numNeurons,length(bigLambda),DIM(1)*DIM(2));
        
        train = round(numIms*0.75);
        for kk=1:numNeurons
            r = reducedSpikeData(kk,evokedSlots(1):evokedSlots(1)+train-1)';
            constraints = [r;zeros(DIM(1)*DIM(2),1)];
            tempA = convertedStack(1:train,:);
            
            rTest = reducedSpikeData(kk,evokedSlots(1)+train:evokedSlots(1)+numIms-1)';
            testA = convertedStack(train+1:numIms,:);
            for ll=1:length(bigLambda)
                A = [tempA;bigLambda(ll).*L];
                F(kk,ll,:) = A\constraints;
                
                RMS(kk,ll) = norm(rTest-testA*squeeze(F(kk,ll,:)))./sqrt(DIM(1)*DIM(2));
            end
        end
        
%         bigLambda_Unbiased = logspace(0,3,15);
%         RMS_Unbiased = zeros(numNeurons,length(bigLambda_Unbiased));
%         F_Unbiased = zeros(numNeurons,length(bigLambda_Unbiased),DIM(1)*DIM(2));
%         for kk=1:numNeurons
%             r = reducedSpikeData(kk,evokedSlots(1):evokedSlots(1)+train-1)';
%             constraints = [r;zeros(DIM(1)*DIM(2),1)];
%             tempA = unbiased_convertedStack(1:train,:);
%             
%             rTest = reducedSpikeData(kk,evokedSlots(1)+train:evokedSlots(1)+numIms-1)';
%             testA = unbiased_convertedStack(train+1:numIms,:);
%             for ll=1:length(bigLambda_Unbiased)
%                 A = [tempA;bigLambda_Unbiased(ll).*L];
%                 F_Unbiased(kk,ll,:) = A\constraints;
%                 
%                 RMS_Unbiased(kk,ll) = norm(rTest-testA*squeeze(F_Unbiased(kk,ll,:)))./sqrt(DIM(1)*DIM(2));
%             end
%         end
        
        bigLambda_UnbiasedPink = logspace(2,5,10);
        RMS_UnbiasedPink = zeros(numNeurons,length(bigLambda_UnbiasedPink));
        F_UnbiasedPink = zeros(numNeurons,length(bigLambda_UnbiasedPink),DIM(1)*DIM(2));
        for kk=1:numNeurons
            r = reducedSpikeData(kk,evokedSlots(1):evokedSlots(1)+train-1)';
            constraints = [r;zeros(DIM(1)*DIM(2),1)];
            tempA = unbiasedPink_convertedStack(1:train,:);
            
            rTest = reducedSpikeData(kk,evokedSlots(1)+train:evokedSlots(1)+numIms-1)';
            testA = unbiasedPink_convertedStack(train+1:numIms,:);
            for ll=1:length(bigLambda_UnbiasedPink)
                A = [tempA;bigLambda_UnbiasedPink(ll).*L];
                F_UnbiasedPink(kk,ll,:) = A\constraints;
                
                RMS_UnbiasedPink(kk,ll) = norm(rTest-testA*squeeze(F_UnbiasedPink(kk,ll,:)))./sqrt(DIM(1)*DIM(2));
            end
        end
        
        newFileName = strcat(fileName{ii},'-',num2str(jj));
        save(newFileName,'F','bigLambda','DIM','RMS','reducedSpikeData',...
            'RMS_UnbiasedPink','F_UnbiasedPink','bigLambda_UnbiasedPink');
    end
end