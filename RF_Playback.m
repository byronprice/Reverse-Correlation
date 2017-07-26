function [] = RF_Playback(AnimalName,Date,NoiseType)
% RF_Playback.m

fileName = sprintf('NoiseData%s%d_%d-PseudoInvResults.mat',NoiseType,Date,AnimalName);
load(fileName,'F','DIM','totalUnits');

load('NoiseVars.mat');
numStimuli = round(totalUnits*10+totalUnits*10*0.1);
numNoise = numStimuli-totalUnits*10;


startEXP = 254;
endEXP = 255;

tcpipClient = tcpip('128.197.59.169',30000,'NetworkRole','client');
bufferSize = 50000; % bytes, (we won't need this much)
set(tcpipClient,'InputBufferSize',bufferSize);
set(tcpipClient,'Timeout',5);
fopen(tcpipClient);

% Acquire a handle to OpenGL, so we can use OpenGL commands in our code:
global GL;
usb = usb1208FSPlusClass;

% Make sure this is running on OpenGL Psychtoolbox:
AssertOpenGL;

TimeEstimate = numStimuli*(flipInterval+WaitTime+0.5)/60+4*0.5;
fprintf('\nEstimated time is %3.2f minutes.',TimeEstimate);
WaitSecs(1);

% Choose screen with maximum id - the secondary display:
screenid = max(Screen('Screens'));

% % Open a fullscreen onscreen window on that display, choose a background
% % color of 127 = gray with 50% max intensity; 0 = black;255 = white
background = 127;
[win,~] = Screen('OpenWindow', screenid,background);

gammaTable = makeGrayscaleGammaTable(gama,0,255);
Screen('LoadNormalizedGammaTable',win,gammaTable);

% Query window size in pixels
[w_pixels,h_pixels] = Screen('WindowSize', win);
minPix = min(w_pixels,h_pixels);
maxPix = max(w_pixels,h_pixels)-200;

% screen size in millimeters and a conversion factor to get from mm to pixels
[w_mm,h_mm] = Screen('DisplaySize',screenid);

% a bit confusing, we want the stimuli produced to have a certain number of
%  effective pixels, which project to larger squares of on-screen pixels
maxPix = maxPix-mod(maxPix,screenPix_to_effPix);
minPix = minPix-mod(minPix,screenPix_to_effPix);

conv_factor = (w_mm/w_pixels+h_mm/h_pixels)/2;


% GENERATION OF NOISE
if strcmp(NoiseType,'white') == 1
    beta = 0;
%     spaceExp = 2;
elseif strcmp(NoiseType,'pink') == 1
    beta = -1.5;
%     spaceExp = 2;
elseif strcmp(NoiseType,'brown') == 1
    beta = -2;
%     spaceExp = 2;
else 
    display('NoiseType must be ''white'', ''pink'' or ''brown'' as a string.')
    return;
end

Grey = 127;

shuffled = randperm(numStimuli);
greyIms = shuffled(1:numNoise);
rfIms = shuffled(numNoise+1:end);clear shuffled;

stimulusNumber = zeros(numStimuli,1);
S = zeros(numStimuli,DIM(1)*DIM(2),'uint8');
% below pink noise from Jon Yearsley, 1/f noise generate spatial data

for ii=1:numNoise
    X = spatialPattern([DIM(2),DIM(1)],beta);
    X = X-min(min(X));
    X = (X./max(max(X))).*255;
    Y = X(:);
    meanVal = mean(Y);difference = meanVal-Grey;
    %figure();imagesc(reshape(Y-difference,DIM));
    S(greyIms(ii),:) = Y-difference;
    stimulusNumber(greyIms(ii),1) = 1;
end

matrix = zeros(DIM(1),DIM(2));
matrix(1,:) = 1;matrix(end,:) = 1;matrix(:,1) = 1;matrix(:,end) = 1;
edgeInds = matrix(:)==1;clear matrix;
a = 0;b = 255;
for ii=1:totalUnits*10
    index = mod(ii,totalUnits)+1;
    temprf = F(index,:);
    temprf = reshape(temprf,[DIM(1)-2,DIM(2)-2]);
    
    currentMin = min(temprf(:));currentMax = max(temprf(:));
    temprf = ((b-a).*(temprf-currentMin))/(currentMax-currentMin)+a;
    
    fullIm = Grey.*ones(DIM(1),DIM(2));
    fullIm(~edgeInds) = temprf(:);
    
    if mod(ii,2) == 1
        fullIm = b-fullIm;
        stimulusNumber(rfIms(ii)) = (index+1)*2;
    end
    S(rfIms(ii),:) = fullIm(:);
    stimulusNumber(rfIms(ii)) = index+1;
end

S = uint8(S);

wLow = round((w_pixels-maxPix)/2);
wHigh = round(w_pixels-wLow);
hLow = round((h_pixels-minPix)/2);
hHigh = round(h_pixels-hLow);
destRect = [wLow hLow wHigh hHigh];

Priority(9);
% Retrieve monitor refresh duration
ifi = Screen('GetFlipInterval', win);

WaitTimes = WaitTime+unifrnd(0,1,[numStimuli,1]);

usb.startRecording;WaitSecs(1);usb.strobeEventWord(0);
WaitSecs(30);
usb.strobeEventWord(startEXP);WaitSecs(1);
tt = 1;
vbl = Screen('Flip',win);
while tt <= numStimuli
    if tcpipClient.BytesAvailable > 0
        data = fread(tcpipClient,tcpipClient.BytesAvailable/8,'double');
        if sum(data) > 0
            WaitSecs(4);
        else
            % Convert it to a texture 'tex':
            Img = reshape(S(tt,:),[DIM(1),DIM(2)]);
            tex = Screen('MakeTexture',win,Img);
            Screen('DrawTexture',win, tex,[],destRect,[],0); % 0 is nearest neighbor
            % 1 is bilinear filter
            vbl = Screen('Flip',win);usb.strobeEventWord(1);
            vbl = Screen('Flip',win,vbl-ifi/2+flipInterval);
            usb.strobeEventWord(stimulusNumber(tt));
            vbl = Screen('Flip',win,vbl-ifi/2+WaitTimes(tt));
            Screen('Close',tex);
            tt = tt+1;
        end
    end
    if mod(tt,round(numStimuli/3)) == 0
        Screen('Flip',win);usb.strobeEventWord(0);
        WaitSecs(30);
    end
end

Screen('Flip',win);usb.strobeEventWord(0);
usb.strobeEventWord(endEXP);
WaitSecs(1);
usb.stopRecording;
% Close window
Screen('CloseAll');
Priority(0);

Date = datetime('today','Format','yyyy-MM-dd');
Date = char(Date); Date = strrep(Date,'-','');Date = str2double(Date);
filename = sprintf('RFPlaybackStim%s%d_%d.mat',NoiseType,Date,AnimalName);
save(filename,'S','numStimuli','flipInterval','numNoise',...
    'DistToScreen','screenPix_to_effPix','minPix','NoiseType',...
    'conv_factor','WaitTimes','beta','DIM','maxPix','stimulusNumber');
end

function gammaTable = makeGrayscaleGammaTable(gamma,blackSetPoint,whiteSetPoint)
% Generates a 256x3 gamma lookup table suitable for use with the
% psychtoolbox Screen('LoadNormalizedGammaTable',win,gammaTable) command
% 
% gammaTable = makeGrayscaleGammaTable(gamma,blackSetPoint,whiteSetPoint)
%
%   gamma defines the level of gamma correction (1.8 or 2.2 common)
%   blackSetPoint should be the highest value that results in a non-unique
%   luminance value on the monitor being used (sometimes values 0,1,2, all
%   produce the same black pixel value; set to zero if this is not a
%   concern)
%   whiteSetPoint should be the lowest value that returns a non-unique
%   luminance value (deal with any saturation at the high end)
% 
%   Both black and white set points should be defined on a 0:255 scale

gamma = max([gamma 1e-4]); % handle zero gamma case
gammaVals = linspace(blackSetPoint/255,whiteSetPoint/255,256).^(1./gamma);
gammaTable = repmat(gammaVals(:),1,3);
end