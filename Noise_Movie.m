function [] = Noise_Movie(AnimalName)
%Noise_Movie.m
%   Display a movie of noise stimuli to infer the receptive fields
%    of neurons using reverse correlation or a GLM point process model
%    See Smyth et al. 2003 Receptive Field Organization for evaluation of 
%     spatial (not spatiotemporal) receptive fields using natural images
%
%
%INPUT: AnimalName - unique identifier for the animal as a number, e.g.
%            12345
%
%       see file MovieNoiseVars.mat for more changeable presets
%         movies previously made and saved as 10Sec_PinkNoiseMovie1.mat
%         etc.
%
%OUTPUT: file named 'NoiseStimNoiseTypeDate_AnimalName.mat' , e.g. 
%          NoiseStimpink20160718_12345.mat
%
%        S - matrix sized numStimuli-by-numPixels that represent each of
%          the stimuli presented, try 
%          image = reshape(S(1,:),[width,height]); to view one of the white
%          noise stimuli
%        effectivePixels - effective width and height of the display window 
%          in pixels, first chosen as the minimum of the width and height 
%          of the current display ... the display is forced to be square
%          and the width in effectivePixels is 1/4 the width in true screen
%          pixels
%        DistToScreen - 25 cm for now
%
% Created: 2017/03/24, 24 Cummington, Boston
%  Byron Price
% Updated: 2017/06/22
% By: Byron Price

cd('~/CloudStation/ByronExp/NoiseRetino')
load('MovieNoiseVars.mat');

% for motion-contingent display / interaction with recording computer
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

fprintf('\nEstimated time is %3.2f minutes.',(numMoviesToDisplay*movieTime_Seconds*2+3*30)/60);
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

% screen size in millimeters and a conversion factor to get from mm to pixels
[w_mm,h_mm] = Screen('DisplaySize',screenid);
conv_factor = (w_mm/w_pixels+h_mm/h_pixels)/2;

wLow = round((w_pixels-maxPix)/2);
wHigh = round(w_pixels-wLow);
hLow = round((h_pixels-minPix)/2);
hHigh = round(h_pixels-hLow);
destRect = [wLow hLow wHigh hHigh];

Priority(9);
% Retrieve monitor refresh duration
ifi = Screen('GetFlipInterval', win);deltaIFI = ifi/2;
% flipIntervals = flipInterval+exprnd(0.1,[numStimuli,1]);
% WaitTimes = WaitTime+exprnd(0.25,[numStimuli,1]);
flipInterval = 1/movie_FrameRate-deltaIFI;

movieNums = random('Discrete Uniform',250,[numMoviesToDisplay,1]);

usb.startRecording;WaitSecs(1);usb.strobeEventWord(0);
WaitSecs(30);
usb.strobeEventWord(startEXP);WaitSecs(1);
vbl = Screen('Flip',win);
count = 1;
while count<=numMoviesToDisplay
    if tcpipClient.BytesAvailable > 0
        data = fread(tcpipClient,tcpipClient.BytesAvailable/8,'double');
        if sum(data) > 0
            WaitSecs(5);
            fprintf('Motion Pause ...\n');
        else
            index = movieNums(count);
            load(sprintf('12Sec_PinkNoiseMovie%d.mat',index),'S','numStimuli');
            %     S(S<60) = 0;S(S>=60 & S<196) = 127;S(S>=196) = 255;
            tt = 1;
            while tt <= numStimuli
                %     Img = uint8(kron(single(S(:,:,tt)),ones(screenPix_to_effPix)));
                % Img = S(:,:,tt);Img = Img(newInds);
                % Convert it to a texture 'tex':
                tex = Screen('MakeTexture',win,S(:,:,tt));
                Screen('DrawTexture',win, tex,[],destRect,[],0); % 0 is nearest neighbor
                % 1 is bilinear filter
                vbl = Screen('Flip',win,vbl+flipInterval);usb.strobeEventWord(index);
                Screen('Close',tex);
                tt = tt+1;
            end
            vbl = Screen('Flip',win);usb.strobeEventWord(0);
            if mod(count,30) == 0
                vbl = Screen('Flip',win,vbl+30);
            end
            waitTime = 10+normrnd(0,1);
            vbl = Screen('Flip', win,vbl+waitTime);
            count = count+1;
        end
    end
end
Screen('Flip',win);usb.strobeEventWord(0);
usb.strobeEventWord(endEXP);
WaitSecs(1);

usb.stopRecording;
% Close window
Screen('CloseAll');
Priority(0);

movieType = '12Sec_PinkNoise';

Date = datetime('today','Format','yyyy-MM-dd');
Date = char(Date); Date = strrep(Date,'-','');Date = str2double(Date);
filename = sprintf('NoiseMovieStim%s%d_%d.mat','pink',Date,AnimalName);
save(filename,'numMoviesToDisplay','movie_FrameRate','movieTime_Seconds',...
    'DistToScreen','screenPix_to_effPix','minPix','maxPix',...
    'conv_factor','movieNums','destRect','movieType');
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
