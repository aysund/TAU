function p = ad_endoTAUparams_v4(sessionType, temporalDistType)

% Written by Rachel Denison
% Updated by Antonio Fernandez - 2016 :)
% Updates by Aysun Duyar - 2020, 2022

% mkdir('staircase')
% addpath(genpath('staircase'))

if nargin<1
    sessionType = 'practice';
end

p.testingLocation = 'CarrascoR1'; % 'CarrascoL1','laptop','desk'

switch p.testingLocation
    case {'laptop','desk'}
        p.refRate = 1/60;
        p.screenSize = [9 13]; % (in)
        p.screenRes = [900 1440];
        p.viewDist = 21; % (in)
        p.eyeTracking = 1;
    case 'CarrascoR1'
        p.refRate = 1/100;
        p.screenSize = [40 30];
        p.screenRes = [1280 960];
        p.viewDist = 57;
        p.eyeTracking = 1;
    otherwise
        error('Testing location not found in temporalAttentionExoParams.m ')
end

switch sessionType
    case 'practice'
        % [1=valid -1=invalid 2=neutral(both) 0=neutral(none)]
        p.cueValidity = [1 3];
        p.cueValidityFactor = [1 3]; % eg. [1 1 2 3] is 50% valid, 25% invalid, 25% neutral (1=valid, 2=invalid, 3= neutral. repeat for the proportion)
        
        p.staircase = 0;
        p.isPractice = 1;        
        
        currStairt1 = 3; % center target version
        currStairt2 = 3;
        
        p.nBlocks = 1;
        p.eyeTracking = 0;
        p.nReps = 2; % should be 2
        
    case 'thresholding'
        p.cueValidity = 3;
        p.cueValidityFactor = 3;
        
        p.staircase = 1;        
        
        currStairt1 = 3.5; % center target version
        currStairt2 = 3.5;
        
        p.isPractice = 0;
        p.nReps = 4;
        p.nBlocks = 2;

    case 'experimental'
        p.cueValidity = [1 3];
        p.cueValidityFactor = [1 3]; % eg. [1 1 2 3] is 50% valid, 25% invalid, 25% neutral (1=valid, 2=invalid, 3= neutral. repeat for the proportion)
        
        currStairt1 = [];
        while isempty(currStairt1)
            currStairt1 = input('\n Enter stairt1 (value) obtained from the thresholding session: \n');
            currStairt2 = input('\n Enter stairt2 (value) obtained from the thresholding session: \n');
        end

        p.staircase = 0;
        p.isPractice = 0;

        p.nBlocks = 10;
        p.nReps = 2; % repeat per block
        
        if strcmp(temporalDistType, 'certain')
            p.nBlocks = 8;
        end
end
if p.isPractice
    fprintf(['!!!!!!! THIS IS A PRACTICE SESSION !!!!!!! \n PRACTICE PRACTICE PRACTICE PRACTICE \n',...
        'PRACTICE PRACTICE PRACTICE PRACTICE PRACTICE \n', ...
        'PRACTICE PRACTICE PRACTICE PRACTICE PRACTICE PRACTICE \n\n']);
end

p.keyNames = {'1!','2@','p','space', 'n', 'y'};
p.keyCodes = KbName(p.keyNames);
p.backgroundColor = 0.5;
p.fixationColor = 0.75; % aysun changed this to make fixation dimmer
p.goCueColor = 0.25; % aysun changed this because she changed the fixation color

p.font = 'Verdana';
p.fontSize = 20; % Rachel's was 24
p.showPlaceholders = 1;
p.phLineWidth = 2; % (pixels)
p.eyeRad = 1.5; % allowed fixation radius (degr22122212221222ppees)

% Condition
p.targetContrasts = 1; % [.64 1]; aysun changed
p.respInterval = [1 2]; % [1 = early 2 = late]
p.nTimingFactor   = 7; if strcmp(temporalDistType,'uniform'),p.nTimingFactor = 6; end
p.nTrialsPerBlock = p.nTimingFactor*numel(p.cueValidity)*p.nReps*numel(p.respInterval);

if (~strcmp(temporalDistType,'uniform') && p.nTrialsPerBlock ~=56) || (strcmp(temporalDistType,'uniform') && p.nTrialsPerBlock ~=48)  
    fprintf('CHECK TRIAL - BLOCK STRUCTURE!!');
    return
end

% Timing
p.readyT1SOAs = [1275 1400 1525]/1000;
p.T1T2SOA = 250/1000;
p.toneDur = 0.2;
p.targetDur = 0.03; % Antonio and Rachel has 30 ms
% p.readyRespCueSOA = 2.3; % max T2 onset -resp cue is 500 ms
% p.respGoSOA = p.T1respCueSOA + .5; % 
p.T1respcueSOA = 750/1000;
p.goCueON = 1;
p.respGoSOA = .5;
% p.respGoSOA = 1.0;
p.feedbackDur = 0.5; % inter-trial interval (also, the duration of the feedback symbol)
p.eyeSlack = 0.12; % cushion between last fixation check and next stimulus presentation

%%%%%%% Generate Trials
trialStructure = repmat(fullfact([ numel(p.cueValidity), numel(p.respInterval), 1, p.nTimingFactor]),p.nReps,1);

trials_1block = nan(size(trialStructure,1),4);
% col1: precue
trials_1block(trialStructure(:,1)==1, 1) = 3; % neutral
trials_1block(trialStructure(:,1)==2 & trialStructure(:,2)==1, 1) = 1; % T1 valid
trials_1block(trialStructure(:,1)==2 & trialStructure(:,2)==2, 1) = 2; % T2 valid
%col2: resp cue
trials_1block(:,2) = trialStructure(:,2);
% col3: t1 onset
% trials_1block(:,3) = p.readyT1SOA;
% col4: t2 onset
trials_1block(:,4) = p.T1T2SOA;

switch temporalDistType
    case 'uniform'
        p.timingDistr = [16, 16, 16]; % number of trials
        
        earlyInd = trialStructure(:,4)==1 | trialStructure(:,4)==2;
        midInd   = trialStructure(:,4)==3 | trialStructure(:,4)==4;
        lateInd  = trialStructure(:,4)==5 | trialStructure(:,4)==6;
        
    case 'wide'
        p.timingDistr = [16, 24, 16]; % number of trials
        
        earlyInd = trialStructure(:,4)==1 | trialStructure(:,4)==2;
        midInd   = trialStructure(:,4)==3 | trialStructure(:,4)==4 | trialStructure(:,4)==5;
        lateInd  = trialStructure(:,4)==6 | trialStructure(:,4)==7;
        
    case 'narrow'
        p.timingDistr = [4, 48, 4]; % number of trials
        
        earlyInd = trialStructure(:,4)==1 & find(trialStructure(:,4))<28;
        midInd   = trialStructure(:,4)~=1;
        lateInd  = trialStructure(:,4)==1 & find(trialStructure(:,4))>28;
        
    case 'certain'
        p.timingDistr = [0, 56, 0]; % number of trials
        
        earlyInd = false(size(trialStructure(:,4)));
        midInd   = true(size(trialStructure(:,4)));
        lateInd  = false(size(trialStructure(:,4)));
end

trials_1block(earlyInd,3) = p.readyT1SOAs(1);
trials_1block(midInd,3)   = p.readyT1SOAs(2);
trials_1block(lateInd,3)  = p.readyT1SOAs(3);

% introduce uncertainty for T1 (75 ms +-mid uncertainty point
% (such that it's 150 ms uncertainty in the early and 150 ms uncertainty in
% the late pts than the "expected" point))
earlyJitter = [(-1).*randi(75,sum(earlyInd)/2,1,p.nBlocks); randi(75,sum(earlyInd)/2,1,p.nBlocks)]./1000;
lateJitter = [(-1).*randi(75,sum(earlyInd)/2,1,p.nBlocks); randi(75,sum(earlyInd)/2,1,p.nBlocks)]./1000;

myRandomizer1 = zeros(sum(earlyInd),p.nBlocks); myRandomizer2 = zeros(sum(lateInd),p.nBlocks);
myRandomizer1(:,1) = randperm(sum(earlyInd))';  myRandomizer2(:,1) = randperm(sum(lateInd))';
for ii = 2: p.nBlocks
    myRandomizer1(:,ii) = randperm(sum(earlyInd))';
    myRandomizer2(:,ii) = randperm(sum(lateInd))';
end
myRandomizer1 = reshape(myRandomizer1,sum(earlyInd),1,p.nBlocks);
myRandomizer2 = reshape(myRandomizer2,sum(lateInd),1,p.nBlocks);
earlyJitter = earlyJitter(myRandomizer1);
lateJitter  = lateJitter(myRandomizer2);

p.randTimings = nan(size(trialStructure,1),1,p.nBlocks);
p.randTimings(midInd,:,:)   = 0; % mid timing, no jitter
p.randTimings(earlyInd,:,:) = earlyJitter; % earlier than mid
p.randTimings(lateInd,:,:)  = lateJitter; % later than mid

% set target axis and cw/ccw
targOris = repmat(fullfact([2,2,2,2]), 3 ,1);
if ~strcmp(temporalDistType,'uniform'), targOris = [targOris; targOris(1:8,:)]; end
targOris(targOris(:,1)==1,1)=0; targOris(targOris(:,2)==1,2)=0;
targOris(targOris(:,1)==2,1)=90; targOris(targOris(:,2)==2,2)=90;
targOris(targOris(:,3)==1,3)=-1; targOris(targOris(:,4)==1,4)=-1; 
targOris(targOris(:,3)==2,3)=1; targOris(targOris(:,4)==2,4)=1;  

trials_allblocks = cat(2, repmat(trials_1block,[1,1,p.nBlocks]), repmat(targOris,[1,1,p.nBlocks]));
if strcmp(sessionType, 'thresholding')
    trials_allblocks(:,3,:) = p.readyT1SOAs(2);
else
    trials_allblocks(:,3,:) = trials_allblocks(:,3,:) + p.randTimings;
end


% trials_allblocks(:,[5 6 7 8],:) = Shuffle(trials_allblocks(:,[5 6 7 8],:),[2]);
for ii = 1:p.nBlocks
    randomizer = randperm(p.nTrialsPerBlock);
    trials_allblocks(:,[5 6 7 8],ii) = trials_allblocks(randomizer,[5 6 7 8],ii);
    clear randomizer
end
p.allTrials = trials_allblocks;

% Images
p.imSize = [4 4]; % this is the size of the image container that holds the stim
% p.targetSize = 0.25; % 0.5 sigma of gaussian / 1.5 side length of T/L / 1.5 width of triangle
p.imPos = [0 0]; %jun21 - central target
p.targetSize = .3; %AF
p.spatialFrequency = 4; % jun21 - central target

p.rotDirs  = [-1 1];  % -1: CCW, 1: CW
p.targetAx = [0, 90]; % 0: horizontal, 90: vertical
% p.targetAx = [45, 135]; % 0: horizontal, 90: vertical

p.plHoldSize = 0.2;
p.plHoldDist = 0.1;
imdist = 1.2; %jun21
p.plHoldPos = [-imdist-p.plHoldDist, -imdist-p.plHoldDist; imdist+p.plHoldDist, -imdist-p.plHoldDist; imdist+p.plHoldDist, imdist+p.plHoldDist; -imdist-p.plHoldDist, imdist+p.plHoldDist]'; % square shape
p.plHoldColor = [1 1 1];
p.plHoldLineWidth = 2; % if circle around target

p.fixOut = 9;
p.fixIn = 7;
p.feedbackSize = p.fixOut;

p.targetOrientation = [-1 1]; %AF

if p.staircase %AF
    stairParams.whichStair = 1; %AF
    %     stairParams.alphaRange = [0.5 1 1.5 2 3 4 5 6 7 8]; %AF % steps of the staircase
%     stairParams.alphaRange = [0.5 1 1.5 2 2.5 3 4 5 6 7]; %AF % steps of the staircase
    stairParams.alphaRange = .25:.5:3.5; % target center version
%     switch stimLoc
%         case 'fovea'
%             stairParams.alphaRange = .25:.5:3.5; % target center version
%         case 'lowerRight'
%             stairParams.alphaRange = [.25:.5:4.5 4.5]; % lower right version
%     end
    
    stairParams.fitBeta = 2; %AF params for the psychoetric fn
    stairParams.fitLambda = 0.01; %AF
    stairParams.fitGamma = 0.5; %AF
    stairParams.perfLevel =0.75; %AF
    
    % T1
    stairt1Params = stairParams; %AF
    stairt1Params.useMyPrior = []; %AF
    p.stairst1 = usePalamedesStaircase(stairt1Params); %AF
    p.stairst1.xCurrent = currStairt1;
    
    % T2
    stairt2Params = stairParams; %AF
    stairt2Params.useMyPrior = []; %AF
    p.stairst2 = usePalamedesStaircase(stairt2Params); %AF
    p.stairst2.xCurrent = currStairt2;
    
    fprintf('\nStaircase is ON\n'); %AF
else
    p.stairst1.xCurrent = currStairt1;
    p.stairst2.xCurrent = currStairt2;
end

% staircase of T2 trial to trial adjustment flag
% We care about mid point performance
p.isT1mid = 0; 

% Sounds
p.Fs = 44100;
p.cueFreqs = [1000, 400];
for iTone = 1:numel(p.cueFreqs)
    tone = MakeBeep(p.cueFreqs(iTone), p.toneDur, p.Fs);
    p.cueTones(iTone,:) = applyEnvelope(tone, p.Fs);
end

tb = [0 .1 .2];
fb = [50 200 400];
tFinal = .2;
t=linspace(0,tFinal,p.Fs*tFinal);
f=interp1(tb,fb,t);
y=sin(f*2*pi.*t); y = [0,y];

p.cueTones= [ p.cueTones; y];
 

 
% 10^0.5 for every 10dB