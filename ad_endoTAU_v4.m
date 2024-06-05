function ad_endoTAU_v4()
close all
clear all

addpath(genpath('/Users/purplab/Desktop/Aysun/endoTAU/'));

% Written by Rachel Denison
% Updates by Aysun Duyar
     
Screen('CloseAll')
%% Setup
% block = 1; % this is only to still display a break message if the last trial in a block is skipped due to fixation break
instruct = input('\n Do you want instructions? (0 or 1):\n');

c = clock;  
expt.date = c;
randomizeRange = round(sum(c(3:end)));
rng(randomizeRange);

testNum = 999;
subjectID = [];
while isempty(subjectID)
    Prompt        = {'Subject ID (testNum = 999):', 'Session Type (99-practice 0-thresholding 1-experimental):','Session Number:', 'Expectation (0-certain, 1-narrow, 2-wide, 3-uniform)'};
    Answer        = inputdlg(Prompt,'Info',1);
    subjectID     = str2double(Answer{1});
    sesT          = str2double(Answer{2});
    sessionNum    = str2double(Answer{3});
    expectationType = str2double(Answer{4});
end

switch sesT
    case 99
        sessionType = 'practice';
    case 0
        sessionType = 'thresholding';
    case 1
        sessionType = 'experimental';
    otherwise
        error('Wrong session type input!')
end

switch expectationType
    case 0
        temporalDistType = 'certain';
    case 1
        temporalDistType = 'narrow';
    case 2
        temporalDistType = 'wide';
    case 3
        temporalDistType = 'uniform';
    otherwise
        error('Wrong expectation input!');
end

p = ad_endoTAUparams_v4(sessionType, temporalDistType);

sum(sum(sum(p.allTrials)))

if isequal(subjectID,testNum)
    p.eyeTracking = 0;
end

if subjectID == 11 %RC
    switch sessionType
        case 'practice'
            currStairt1 = 1; % center target version
            currStairt2 = 1;
        case 'thresholding'
            currStairt1 = 2.25; % center target version
            currStairt2 = 2.25;           
    end
    stairParams.alphaRange = linspace(.25,2.25,7); % target center version
end

slack = p.refRate/2;
rad = round(ang2pix(p.eyeRad, p.screenSize(1), p.screenRes(1), p.viewDist, 'central'));

% Running on PTB-3? Abort otherwise.
AssertOpenGL;

%% Display key settings to the experimenter

fprintf('\nExperiment settings:\n')
fprintf('tilts (T1, T2) = [%1.1f, %1.1f]\n', p.stairst1.xCurrent,p.stairst2.xCurrent)
fprintf('expectation = %s \n', temporalDistType)
fprintf('ID = %.f \n\n',subjectID)

ok = input('\n Settings ok? [n if not]','s');
if strcmp(ok,'n')
    error('Ok, check parameters')
    return
end

%% Check for existing data file
cd '/users/purplab/Desktop/Aysun/endoTAU/data'
eyeDataDir = 'eyedata';
switch sessionType
    case 'experimental'
        fname = sprintf('tau_s%.f_%d_v3.txt', subjectID,sessionNum);            % Name of the .txt file, according to subject number
        eyeFile = sprintf('eTAU');
        matFileName = sprintf('ex_s%.0f_%d_v3.mat', subjectID,sessionNum);
    case 'practice'
        fname = sprintf('PRtau_s%.f_v3.txt', subjectID);
        eyeFile = sprintf('pTAU');
        matFileName = sprintf('pr_s%.0f_%d_v3.mat', subjectID, sessionNum);
    case 'thresholding'
        fname = sprintf('THtau_s%.f_%d_v3.txt',subjectID, sessionNum);
        eyeFile = sprintf('tTAU');
        matFileName = sprintf('th_s%.0f_%d_v3.mat', subjectID, sessionNum);
end

if ~exist(fname, 'file')
    fp = fopen(fname, 'w');                                                 % Open the file
    fprintf(fp, 'random point = %.f \r\n', randomizeRange');
    fprintf(fp, ['trialNum\t blockNum\t stairT1\t stairT2\t respInterval\t'...
        'cueValidity\t cuedInterval\t t1or\t t2or\t t1ax\t t2ax\t t1phase\t t2phase\t',...
        ' RT\t RespKey\t Response\t correct\t currITI\t readyT1SOA\t readyT2SOA \n']);   % Prints the headers
else
    warn = sprintf('WARNING: File Exists. Are you sure about the subject ID = %d?', subjectID);
    ButtonName = questdlg(warn, 'Log File Exists','Yes','No','No');
    switch ButtonName
        case 'No'
            return;
    end
    fp = fopen(fname, 'a+');
end

%% Eye data i/o
% Check to see if this eye file already exists
existingEyeFile = dir(sprintf('%s/%s.edf', eyeDataDir, eyeFile));
if ~isempty(existingEyeFile) && p.eyeTracking
    %     error('eye file already exists! please choose another name.')
    fprintf('EYE FILE ALREADY EXISTS!! GENERATNING A NEW ONE...');
    eyeFile = sprintf('%s0', eyeFile);
end

%% Keyboard
% Find keyboard device number
devNum = PsychHID('devices',-1);
% devNum = 3;

if isempty(devNum)
    error('Could not find Keypad!')
end

%% Sound
% Perform basic initialization of the sound driver
InitializePsychSound(1); % 1 for precise timing

% Open audio device for low-latency output
reqlatencyclass = 2; % Level 2 means: Take full control over the audio device, even if this causes other sound applications to fail or shutdown.
pahandle = PsychPortAudio('Open', [], [], reqlatencyclass, p.Fs, 1); % 1 = single-channel
PsychPortAudio('Volume', pahandle, .6); % 1 denotes 100% volume.

%% Screen
% Set up screen
screenNumber = max(Screen('Screens'));

% Check screen resolution and refresh rate
scr = Screen('Resolution', screenNumber);
if ~all([scr.width scr.height scr.hz] == [p.screenRes round(1/p.refRate)]) && ~isequal(subjectID,testNum)
    error('Screen resolution and/or refresh rate has not been set correctly by the experimenter!')
end

% Set resolution and refresh rate
if strcmp(p.testingLocation, 'CarrascoR1')
    Screen('Resolution',screenNumber, p.screenRes(1), p.screenRes(2), round(1/p.refRate));
end

% Set up window
[window, rect] = Screen('OpenWindow', screenNumber);
% [window rect] = Screen('OpenWindow', screenNumber, [], [0 0 800 600]);
white = WhiteIndex(window);  % Retrieves the CLUT color code for white.
black = BlackIndex(window);
[cx, cy] = RectCenter(rect);
center = [cx, cy];
Screen('TextSize', window, p.fontSize);
Screen('TextColor', window, white);
Screen('TextFont', window, p.font);
p.fixationITIColor = black;

% Check screen size
[sw, sh] = Screen('WindowSize', window); % height and width of screen (px)
if ~all([sw sh] == p.screenRes) && ~isequal(subjectID,testNum)
    error('Screen resolution is different from requested!')
end

% Check refresh rate
flipInterval = Screen('GetFlipInterval', window); % frame duration (s)
if abs(flipInterval - p.refRate) > 0.001 && ~isequal(subjectID,testNum)
    error('Refresh rate is different from requested!')
end

% Check font
if ~strcmp(p.font, Screen('TextFont', window))
    error('Font was not set to requested: %s', p.font)
end

% Load calibration file
switch p.testingLocation
    case 'CarrascoR1'
        %load('/users/purplab/Desktop/Aysun/temporal-attention-master/Displays/GammaTable_R1.mat', GammaTable);
        load('/users/purplab/Desktop/MonitorCalibrationScripts/GammaTables/GammaTable_R1_IBM_20220814.mat');
        Screen('LoadNormalizedGammaTable', window,GammaTable);
        % check gamma table
        gammatable = Screen('ReadNormalizedGammaTable', window);
        if nnz(abs(gammatable-GammaTable)>0.0001)
            error('Gamma table not loaded correctly! Perhaps set screen res and retry.')
        end
    otherwise
        fprintf('\nNot loading gamma table ...\n')
end

%% Stimulus Params
% Calculate stimulus dimensions (px) and position
imPos = round(ang2pix(p.imPos, p.screenSize(1), p.screenRes(1), p.viewDist, 'radial')); % from screen center
targetSize = round(ang2pix(p.targetSize, p.screenSize(1), p.screenRes(1), p.viewDist, 'central'));
plHoldSize = round(ang2pix(p.plHoldSize, p.screenSize(1), p.screenRes(1), p.viewDist, 'central'));
plHoldPos = round(ang2pix(p.plHoldPos, p.screenSize(1), p.screenRes(1), p.viewDist, 'central')); % from target center

pixelsPerDegree = round(ang2pix(1, p.screenSize(1), p.screenRes(1), p.viewDist, 'central'));

% Make the rects for placing the images
imSize = pixelsPerDegree*p.imSize(1)+1;
imRect = CenterRectOnPoint([0 0 imSize imSize], cx+imPos(1), cy+imPos(2));

%% Generate trials
% Construct trials matrix
% a "state" is one of two states the target can take, eg left vs. right
% orientation, or low vs. high SF
 trials_headers = {'trialNum','blockNum','stairLevelT1H','stairLevelT1V', 'stairLevelT2H', 'stairLevelT2V',...
    'respInterval','cueValidity','cuedInterval',...
    't1or','t2or','t1ax','t2ax','t1phase','t2phase',...
    'rt','responseKey','response','correct', 'ITI', 'readyT1SOA', 'readyT2SOA'};

% make sure column indices match trials headers
% trialNumIdx = strcmp(trials_headers,'trialNum');
% blockNumIdx = strcmp(trials_headers,'blockNum');
% stairLevelT1Idx = strcmp(trials_headers,'stairLevelT1');
% stairLevelT2Idx = strcmp(trials_headers,'stairLevelT1');
% respIntervalIdx = strcmp(trials_headers,'respInterval');
% cueValidityIdx = strcmp(trials_headers,'cueValidity');
% cuedIntervalIdx = strcmp(trials_headers,'cuedInterval');
% t1orIdx = strcmp(trials_headers,'t1or');
% t2orIdx = strcmp(trials_headers,'t2or');
% t1axIdx = strcmp(trials_headers,'t1ax');
% t2axIdx = strcmp(trials_headers,'t2ax');
% t1phIdx = strcmp(trials_headers,'t1phase');
% t2phIdx = strcmp(trials_headers,'t2phase');
% rtIdx = strcmp(trials_headers,'rt');
% responseKeyIdx = strcmp(trials_headers,'responseKey');
% responseIdx = strcmp(trials_headers,'response');
% correctIdx = strcmp(trials_headers,'correct');
% ITIIdx = strcmp(trials_headers,'ITI');

% randomize the order of trials, for each block.
trialOrders = zeros(p.nBlocks,size(p.allTrials,1));
for tt = 1:p.nBlocks
    trialOrders(tt,:) = randperm(size(p.allTrials,1));
end
nTrials = numel(trialOrders);

% show trials and blocksk
fprintf('\n%s\n\n%d trials, %1.2f blocks, %f \n\n', datestr(now), nTrials, p.nBlocks);

%% Present Instructions
if instruct
    sampTarget = buildColorGrating(pixelsPerDegree, p.imSize, ...
        p.spatialFrequency, -1, 90, 1, 0, 'bw',[],[],p.backgroundColor);
    % mask with an aperture (eg. 2d gaussian)
    sampStim = maskWithGaussian(sampTarget, size(sampTarget,1), targetSize);
    % Make textures
    sampTex = Screen('MakeTexture', window, sampStim*white);
%     ad_temporalExoInstructions(sessionType, window, devNum, white, p, imRect, sampTex, pahandle);
end

%% Eyetracker
if p.eyeTracking
    % Initialize eye tracker
    [el, exitFlag] = rd_eyeLink('eyestart', window, eyeFile);
    if exitFlag
        return
    end
    
    % Write subject ID into the edf file
    Eyelink('message', 'BEGIN DESCRIPTIONS');
    Eyelink('message', sprintf('Subject code: %.f', subjectID));
    Eyelink('message', 'END DESCRIPTIONS');
    
    % No sounds indicating success of calibration
    %     el.targetbeep = false;
    %     el.calibration_failed_beep = [0 0 0];
    %     el.calibration_success_beep = [0 0 0];
    el.drift_correction_target_beep = [0 0 0];
    el.drift_correction_failed_beep = [0 0 0];
    el.drift_correction_success_beep = [0 0 0];
    
    % Accept input from all keyboards
    el.devicenumber = -1; % see KbCheck for details of this value
    
    % Update with custom settings
    EyelinkUpdateDefaults(el);
    
    % Calibrate eye tracker
    [~, exitFlag] = rd_eyeLink('calibrate', window, el);
    if exitFlag
        return
    end
else
    fprintf('EYETRACKING IS OOOOOOOFFFFFFFFFF!!!\n EYETRACKING IS OOOOOOOFFFFFFFFFF!!!\nEYETRACKING IS OOOOOOOFFFFFFFFFF!!!\n')
end

%% Present trials
numberDur = 1; % duration of the numbers during countdown

% Show fixation and wait for a button press
Screen('FillRect', window, white*p.backgroundColor);
DrawFormattedText(window, 'Press Space to start the experiment', 'center', 'center', [1 1 1]*white);
Screen('Flip', window);
[getsc, ~, ~] = KbWait(devNum);

% Trials
eyeSkip = zeros(size(p.allTrials,1),1); % trials skipped due to an eye movement, same size as trials matrix

% workspaceFile = [];
% % Option to load a previous run from a saved workspace file (the TEMP.mat file)
% % note this will overwrite most of the settings generated above
% if ~isempty(workspaceFile)
%     eyeFile0 = eyeFile;
%     load(workspaceFile)
%     fprintf('\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n !!!!!!!!!!!!!!! LOADING FROM WORKSPACE FILE:\n%s\n\n', workspaceFile)
%     subjectID(end+1) = 'W';
%     eyeFile = eyeFile0; % use a new eye file
%     p.forwardMask = 0; % ad only
%     block = floor(iTrial/p.nTrialsPerBlock); % ad only
% end

Screen('FillRect', window, white*p.backgroundColor);
timeFeedback = Screen('Flip', window, getsc + 0.5);
stairValues = zeros(nTrials,2);
% % stairValues = zeros(nTrials,4);
iTrialAll   = 1;
timing.startTime = GetSecs;

overallAcc = [0 0];
neutrAccPerc_titrate = [0 0];
timeGoCue=0;

% Present trials
for bb = 1:p.nBlocks
    trialCounter = 1;  % reset the counter for each block
    blockAcc = [0, 0]; % reset every block
    neutrAcc = [0, 0];
    neutrnumMat = [0, 0];
    
    clear trials trialOrder
    trialOrder = trialOrders(bb, :);
    trials = p.allTrials(:,:,bb);
    nTrialsPerBlock = size(trials,1); % reset the number of trials within block
    
    while trialCounter <= nTrialsPerBlock
        % Initialize for eye tracking trial breaks
        if trialCounter>1
            eyeSkip(trialIdx) = stopThisTrial; % this is for the previous trial
        end
        stopThisTrial = 0;
        
        % Current trial number
        iTrial       = trialCounter;       % the trial we're on now
        trialIdx     = trialOrder(iTrial); % the current index into trials
        trialCounter = trialCounter+1;     % update the trial counter so that we will move onto the next trial, even if there is a fixation break
        
        % Get conditions for this trial
        cuedInterval = trials(trialIdx,1); % 1, 2, or 3 (neutral)
        respInterval = trials(trialIdx,2); % 1 or 2
        readyT1SOA   = trials(trialIdx,3); % randomized
        readyT2SOA   = readyT1SOA + trials(trialIdx,4); % T1+250 ms
        rots         = trials(trialIdx,[5, 6]); % rotations (H or V)
        
        % set whether T2 should be adjusted
        if readyT1SOA==p.readyT1SOAs(2)
            p.isT1mid = true;
        else
            p.isT1mid = false;
        end
        
        angles = [0, 0];
        % T1
        angles (1) = (trials(trialIdx,7)).*p.stairst1.xCurrent; %angles, based on the stair level
        % T2
        angles (2) = (trials(trialIdx,8)).*p.stairst2.xCurrent; %angles, based on the stair level
        tilts = rots + angles;
        
        % random phases
        ph = round(rand(2,1)*360);
        
        % select the random intertrial interval
        % currITI = .9 + (1.4-.9)*rand;
        currITI = 1.1 + (1.7-1.1)*rand;
        
        % Generate target textures
        % make big grating
        t1 = buildColorGrating(pixelsPerDegree, p.imSize, ...
            p.spatialFrequency, tilts(1), ph(1), p.targetContrasts, 0, 'bw',[],[],p.backgroundColor);
        t2 = buildColorGrating(pixelsPerDegree, p.imSize, ...
            p.spatialFrequency, tilts(2), ph(2), p.targetContrasts, 0, 'bw',[],[],p.backgroundColor);
        
        % mask with an aperture (eg. 2d gaussian)
        stim1 = maskWithGaussian(t1, size(t1,1), targetSize);
        stim2 = maskWithGaussian(t2, size(t2,1), targetSize);
        % Make textures
        tex1 = Screen('MakeTexture', window, stim1*white);
        tex2 = Screen('MakeTexture', window, stim2*white);
        
        % Select tones
        precueTone = p.cueTones(cuedInterval,:);
        respTone = p.cueTones(respInterval,:);
        
        % Present ITI - jittered
        Screen('DrawDots', window, plHoldPos, plHoldSize, black, center, 2); % placeholders
        putFixExo(window,[cx cy], p.fixOut, p.fixIn, white*p.fixationITIColor, p.backgroundColor*white)
        timeITI = Screen('Flip', window, timeFeedback+p.feedbackDur);
        
        % Present fixation
        Screen('DrawDots', window, plHoldPos, plHoldSize, black, center, 2); % placeholders
        putFixExo(window,[cx cy], p.fixOut, p.fixIn, white*p.fixationColor, p.backgroundColor*white)
        timeFix = Screen('Flip', window, timeITI + currITI);
        
        % Check fixation hold
        if p.eyeTracking
            driftCorrected = rd_eyeLink('trialstart', window, {el, iTrialAll, cx, cy, rad});
            if driftCorrected
                % restart trial
                Screen('DrawDots', window, plHoldPos, plHoldSize, black, center, 2); % placeholders
                putFixExo(window,[cx cy], p.fixOut, p.fixIn, white*p.fixationITIColor, p.backgroundColor*white)
                timeFix = Screen('Flip', window, timeFeedback+p.feedbackDur);
            end
        end
        
        % Present precue
        PsychPortAudio('FillBuffer', pahandle, precueTone);
        timePrecue = PsychPortAudio('Start', pahandle, [], timeITI + currITI, 1);
        
        if p.eyeTracking
            Eyelink('Message', 'EVENT_PRECUE');
        end
        
        timeStart = timePrecue; % deneme
        
        % Check for eye movements
        if p.eyeTracking
            while GetSecs < timeStart + readyT1SOA - p.eyeSlack && ~stopThisTrial
                WaitSecs(.01);
                fixation = rd_eyeLink('fixcheck', window, {cx, cy, rad});
                [stopThisTrial, trialOrder, nTrialsPerBlock] = fixationBreakTasks(...
                    fixation, window, white*p.backgroundColor, trialOrder, iTrial, nTrialsPerBlock);
            end
            fixCue(iTrialAll) = fixation;
            if stopThisTrial
                continue
            end
        end
                
        % Present images
        %%% T1 %%%
        Screen('DrawTexture', window, tex1, [], imRect);
        Screen('DrawDots', window, plHoldPos, plHoldSize, black, center, 2); % placeholders
        timeIm1 = Screen('Flip', window, timeStart + readyT1SOA - slack);
        if p.eyeTracking
            Eyelink('Message', 'EVENT_T1');
        end
        
        % blank
        Screen('FillRect', window, white*p.backgroundColor);
        % ad_drawPlaceholders(window, p.exoColor*white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders, [cx, cy], exoPos, exoSize);
        putFixExo(window,[cx cy], p.fixOut, p.fixIn, white*p.fixationColor, p.backgroundColor*white)
        Screen('DrawDots', window, plHoldPos, plHoldSize, black, center, 2); % placeholders
        timeBlank1 = Screen('Flip', window, timeIm1 + p.targetDur - slack);
        
        if p.eyeTracking
            Eyelink('Message', 'ISI1');
        end
        
        % Check for eye movements
        if p.eyeTracking
            while GetSecs < timeStart + readyT2SOA - p.eyeSlack && ~stopThisTrial
                WaitSecs(.01);
                fixation = rd_eyeLink('fixcheck', window, {cx, cy, rad});
                [stopThisTrial, trialOrder, nTrialsPerBlock] = fixationBreakTasks(...
                    fixation, window, white*p.backgroundColor, trialOrder, iTrial, nTrialsPerBlock);
            end
            fixT1(iTrialAll) = fixation;
            if stopThisTrial
                continue
            end
        end
        
        %%% T2 %%%
        Screen('DrawTexture', window, tex2, [], imRect);
        Screen('DrawDots', window, plHoldPos, plHoldSize, black, center, 2); % placeholders
        timeIm2 = Screen('Flip', window, timeStart + readyT2SOA - slack);
        
        if p.eyeTracking
            Eyelink('Message', 'EVENT_T2');
        end
        
        % blank
        Screen('FillRect', window, white*p.backgroundColor);
        putFixExo(window,[cx cy], p.fixOut, p.fixIn, white*p.fixationColor, p.backgroundColor*white)
        Screen('DrawDots', window, plHoldPos, plHoldSize, black, center, 2); % placeholders
        timeBlank2 = Screen('Flip', window, timeIm2 + p.targetDur - slack);
        
        if p.eyeTracking
            Eyelink('Message', 'EVENT_ISI2');
        end
        
        % Check for eye movements
        if p.eyeTracking
            while GetSecs < timeStart + readyT1SOA + p.T1respcueSOA - p.eyeSlack && ~stopThisTrial
                WaitSecs(.01);
                fixation = rd_eyeLink('fixcheck', window, {cx, cy, rad});
                [stopThisTrial, trialOrder, nTrialsPerBlock] = fixationBreakTasks(...
                    fixation, window, white*p.backgroundColor, trialOrder, iTrial, nTrialsPerBlock);
            end
            fixT2(iTrialAll) = fixation;
            if stopThisTrial
                continue
            end
        end
        
        % Present response cue
        PsychPortAudio('FillBuffer', pahandle, respTone);
        timeRespCue = PsychPortAudio('Start', pahandle, [], timeStart + readyT1SOA + p.T1respcueSOA, 1);
        if p.eyeTracking
            Eyelink('Message', 'EVENT_RESPCUE');
        end
        
        
        % Check for eye movements
        if p.eyeTracking
            while GetSecs < timeStart + readyT1SOA + p.T1respcueSOA && ~stopThisTrial
                WaitSecs(.01);
                fixation = rd_eyeLink('fixcheck', window, {cx, cy, rad});
                [stopThisTrial, trialOrder, nTrialsPerBlock] = fixationBreakTasks(...
                    fixation, window, white*p.backgroundColor, trialOrder, iTrial, nTrialsPerBlock);
            end
            fixRC(iTrialAll) = fixation;
            if stopThisTrial
                continue
            end
        end
        
        if p.goCueON
            % Present go cue (indicating you're allowed to make a response)
            Screen('FillRect', window, white*p.backgroundColor);
            % ad_drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders)
            Screen('DrawDots', window, plHoldPos, plHoldSize, black, center, 2); % placeholders
            putFixExo(window,[cx cy], p.fixOut+1.5, p.fixIn, white*p.goCueColor, p.backgroundColor*white)
            
            WaitSecs(p.respGoSOA - slack);
            timeGoCue = Screen('Flip', window, timeRespCue + p.respGoSOA - slack);
            
            if p.eyeTracking
                Eyelink('Message', 'EVENT_GOSIGNAL');
            end
        end
        
        % Collect response
        responseKey = [];
        while isempty(responseKey) % record wrong key as missed trial
            [secs, keyCode] = KbWait(devNum);
            rt = secs - timeRespCue;
            
            if numel(find(keyCode))>1 % if more than one key was pressed simultaneously
                responseKey = [];
            else
                responseKey = find(p.keyCodes(1:3)==find(keyCode));
                if responseKey == 3 % if p pressed
                    ShowCursor
                    save(['/users/purplab/Desktop/Aysun/endoTAU/data/',matFileName]) % saves the workspace on each trial
                    % Save eye data and shut down the eye tracker
                    if p.eyeTracking
                        rd_eyeLink('eyestop', window, {eyeFile, eyeDataDir});
                        
                        % rename eye file
                        eyeFileFull = sprintf('%s/%.fs%.f_%.f_v3.edf', eyeDataDir, sesT, subjectID,sessionNum);
                        copyfile(sprintf('%s/%s.edf', eyeDataDir, eyeFile), eyeFileFull)
                    end
                    
                    % Clean up
                    PsychPortAudio('Stop', pahandle);
                    PsychPortAudio('Close', pahandle);
                    Screen('CloseAll')
                    ShowCursor
                    break
                end
            end
        end
        response = responseKey;
        
        if p.eyeTracking
            Eyelink('Message', 'TRIAL_END');
        end
        
        % Feedback
        correctResp = trials(trialIdx,respInterval+6); % 1 or -1
        correctResp = find(correctResp==[-1, 1]); % cw or ccsw (1 or 2)
        if response == correctResp
            correct = 1;
            feedbackColor = [0 1 0]*white;
            quadrants = [1 2 3 4];
        else
            correct = 0;
            feedbackColor = [1 0 0]*white;
            quadrants = [3 4];
        end
        
        % update accuracy based on target
        % for display
        blockAcc(respInterval) = blockAcc(respInterval)+correct;
        
        % for titration
        % titration flag
        if cuedInterval==3 && p.isT1mid
            neutrAcc(respInterval) = neutrAcc(respInterval)+correct;
            neutrnumMat(respInterval) = neutrnumMat(respInterval) + 1;
        end
        
        Screen('DrawDots', window, plHoldPos, plHoldSize, black, center, 2); % placeholders
        putFeedback(window, [cx, cy], quadrants, feedbackColor, p.feedbackSize, 2)
        timeFeedback = Screen('Flip', window);
        
        if p.eyeTracking
            Eyelink('Message', 'EVENT_FEEDBACK');
        end
        
        if cuedInterval==3
            cueValidity = 0;
        else
            cueValidity = 1;
        end
        
        
        % save to the txt file
        fprintf(fp, ['%.f\t %.f\t %.2f\t %.2f\t',...
            '%.0f\t %.0f\t %.0f\t',...
            '%.2f\t %.2f\t %.0f\t %.0f\t %.0f\t %.0f\t',...
            '%.7f\t %.0f\t %.0f\t %.0f\t %.5f\t %.5f\t %.5f \n'],...
            iTrial, bb, p.stairst1.xCurrent, p.stairst2.xCurrent,...
            respInterval, cueValidity, cuedInterval,...
            angles(1), angles(2), rots(1), rots(2), ph(1), ph(2),...
            rt, responseKey, response, correct, currITI, readyT1SOA, readyT2SOA);
        
        % Store timing
        timing.readyT1SOA(iTrialAll)    = readyT1SOA;
        timing.readyT2SOA(iTrialAll)    = readyT2SOA;
        timing.timeFix(iTrialAll)       = timeFix;
        timing.timeCue(iTrialAll)       = timeStart;
        timing.timeIm1(iTrialAll)       = timeIm1;
        timing.timeBlank1(iTrialAll)    = timeBlank1;
        timing.timeIm2(iTrialAll)       = timeIm2;
        timing.timeBlank2(iTrialAll)    = timeBlank2;
        timing.timeRespCue(iTrialAll)   = timeRespCue;
        timing.timeGoCue(iTrialAll)     = timeGoCue;
        timing.timeFeedback(iTrialAll)  = timeFeedback;
        timing.timeITI(iTrialAll)       = timeITI;
        timing.timePrecue(iTrialAll)    = timePrecue;
        
        save(['/users/purplab/Desktop/Aysun/endoTAU/data/',matFileName]) % saves the workspace on each trial
        
        %%%%%%%%%%%%%%%%%
        stairValues(iTrialAll,:) = [p.stairst1.xCurrent, p.stairst2.xCurrent];
        
        % Adjust staircase level
        if p.staircase && p.isT1mid
            if (respInterval == 1) 
                p.stairst1 = usePalamedesStaircase(p.stairst1, correct); %AF
            elseif respInterval == 2 %AF
                p.stairst2 = usePalamedesStaircase(p.stairst2, correct); %AF
            end
        end
%         fprintf('%f %f \n',p.stairst1.xCurrent, p.stairst2.xCurrent)
        
        % Update the overall trial number
        iTrialAll = iTrialAll+1;
        Screen('Close', tex1)
        Screen('Close', tex2)
    end % end of the block
    
    %%%%%%%%% Before moving onto next block: %%%%%%%%%
    blockAccPerc = sum(blockAcc,1)./size(trials,1)*2*100; % T1 & T2 accuracy percentage in this block
    neutrAccPerc = neutrAcc./neutrnumMat
    neutrAccPerc_titrate = neutrAccPerc_titrate + neutrAccPerc;
    
    overallAcc = overallAcc + neutrAccPerc;
    
    accMessage = sprintf('Accuracy: %.2f%%', mean(blockAccPerc));
    blockMessage = sprintf('%s You''ve completed %d of %d blocks.', highpraise, bb, p.nBlocks);
    
    if bb==p.nBlocks
        keyMessage = '';
    else
        keyMessage = 'Please take a break! Press any key to go on when you are ready.';
    end
    
    % Save workspace if practice block (since we might quit after just
    % one block)
    if p.isPractice
        save(sprintf('/users/purplab/Desktop/Aysun/endoTAU/data/WORKSPACE_%.f', subjectID))
        breakMessage = sprintf('%s\n\n%s', accMessage, keyMessage);
        DrawFormattedText(window, breakMessage, 'center', 'center', [1 1 1]*white);
        Screen('Flip', window);
    else
        breakMessage = sprintf('%s\n%s\n\n%s', blockMessage, accMessage, keyMessage);
        DrawFormattedText(window, breakMessage, 'center', 'center', [1 1 1]*white);
        Screen('Flip', window);
    end
    
    if (sesT == 1) && (bb==5) % experimental session
        neutrAccPerc_titrate = neutrAccPerc_titrate/5
        
        % copied from Antonio Fernandez - temporal fields experiment
        labelStairs = {'stairst1', 'stairst2'};
        fprintf('before = %.2f %.2f', p.stairst1.xCurrent, p.stairst2.xCurrent)
        for targIndex = 1:2
            adjustFlag = 1; % AD - to overcome double adjustments
            
            if adjustFlag == 1
                if p.(labelStairs{targIndex}).xCurrent <= 1.5
                    if neutrAccPerc_titrate(targIndex) > .90
                        p.(labelStairs{targIndex}).xCurrent = p.(labelStairs{targIndex}).xCurrent-0.5;
                        adjustFlag = 0;
                    elseif neutrAccPerc_titrate(targIndex) > .82
                        p.(labelStairs{targIndex}).xCurrent = p.(labelStairs{targIndex}).xCurrent -0.25;
                        adjustFlag = 0;
                    elseif neutrAccPerc_titrate(targIndex) < .63
                        % p.(labelStairs{targIndex}).xCurrent = p.(labelStairs{targIndex}).xCurrent +0.2;
                        p.(labelStairs{targIndex}).xCurrent = p.(labelStairs{targIndex}).xCurrent +0.25; % aysun
                        adjustFlag = 0;
                    else
                        p.(labelStairs{targIndex}).xCurrent = p.(labelStairs{targIndex}).xCurrent;
                        adjustFlag = 0;
                    end
                end
            end
            if adjustFlag == 1
                if p.(labelStairs{targIndex}).xCurrent > 1.5 &&  p.(labelStairs{targIndex}).xCurrent <= 3
                    if neutrAccPerc_titrate(targIndex) > .90
                        p.(labelStairs{targIndex}).xCurrent =  p.(labelStairs{targIndex}).xCurrent -1;
                        adjustFlag = 0;
                    elseif neutrAccPerc_titrate(targIndex) > .82
                        p.(labelStairs{targIndex}).xCurrent =  p.(labelStairs{targIndex}).xCurrent -.5;
                        adjustFlag = 0;
                    elseif neutrAccPerc_titrate(targIndex) < .63
                        p.(labelStairs{targIndex}).xCurrent =  p.(labelStairs{targIndex}).xCurrent +.25;
                        adjustFlag = 0;
                    else
                        p.(labelStairs{targIndex}).xCurrent =  p.(labelStairs{targIndex}).xCurrent;
                        adjustFlag = 0;
                    end
                end
            end
            if adjustFlag == 1
                if  p.(labelStairs{targIndex}).xCurrent > 3
                    if neutrAccPerc_titrate(targIndex) > .90
                        p.(labelStairs{targIndex}).xCurrent =  p.(labelStairs{targIndex}).xCurrent -1.25;
                    elseif neutrAccPerc_titrate(targIndex) > .82
                        p.(labelStairs{targIndex}).xCurrent =  p.(labelStairs{targIndex}).xCurrent -0.75;
                    elseif neutrAccPerc_titrate(targIndex) < .60
                        p.(labelStairs{targIndex}).xCurrent =  p.(labelStairs{targIndex}).xCurrent + .5;
                    elseif neutrAccPerc_titrate(targIndex) < .65  && neutrAccPerc_titrate(targIndex) > .60
                        p.(labelStairs{targIndex}).xCurrent =  p.(labelStairs{targIndex}).xCurrent +0.25;
                    else
                        p.(labelStairs{targIndex}).xCurrent  = p.(labelStairs{targIndex}).xCurrent;
                    end
                end
                adjustFlag = 0;
            end
            
            % AD - to prevent titration being off the normal range
            if p.(labelStairs{targIndex}).xCurrent > 5
                p.(labelStairs{targIndex}).xCurrent = 4;
            end
            
            % AD - to prevent titration being off the normal range
            if p.(labelStairs{targIndex}).xCurrent < 0.25
                p.(labelStairs{targIndex}).xCurrent = 0.75;
            end
        end
        
        if p.stairst1.xCurrent == p.stairst2.xCurrent
            p.stairst1.xCurrent = p.stairst1.xCurrent - 0.15;
        end
        fprintf('after = %.2f %.2f', p.stairst1.xCurrent, p.stairst2.xCurrent)
        
        neutrAccPerc_titrate = [0 0];
    end
    WaitSecs(1);
    if bb < p.nBlocks
       
        % Collect response
        responseKey = [];
        while isempty(responseKey) % record wrong key as missed trial
            [secs, keyCode] = KbWait(devNum);
            rt = secs - timeRespCue;
            
            if numel(find(keyCode))>1 % if more than one key was pressed simultaneously
                responseKey = [];
            else
                responseKey = find(p.keyCodes(1:4)==find(keyCode));
                if responseKey == 3 % if p pressed
                    ShowCursor
                    save(['/users/purplab/Desktop/Aysun/endoTAU/data/',matFileName]) % saves the workspace on each trial
                    % Save eye data and shut down the eye tracker
                    if p.eyeTracking
                        rd_eyeLink('eyestop', window, {eyeFile, eyeDataDir});
                        
                        % rename eye file
                        eyeFileFull = sprintf('%s/%.fs%.f_%.f_v3.edf', eyeDataDir, sesT, subjectID,sessionNum);
                        copyfile(sprintf('%s/%s.edf', eyeDataDir, eyeFile), eyeFileFull)
                    end
                    
                    % Clean up
                    PsychPortAudio('Stop', pahandle);
                    PsychPortAudio('Close', pahandle);
                    Screen('CloseAll')
                    ShowCursor
                    break
                end
            end
        end
        
        
    end
end
timing.endTime = GetSecs;

DrawFormattedText(window, 'All done! Thanks for your effort!', 'center', 'center', [1 1 1]*white);
Screen('Flip', window);
WaitSecs(2);
fprintf('\noverall neutral accuracy is:\n')
overallAcc = overallAcc./p.nBlocks

%% Save eye data and shut down the eye tracker
if p.eyeTracking
    rd_eyeLink('eyestop', window, {eyeFile, eyeDataDir});
    
    % rename eye file
    eyeFileFull = sprintf('%s/%.fs%.f_%.f_v3.edf', eyeDataDir, sesT, subjectID,sessionNum);
    copyfile(sprintf('%s/%s.edf', eyeDataDir, eyeFile), eyeFileFull)
end

%% Clean up
PsychPortAudio('Stop', pahandle);
PsychPortAudio('Close', pahandle);
Screen('CloseAll')
ShowCursor

%% Plot for thresholding

figure
subplot(2,1,1)
plot(1:numel(stairValues)/2, stairValues(:,1), 'm*')
set(gca,'YTick', [0.25:0.25:3, 3.5:0.5:7])
grid on
legend('T1')
title('Stair Values - check if they are asymptoting')

subplot(2,1,2)
plot(1:numel(stairValues)/2, stairValues(:,2), 'm*')
set(gca,'YTick', [0.25:0.25:3, 3.5:0.5:7])
grid on
legend('T2')

xlabel('time')
ylabel('stairs (angle)')
figName = sprintf('tilts_s%.0f_%.f_ses%.0f.jpg', subjectID, sesT,sessionNum);
saveas(gcf,figName);

%% Store expt info
expt.subjectID = subjectID;
expt.p = p;
expt.timing = timing;
expt.trialOrders = trialOrder;
expt.trials_headers = trials_headers;
expt.trials = trials;

if p.staircase
    expt.staircase.stairValues = stairValues;
    expt.staircase.threshold = [p.stairst1.xCurrent, p.stairst2.xCurrent];
end

if p.eyeTracking
    expt.eye.fixCue = fixCue;
    expt.eye.fixT1 = fixT1;
    expt.eye.fixT2 = fixT2;
    expt.eye.fixRC = fixRC;
end

%% Save data
fclose(fp);
% save(matFileName, 'expt')
save(['/users/purplab/Desktop/Aysun/endoTAU/data/',matFileName]) % saves the workspace on each trial
        
%% Plot timings
figure
subplot(1,2,1), plot(sort(expt.timing.readyT1SOA)), title('T1 onset')
subplot(1,2,2), plot(sort(expt.timing.readyT2SOA)), title('T2 onset')

timinzplotname = sprintf('targOnsets_s%.f_ses%.f.png', subjectID,sessionNum);
saveas(gcf,timinzplotname)

%% Output staircase
if p.staircase
    fprintf('Stair level is: %f %f\n ( %s)', expt.staircase.threshold, accMessage);
end

end