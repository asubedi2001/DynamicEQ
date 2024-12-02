clc; clear; close all

% MATLAB toolkits required:
% - Communications Toolbox
% - Audio Toolbox

% human hearing limits: 20 hz to 20khz
minFrequency = 0;
maxFrequency = 22050;

% how many samples to process at a time
frameLength = 1024;

% audio reader:
songReader = dsp.AudioFileReader('../audio/porter.wav', 'SamplesPerFrame', frameLength);
sampleRate = songReader.SampleRate;
disp(sampleRate)
% writer -- outputs to speakers:
deviceWriter = audioDeviceWriter('SampleRate', songReader.SampleRate);

% create frame buffer
% 43.07 frames a second * 2 -> buffer is around ~1.857s long (if bufferSize = 80 and doesnt change)
bufferSize = 80;
inputFrameBuffer = zeros(frameLength, bufferSize);
outputFrameBuffer = zeros(frameLength, bufferSize);
bufferIdx = 1;

% used to display sound waves while playing:
scope = timescope( ...
    'SampleRate', songReader.SampleRate, ...
    'TimeSpan', 2, ...
    'BufferLength', songReader.SampleRate*2*2, ...
    'YLimits', [-1,1], ...
    'TimeSpanOverrunAction', "Scroll");

% plot
freqFigure = figure();
outputFreqX = zeros(1);
outputFreqY = zeros(1);
outputPlot = plot(outputFreqX, outputFreqY, 'DisplayName', 'Output');
hold on
inputFreqX = zeros(1);
inputFreqY = zeros(1);
inputPlot = plot(inputFreqX, inputFreqY, 'DisplayName', 'Input');
xscale log
title("Frequency Domain Graph")
xlabel("Frequency (Hz)")
ylabel("Amplitude")
outputPlot.XDataSource = 'outputFreqX';
outputPlot.YDataSource = 'outputFreqY';
inputPlot.XDataSource = 'inputFreqX';
inputPlot.YDataSource = 'inputFreqY';
legend
hold off

% code from 
% https://stackoverflow.com/questions/61369933/how-to-build-a-bandpass-filter-in-matlab-with-the-butter-function
% https://www.mathworks.com/matlabcentral/answers/1890802-quickest-way-to-do-a-bandpass-filter
f_low = 1; % Hz
f_high = 800; % Hz
f_sampling = songReader.SampleRate;
f_nrm_low   = f_low /(f_sampling/2);
f_nrm_high  = f_high /(f_sampling/2);
% determine filter coefficients:
[z,p,k] = butter(4,[f_nrm_low f_nrm_high],'bandpass');
% convert to zero-pole-gain filter parameter (recommended)
sos = zp2sos(z,p,k);

inputEq = graphicEQ('Structure','Cascade', 'SampleRate', songReader.SampleRate);
inputEq.Gains = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

gainsFigure = figure();
gains = zeros(10, 1);
gainsPlot = plot(gains, 'DisplayName', 'Gains');
hold on
title("Equalizer Gains")
gainsPlot.YDataSource = 'inputEq.Gains';
ylim([-1 inf])
hold off

bassBand = struct( ...
    'octaves', logical([1, 1, 1, 1, 0, 0, 0, 0, 0, 0]), ...
    'freqRange', [22, 355], ...
    'inOutFreqDiffs', [], ...
    'outFreqSums', [], ...
    'weightedAvg', [], ...
    'numFrames', 2, ...
    'decimation', 4);

midBand = struct( ...
    'octaves', logical([0, 0, 0, 0, 1, 1, 1, 1, 0, 0]), ...
    'freqRange', [355, 5623], ...
    'inOutFreqDiffs', [], ...
    'outFreqSums', [], ...
    'weightedAvg', [], ...
    'numFrames', 1, ...
    'decimation', 2);

trebBand = struct( ...
    'octaves', logical([0, 0, 0, 0, 0, 0, 0, 0, 1, 1]), ...
    'freqRange', [5623, 22387], ...
    'inOutFreqDiffs', [], ...
    'outFreqSums', [], ...
    'weightedAvg', [], ...
    'numFrames', 1, ...
    'decimation', 1);

bands = [bassBand midBand trebBand];

frameCount = 0;
framePeriod = 20;
redrawPeriod = 20;
while ~isDone(songReader) && frameCount < 2000 % use frame count for early stopping TODO: remove
    audioFrame = songReader();
    outputFrame = inputEq(audioFrame);
    noiseFrame = sosfilt(sos, wgn(frameLength,1, 2)); % limit bandwidth of white noise
    %noiseFrame = wgn(frameLength,1,-25);
    inputFrame = outputFrame + noiseFrame;

    % update buffers (and index) according to current frame
    inputFrameBuffer(:, bufferIdx) = inputFrame;
    outputFrameBuffer(:, bufferIdx) = outputFrame;
    bufferIdx = mod(bufferIdx, bufferSize) + 1;

    deviceWriter(inputFrame);
    scope(inputFrame)

    % these FFTs are not necessary; only used for plotting:
    [outputFreqX, outputFreqY] = perform_and_plot_fft(outputFrame, sampleRate, maxFrequency);
    [inputFreqX, inputFreqY] = perform_and_plot_fft(inputFrame, sampleRate, maxFrequency);

    for i = 1:numel(bands)
        % some of the combined frames may be 0 to begin with:
        combinedInputFrame = [];
        for j = 1:bands(i).numFrames
            combinedInputFrame = cat(1, inputFrameBuffer(:, mod(bufferIdx - j, bufferSize) + 1), combinedInputFrame);
        end

        combinedOutputFrame = [];
        for j = 1:bands(i).numFrames
            combinedOutputFrame = cat(1, outputFrameBuffer(:, mod(bufferIdx - j, bufferSize) + 1), combinedOutputFrame);
        end
    
        % decimate performs low-pass filter then downsamples:
        decimatedInputFrame = decimate(combinedInputFrame, bands(i).decimation);
        decimatedOutputFrame = decimate(combinedOutputFrame, bands(i).decimation);
        decimatedSampleRate = sampleRate/bands(i).decimation;
        decimatedMaxFrequency = decimatedSampleRate/2;

        [bandInputFreqX, bandInputFreqY] = perform_and_plot_fft(decimatedInputFrame, decimatedSampleRate, decimatedMaxFrequency);
        [bandOutputFreqX, bandOutputFreqY] = perform_and_plot_fft(decimatedOutputFrame, decimatedSampleRate, decimatedMaxFrequency);

        bandFreqs = bandOutputFreqX(:) >= bands(i).freqRange(1) & bandOutputFreqX(:) < bands(i).freqRange(2); % outputFreqX should be same as inputFreqX
        if isempty(bands(i).inOutFreqDiffs)
            bands(i).inOutFreqDiffs = zeros(sum(bandFreqs), 1);
            bands(i).outFreqSums = zeros(sum(bandFreqs), 1);
            bands(i).weightedAvg = zeros(sum(bandFreqs), 1);
        end

        bands(i).inOutFreqDiffs = bands(i).inOutFreqDiffs + (bandInputFreqY(bandFreqs) - bandOutputFreqY(bandFreqs));
        bands(i).outFreqSums = bands(i).outFreqSums + bandOutputFreqY(bandFreqs);

        if frameCount == 0
            % initialize first frame value
            bands(i).weightedAvg = bands(i).inOutFreqDiffs;
        elseif frameCount <= bufferSize
            % simple average over the frames in unfilled buffer
            bands(i).weightedAvg = (inputFrameBuffer(bufferIdx) - outputFrameBuffer(bufferIdx)) / frameCount;
            frameAvg = zeros(size(bands(i).inOutFreqDiffs));
            for j = 1:frameCount
                frameAvg = frameAvg + (inputFrameBuffer(j) - outputFrameBuffer(j));
            end
            bands(i).weightedAvg = frameAvg / frameCount;
        else
            % have importance of older values decay over time
            decayFactor = 0.6; % decay of old values
            bands(i).weightedAvg = decayFactor * (bands(i).inOutFreqDiffs) + (1 - decayFactor) * bands(i).weightedAvg;
        end
    end

    if mod(frameCount, framePeriod) == 0 && frameCount ~= 0
        disp('Calculating new gains')
        % TODO: Reduce gains if significantly higher than noise, but ensure
        % always >= 0
        % TODO: Set max limits for gains to prevent distortion

        % define minimum and maximum value thresholds to help nromalize
        LOWER = -0.75; 
        UPPER = 0.75; 

        for i = 1:numel(bands)
            % normalize against 'arbitrary' values (fine tuned through testing)
            bandAvgNormalized = (bands(i).weightedAvg - LOWER) / (UPPER - LOWER);
            bandAvgNormalized = mean(bandAvgNormalized);
            bandGain = bandAvgNormalized * 12;

            disp(['Updating gains for band ' num2str(i) ':'])
            bands(i).gain = bandGain;
            disp(bands(i).gain);
            inputEq.Gains(bands(i).octaves) = bands(i).gain;

            bands(i).inOutFreqDiffs(:) = 0;
            bands(i).outFreqSums(:) = 0;
        end
    end

    if mod(frameCount, redrawPeriod) == 0
        refreshdata
        drawnow
    end

    frameCount = frameCount + 1;
end

function [freqX, freqY] = perform_and_plot_fft(signal, sampleRate, maxFrequency)
    numSamples = numel(signal);
    [freqX, freqY, scale] = perform_fft(signal, numSamples, sampleRate);
    cutoff = maxFrequency/scale;
    freqX = freqX(1:floor(cutoff));
    freqY = freqY(1:floor(cutoff));
end

function [freqX, freqY, scale] = perform_fft(songSignal, numSamples, sampleRate)
    % based on https://www.mathworks.com/help/matlabmobile/ug/acquire-and-analyze-audio-data.html
    fftData = fft(songSignal);
    signalLength = numSamples;
    freqY = abs(fftData/signalLength); % abs = sqrt of sum of squares
    freqY = freqY(1:floor(signalLength/2)+1);
    freqY(2:end-1) = 2*freqY(2:end-1);
    freqX = sampleRate*(0:(signalLength/2))/signalLength;
    scale = sampleRate/signalLength;
end

release(songReader)
release(deviceWriter)
release(scope)