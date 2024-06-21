function [tp] = pickps(seis, fs,sta_win,lta_win,th)
    
% Input:    
    %   - samplingRate: Sampling rate of the seismogram in Hz
    
    % Define parameters
    %sta_win = 0.5;  % Short-Term Average window size in seconds
    %lta_win = 10; % Long-Term Average window size in seconds
    %th = 2.5;     % Threshold for onset detection
    
    % Calculate the number of samples for each window
    sta_n = round(sta_win * fs);
    lta_n = round(lta_win * fs);
    
    % Compute STA/LTA ratio
    sta = movstd(seis, [sta_n 0]);
    lta = movstd(seis, [lta_n 0]);
    ratio = sta ./ lta;
    
    % Thresholding to detect onset
    onsetIndices = find(ratio > th);
    
    % Pick the first onset index as P-wave arrival time
    if ~isempty(onsetIndices)
        tp = onsetIndices(1) / fs;
        %disp(['arrival time: ', num2str(tp), ' seconds : ths ', num2str(th)]);
    else
        %disp('arrival not detected.');
        tp = NaN;
    end
end
