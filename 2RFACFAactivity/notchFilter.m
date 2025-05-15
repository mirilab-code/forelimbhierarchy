function signal_filt = notchFilter(signal, Fs)
    % Fs is the sampling rate of the signal
    Fn = Fs/2;                              
    Wp = [59 61]/Fn;                          % Stopband Frequency (Normalised)
    Ws = [58 62]/Fn;                          % Passband Frequency (Normalised)
    Rp =  1;                                  % Passband Ripple
    Rs = 60;                                  % Passband Ripple (Attenuation)
    [n,Wp] = ellipord(Wp,Ws,Rp,Rs);           % Elliptic Order Calculation
    [z,p,k] = ellip(n,Rp,Rs,Wp,'stop');       % Elliptic Filter Design: Zero-Pole-Gain
    [sos,g] = zp2sos(z,p,k);                  % Second-Order Section For Stability
    % figure
    % freqz(sos, 2^18, Fs)                    % Filter Bode Plot
    % set(subplot(2,1,1), 'XLim',[0 100])
    % set(subplot(2,1,2), 'XLim',[0 100])
    signal_filt = filtfilt(sos, g, signal); 
end