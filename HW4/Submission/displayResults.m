function displayResults(hh, fsamp)
    figure(1);
    subplot(2,1,1)
    [H,f] = freqz(hh,1,1024,fsamp);
    semilogy(f,abs(H)), grid on
    title('Frequency responce F\{h_d[n]\}')
    ylabel('|H_d(\omega)|')
    xlabel('Frequency[Hz]')

    subplot(2,1,2)
    stem(0:length(hh)-1, hh);
    title('Impulse responce')
    ylabel('h_d[n]');
    xlabel('n');

end