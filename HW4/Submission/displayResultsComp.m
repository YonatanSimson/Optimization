function displayResultsComp(h1, h2, fsamp)
    figure;
    subplot(2,1,1)
    [H1, f] = freqz(h1,1,1024,fsamp);
    [H2, f] = freqz(h2,1,1024,fsamp);
    
    semilogy(f,abs(H1), 'b'), grid on
    hold on;
    semilogy(f,abs(H2), 'g'), grid on
    hold off;
    title('Frequency responce F\{h_d[n]\}')
    ylabel('|H_d(\omega)|')
    xlabel('Frequency[Hz]')
    legend('H1', 'H2');
    
    subplot(2,1,2)
    stem(0:length(h1)-1, h1, 'b');
    hold on;
    stem(0:length(h2)-1, h2, 'g');
    hold off;
    title('Impulse responce')
    ylabel('h_d[n]');
    xlabel('n');
    legend('h1', 'h2');

end