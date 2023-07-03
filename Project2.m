%Part 1: Matched filters and correlators in noise free environment

%Initialization of variables
A = 1;        %trasmitted voltage level
bit_time = 1; %pulse width to transmit one bit (Ts = 1 second)
sampling_time = 0.200; %taking 5 samples for each bit, T= Ts/5 = 0.2 sec
num_samples_per_bit = floor (bit_time / sampling_time); %number of samples for 1 bit = 5
num_bits = 10; %number of bits in stream
num_samples = num_bits * num_samples_per_bit; %total number of samples in stream

p = [5 4 3 2 1]/sqrt(55); %pulse shaping filter with unity energy

data_bits = randi([0,1], 1, num_bits); %generating a stream of random bits
data_symbols = (2*data_bits - 1) * A; %Transmitting 1 as 1 V and 0 as -1 V
%each bit is represented by 5 samples,1 of them is either 1 or -1 
%and the rest should be 0 so we used upsample function to get the required
%form
impulses = upsample(data_symbols, num_samples_per_bit);
%define a variable (t) represents the x-axis starts from 0 and has a value with
%each sampling time until reaching the number of samples-1 as we start from 0
t = 0 : sampling_time : (num_samples - 1)*sampling_time;
%plotting the stream of symbols of bits generated before 
figure ('Name', 'Impulse Train');
stem (t, impulses);
title ("Impulse Train");
xlabel("Time(Sec)");
ylabel("Volts(V)");

%applying pulse shaping on symbols to transmit the stream by convoloving
%the pulse shape p (triangular shape) with the symbols
y = conv(impulses, p);
%adjusting the size of transmitted signal because the convolution changed
%the size of stream
y = y(1 : num_samples);
%plotting the transmitted stream after pulse shaping
figure ('Name', 'Transmitted Signal');
plot (t, y);
title ("Tranmitted Signal");
xlabel ("Time(Sec)");
ylabel ("Volts(V)");

%matched filter is the same as p but flipped
matched_filter = fliplr(p);
%received sequence after convolution with matched filter
matched_Rx = conv (y, matched_filter);

%rect filter after energy normalization
rect_filter = [1 1 1 1 1] / sqrt(5);
%received sequence after convolution with rect filter
rect_Rx = conv (y, rect_filter);
%defining a new variable (t_conv) to represent the x-axis starting from 
%0.2 sec with 0.2 sec step until it reaches the size of rect_Rx (the
%convolution result) 
t_conv = 0.2 : sampling_time : ((length(y)) + length(p) - 1)* sampling_time;
%plotting the recieved signal after passing by matched filter
figure ('Name', 'Received Signal');
subplot(2,1,1);
plot (t_conv, matched_Rx,'r');
title ("Matched Filter Output");
xlabel ("Time(Sec)");
ylabel ("Volts(V)");

%plotting the recieved signal after passing by rect filter
subplot(2,1,2);
plot (t_conv, rect_Rx,'b');
title ("Rectangular Filter Output");
xlabel ("Time(Sec)");
ylabel ("Volts(V)");

%Passing the receieved signal by a correlator instead of matched filter
%initialization of received signal after passing by correlator
corr_Rx = zeros(size(matched_Rx)); 
%repeat the pulse shape function by number of bits
p_correlator = repmat(p, 1, num_bits); 
accumlator = 0; %initialization of integrator result

%for loop to multiply the received sequence by the repeated pulse shape
%function and integrating the multiplication for each 5 samples
for i = 1 : num_samples
    if mod((i-1), 5) == 0
        accumlator = 0;
    end
    corr_Rx(i) = accumlator + (p_correlator(i) * y(i));
    accumlator = corr_Rx(i);
end

%plotting both results of matched filter and correlator in the same graph
%with different colours 
figure ('Name', 'Matched Filter Output & Correlator Output');
plot(t_conv, matched_Rx, 'r', t_conv, corr_Rx, 'b');
title ("Matched Filter Output & Correlator Output");
xlabel ("Time(Sec)");
ylabel("Volts(V)");
legend('Matched Filter Output', 'Correlator Output');

%Part 2: Noise analysis

num_bits_noise = 10000; %number of bits in stream
num_samples_noise = num_bits_noise * num_samples_per_bit; %total number of samples in stream
%generating a stream of random bits
data_bits_noise = randi([0,1], 1, num_bits_noise);
%Transmitting 1 as 1 V and 0 as -1 V
data_symbols_noise = (2*data_bits_noise - 1) * A;
%each bit is represented by 5 samples,1 of them is either 1 or -1 
%and the rest should be 0 so we used upsample function to get the required
%form
impulses_noise = upsample(data_symbols_noise, num_samples_per_bit);

%applying pulse shaping on symbols to transmit the stream by convoloving
%the pulse shape p (triangular shape) with the symbols
s = conv(impulses_noise, p); 
%adjusting the size of transmitted signal because the convolution changed
%the size of stream
s = s(1 : num_samples_noise);
%generating normally distributed noise with zero mean and variance = 1
noise = randn(1, num_samples_noise);
%scaled noise initialization
scaled_noise = zeros(8,num_samples_noise);
%Initialization of variables,they are all array of 8 vectors because we
%will try 8 values of SNR [-2,5]
%the transmitted signal with noise added to it after passing by matched filter
matched_Rx_noise = zeros(8,length(s)+ length(matched_filter)-1);
%samples of the matched filter's output at sampling instants
sampled_Rx_matched = zeros(8, num_bits_noise);
%samples of the matched filter's output decoded into bits (0's and 1's)
sampled_Rx_bits_matched = zeros(8, num_bits_noise);
%an array of vectors that identifies the error occured in which bits for
%each stream that passed by the matched filter
error_matched = zeros(8, num_bits_noise);
%bit error rate for each stream that passed by the matched filter
BER_matched = zeros(1,8);

%rect filter after energy normalization
rect_filter_noise = [5 5 5 5 5]/sqrt(125);
%the transmitted signal with noise added to it after passing by rect filter
rect_Rx_noise = zeros(8,length(s)+ length(matched_filter)-1);
%samples of the rect filter's output at sampling instants
sampled_Rx_rect = zeros(8, num_bits_noise);
%samples of the rect filter's output decoded into bits (0's and 1's)
sampled_Rx_bits_rect = zeros(8, num_bits_noise);
%an array of vectors that identifies the error occured in which bits for
%each stream that passed by the rect filter
error_rect = zeros(8, num_bits_noise);
%bit error rate for each stream that passed by the rect filter
BER_rect = zeros(1,8);

Eb = 1; % bit energy which is a unity
No = zeros (1,8); %initialization of No used in the variance of noise (No/2)
%for loop that gets the required BER from the range of SNR [-2,5] db
for SNR = -2 : 5
    No(SNR+3) = Eb/(10^(SNR/10));
    %making the noise with variance No/2 by multipying it by sqrt(No/2)
    scaled_noise(SNR+3,:) = noise*sqrt(No(SNR+3)/2);
    %transmitted signal with noise added to it
    v = s + scaled_noise(SNR+3,:);
    %passing the received signal by both filters (matched and rect)
    matched_Rx_noise(SNR+3,:) = conv (v, matched_filter);
    rect_Rx_noise (SNR+3, :) = conv (v, rect_filter_noise);
    %getting the symbols from both filters (matched and rect),
    %The sampling instants occurs each five symbols 
    for i = 5 : 5 : num_samples_noise
        sampled_Rx_matched(SNR+3,i/5) = matched_Rx_noise(SNR+3,i);
        sampled_Rx_rect(SNR+3,i/5) = rect_Rx_noise (SNR+3, i);
    end
    %Decoding sampled symbols of both filters into bits by comparing with zero to compare it
    %with the original stream of bits that has been transmitted
    for i = 1: num_bits_noise
        if sampled_Rx_matched(SNR+3,i) > 0
            sampled_Rx_bits_matched(SNR+3,i) = 1 ;
        else
            sampled_Rx_bits_matched(SNR+3,i) = 0 ;
        end
        
        if sampled_Rx_rect(SNR+3,i) > 0
            sampled_Rx_bits_rect(SNR+3,i) = 1 ;
        else
            sampled_Rx_bits_rect(SNR+3,i) = 0 ;
        end
    end
    %recognize the bits that differed between transmitted and received streams 
    %for both filters 
    error_matched(SNR+3,:) = sampled_Rx_bits_matched(SNR+3,:) == data_bits_noise;
    %counting the number of flipped bits and get the ratio
    %between num of bits that was flipped and the total number of bits in
    %each stream for each SNR to get bit error rate for each SNR value
    BER_matched(SNR+3) = num_bits_noise - nnz(error_matched(SNR+3,:));
    BER_matched(SNR+3) = BER_matched(SNR+3)/(num_bits_noise) ;
    
    error_rect(SNR+3,:) = sampled_Rx_bits_rect(SNR+3,:) == data_bits_noise;
    BER_rect(SNR+3) = num_bits_noise - nnz(error_rect(SNR+3,:));
    BER_rect(SNR+3) = BER_rect(SNR+3)/(num_bits_noise) ;
    
end
%theoretical value of bit error rate for a matched filter
BER_Theoritical = 0.5 * erfc(sqrt(Eb./No));

%plotting bit error rate with log scale versus SNR for both filters and
%drawing them on the same graph along with the theoretical BER
figure('Name', 'BER Vs SNR');
semilogy(10*log10(Eb./No), BER_Theoritical, 'r', 10*log10(Eb./No), BER_matched, 'b', 10*log10(Eb./No), BER_rect, 'k');
title ("BER Vs SNR");
xlabel("SNR Eb/N0 (db)");
ylabel("BER");
xlim([-2,8]);
legend ('BER Theoretical', 'BER Matched Filter', 'BER Rect Filter')
grid on;

%Part 3: ISI in noise free channel

num_bits_ISI = 100; %number of bits in stream
input_data = randi([0 1], 1, num_bits_ISI); %generating a stream of random bits
Tx = ((2*input_data)-1) * A; %Transmitting 1 as 1 V and 0 as -1 V

%each bit is represented by 5 samples,1 of them is either 1 or -1 
%and the rest should be 0 so we used upsample function to get the required
%form
Tx_Sampeled = upsample(Tx, num_samples_per_bit);

R = [0 0 1 1]; %Different rolloff factors
delay = [2 8 2 8]; %Different delays for ssrc filters

for i = 1:4
    srrc = rcosine(1, num_samples_per_bit, 'sqrt', R(i), delay(i)); %pulse shaping filter
    srrc2 = conv(srrc, srrc);
    %Plot the different 4 Filters
    t_isi = 0 : length(srrc)-1;
    figure ('Name', 'SRRC Filters');
    plot(t_isi, srrc);
    title(strcat("SRRC of R = ", {' '}, string(R(i)), ' delay = ', {' '}, string(delay(i))));
    
    figure ('Name', 'SRRC Convolution');
    plot(srrc2);
    title(strcat("SRRC Conv of R = ", {' '}, string(R(i)), ' delay = ', {' '}, string(delay(i))));
    
    transmitted = conv(Tx_Sampeled, srrc); %Generate the transmitted Signal
    received = conv(transmitted, srrc); %Generarate the received signal
    
    eyediagram(transmitted, 2*num_samples_per_bit); %Eye diagram at point A (Transmitter) 
    title(strcat('transmited siganl when R = ',{' '},string(R(i)),{newline},' delay = ',{' '},string(delay(i))));
    
    eyediagram(received, 2*num_samples_per_bit); %Eye diagram at point B (Receiver)
    title(strcat('received siganl when R = ',{' '},string(R(i)),{newline},' delay = ',{' '},string(delay(i))));
end