% main_progetto1.m - Progetto 1 (MATLAB-only, NO toolbox required)
clear; clc; close all;

% ===== 1) Testo da trasmettere =====
tx_text = "Ciao! Questo e' un test per Progetto 1 - Sistemi di Comunicazione.";

% ===== 2) Testo -> bit =====
bits  = text_to_bits(tx_text);
Nbits = numel(bits);

% ===== 3) Parametri =====
SNRdB = 0:2:14;           % valori SNR da simulare
repN  = 3;                % codifica semplice: ripetizione 3x

Nframes = 50;             % <-- NEW: media su piu' frame (grafici piu' stabili)

modes     = ["AWGN","RAYLEIGH","OSTACOLO"];   % <-- NEW: 3 ambienti
useCoding = [false true];                     % senza / con codifica

BER = zeros(numel(useCoding), numel(modes), numel(SNRdB));

% ===== 4) Simulazione BER (con averaging) =====
for c = 1:numel(useCoding)
    for m = 1:numel(modes)
        for k = 1:numel(SNRdB)

            errSum = 0;
            bitSum = 0;

            for f = 1:Nframes
                b = bits;

                % --- codifica (opzionale) ---
                if useCoding(c)
                    b = rep_encode(b, repN);
                end

                % --- modulazione BPSK ---
                s = bpsk_mod(b);  % +/-1

                % --- canale (ambiente) ---
                if modes(m) == "AWGN"
                    r = channel_awgn(s, SNRdB(k));
                elseif modes(m) == "RAYLEIGH"
                    r = channel_rayleigh(s, SNRdB(k));
                else
                    % OSTACOLO: extra attenuazione (es. muro) -> SNR effettivo piu' basso
                    r = channel_ostacolo(s, SNRdB(k), 6);  % 6 dB loss (puoi cambiare 6)
                end

                % --- demodulazione ---
                bhat = bpsk_demod(r);

                % --- decodifica (opzionale) ---
                if useCoding(c)
                    bhat = rep_decode(bhat, repN);
                end

                bhat = bhat(1:Nbits); % stessa lunghezza del testo originale

                errSum = errSum + sum(bits ~= bhat);
                bitSum = bitSum + Nbits;
            end

            BER(c,m,k) = errSum / bitSum;
        end
    end
end

% ===== 5) Esempio ricostruzione testo (SNR=12 dB per evitare errori) =====
testSNR = 12;
s = bpsk_mod(bits);
r = channel_awgn(s, testSNR);
rx_bits = bpsk_demod(r);
rx_text = bits_to_text(rx_bits(1:Nbits));

fprintf("TX: %s\n", tx_text);
fprintf("RX: %s\n", rx_text);

% ===== 6) Grafico BER =====
figure; hold on; grid on;
for c = 1:numel(useCoding)
    for m = 1:numel(modes)
        semilogy(SNRdB, squeeze(BER(c,m,:)), '-o', 'DisplayName', ...
            sprintf('%s - %s', modes(m), tern(useCoding(c),'con codifica','senza codifica')));
    end
end
xlabel('SNR (dB)'); ylabel('BER');
title('Progetto 1 - BER vs SNR (MATLAB-only)');
legend('Location','southwest');

% ===== 7) Salvataggio risultati =====
saveas(gcf, "BER_vs_SNR.png");
saveas(gcf, "BER_vs_SNR.fig");
save("results.mat","BER","SNRdB","modes","useCoding","repN","Nframes");

%% ===== Funzioni (NO TOOLBOX) =====

function out = tern(cond, a, b), if cond, out=a; else, out=b; end, end

% --- Text <-> Bits (NO de2bi/bi2de) ---
function bits = text_to_bits(txt)
    bytes = uint8(char(txt));          % string -> uint8 bytes
    bits8 = zeros(numel(bytes), 8);    % each byte -> 8 bits
    for i = 1:numel(bytes)
        b = bytes(i);
        for k = 1:8
            bits8(i,k) = bitget(b, 9-k);   % MSB -> LSB
        end
    end
    bits = reshape(bits8.', 1, []);
    bits = logical(bits);
end

function txt = bits_to_text(bits)
    bits = bits(:).';
    n = floor(numel(bits)/8)*8;
    bits = bits(1:n);
    bits8 = reshape(bits, 8, []).';     % each row is 1 byte

    bytes = zeros(size(bits8,1),1,'uint8');
    for i = 1:size(bits8,1)
        b = uint8(0);
        for k = 1:8
            if bits8(i,k)
                b = bitset(b, 9-k, 1);  % MSB -> LSB
            end
        end
        bytes(i) = b;
    end
    txt = char(bytes).';   % row char vector
end

% --- BPSK ---
function s = bpsk_mod(bits)
    s = 2*double(bits) - 1;  % 0->-1, 1->+1
end

function bits = bpsk_demod(r)
    bits = r >= 0;
end

% --- Channels ---
function y = channel_awgn(x, SNRdB)
    Ps = mean(abs(x).^2);
    SNR = 10^(SNRdB/10);
    N0 = Ps/SNR;
    noise = sqrt(N0/2) * randn(size(x));
    y = x + noise;
end

function y = channel_rayleigh(x, SNRdB)
    % Flat Rayleigh fading: r = h*x + n, assume perfect equalization
    h = (randn(size(x)) + 1j*randn(size(x))) / sqrt(2);
    faded = h .* x;

    Ps = mean(abs(faded).^2);
    SNR = 10^(SNRdB/10);
    N0 = Ps/SNR;
    n = sqrt(N0/2) * (randn(size(x)) + 1j*randn(size(x)));
    r = faded + n;

    y = real(r ./ h); % perfect equalization (assumed)
end

function y = channel_ostacolo(x, SNRdB, loss_dB)
    % Ostacolo = extra attenuation -> effective SNR decreases
    y = channel_awgn(x, SNRdB - loss_dB);
end

% --- Simple coding: repetition ---
function enc = rep_encode(bits, N)
    enc = repelem(bits, N);
end

function dec = rep_decode(bits, N)
    bits = bits(:).';
    L = floor(numel(bits)/N)*N;
    bits = bits(1:L);
    M = reshape(bits, N, []);
    dec = sum(M,1) >= (N/2);
    dec = logical(dec);
end