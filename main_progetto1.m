% main_progetto1.m - Progetto 1 (MATLAB-only, NO toolbox required)
% Studente: Kerem Tokatli
%
% Obiettivo:
% - Trasmissione/ricezione di un testo: testo->bit->(codifica)->BPSK->canale->demod->(decod)->testo
% - Confronto codifiche: NONE / REP3 / HAMMING(7,4)
% - Ambienti: AWGN / Rayleigh / Ostacolo
% - Effetti "hardware": timing offset + CFO + phase offset
% - Correzione in ricezione (preamble-based): stima timing + CFO e compensazione

clear; clc; close all;

%% ===== 1) Testo da trasmettere =====
tx_text = "Ciao! Questo e' un test per Progetto 1 - Sistemi di Comunicazione.";
payload_bits = text_to_bits(tx_text);

%% ===== 2) Parametri simulazione =====
SNRdB   = 0:2:14;
Nframes = 50;

codingModes = ["NONE","REP3","HAMMING74"];
modes = ["AWGN","RAYLEIGH","OSTACOLO","AWGN_SYNC","AWGN_SYNC_CORR"];

loss_dB = 6; % ostacolo: perdita (dB)

% Impairments (effetti hardware)
timingDelay     = 2;           % ritardo in campioni
freqOffsetNorm  = 0.001;       % CFO normalizzato (cicli/campione)
phaseOffsetRad  = 10*pi/180;   % offset di fase

%% ===== 3) Preamble per sincronizzazione =====
% Due metà uguali per stimare CFO (idea Schmidl-Cox semplificata)
preambleLenHalf = 64;
preambleBitsHalf = randi([0 1], 1, preambleLenHalf) > 0;
preambleBits = [preambleBitsHalf, preambleBitsHalf];
preambleSym  = bpsk_mod(preambleBits); % complesso (in pratica reale)

%% ===== 4) Risultati =====
BER = zeros(numel(codingModes), numel(modes), numel(SNRdB));

%% ===== 5) Simulazione BER =====
for ci = 1:numel(codingModes)
    for mi = 1:numel(modes)
        for si = 1:numel(SNRdB)

            errSum = 0;
            bitSum = 0;

            for f = 1:Nframes

                % --- Codifica payload ---
                [codedPayload, meta] = encode_bits(payload_bits, codingModes(ci));

                % --- Frame: preamble + payload codificato ---
                txBitsFrame = [preambleBits, codedPayload];
                txSym = bpsk_mod(txBitsFrame); % BPSK (complesso)

                % --- Canale base ---
                if modes(mi) == "AWGN"
                    rx = channel_awgn(txSym, SNRdB(si));

                elseif modes(mi) == "RAYLEIGH"
                    rx = channel_rayleigh(txSym, SNRdB(si));

                elseif modes(mi) == "OSTACOLO"
                    rx = channel_awgn(txSym, SNRdB(si) - loss_dB);

                elseif modes(mi) == "AWGN_SYNC"
                    % AWGN + impairments, senza correzione
                    rx = channel_awgn(txSym, SNRdB(si));
                    rx = apply_timing_offset(rx, timingDelay);
                    rx = apply_cfo_phase(rx, freqOffsetNorm, phaseOffsetRad);

                else
                    % AWGN + impairments, con correzione basata su preamble
                    rx = channel_awgn(txSym, SNRdB(si));
                    rx = apply_timing_offset(rx, timingDelay);
                    rx = apply_cfo_phase(rx, freqOffsetNorm, phaseOffsetRad);

                    rx = correct_sync_with_preamble(rx, preambleSym, preambleLenHalf);
                end

                % --- Demodulazione ---
                rxBitsFrameHat = bpsk_demod(rx);

                % Se il frame è troppo corto (timing), penalizza
                if numel(rxBitsFrameHat) < numel(preambleBits) + numel(codedPayload)
                    errSum = errSum + numel(payload_bits);
                    bitSum = bitSum + numel(payload_bits);
                    continue;
                end

                % --- Estrai payload ---
                rxPayloadCodedHat = rxBitsFrameHat(numel(preambleBits)+1 : numel(preambleBits)+numel(codedPayload));

                % --- Decodifica ---
                rxPayloadHat = decode_bits(rxPayloadCodedHat, codingModes(ci), meta);

                L = min(numel(payload_bits), numel(rxPayloadHat));
                errSum = errSum + sum(payload_bits(1:L) ~= rxPayloadHat(1:L));
                bitSum = bitSum + L;
            end

            BER(ci,mi,si) = errSum / bitSum;
        end
    end
end

%% ===== 6) Demo ricostruzione testo (caso corretto) =====
testSNR = 12;
[demoCoded, metaDemo] = encode_bits(payload_bits, "HAMMING74");
txBitsFrame = [preambleBits, demoCoded];
txSym = bpsk_mod(txBitsFrame);

rx = channel_awgn(txSym, testSNR);
rx = apply_timing_offset(rx, timingDelay);
rx = apply_cfo_phase(rx, freqOffsetNorm, phaseOffsetRad);
rx = correct_sync_with_preamble(rx, preambleSym, preambleLenHalf);

rxBitsFrameHat = bpsk_demod(rx);
rxPayloadCodedHat = rxBitsFrameHat(numel(preambleBits)+1 : numel(preambleBits)+numel(demoCoded));
rxPayloadHat = decode_bits(rxPayloadCodedHat, "HAMMING74", metaDemo);

rx_text = bits_to_text(rxPayloadHat(1:numel(payload_bits)));

fprintf("TX: %s\n", tx_text);
fprintf("RX (corretto - AWGN_SYNC_CORR + Hamming): %s\n", rx_text);

%% ===== 7) Grafici =====
% (A) Confronto sincronizzazione per una codifica (consiglio: HAMMING74)
figure; hold on; grid on;
showCoding = "HAMMING74";
ciShow = find(codingModes == showCoding, 1);

plot_modes = ["AWGN","AWGN_SYNC","AWGN_SYNC_CORR"];
for i = 1:numel(plot_modes)
    miPlot = find(modes == plot_modes(i), 1);
    semilogy(SNRdB, squeeze(BER(ciShow,miPlot,:)), '-o', ...
        'DisplayName', sprintf('%s - %s', plot_modes(i), showCoding));
end
xlabel('SNR (dB)'); ylabel('BER');
title('Confronto effetti di sincronizzazione (con correzione)');
legend('Location','southwest');

% (B) Confronto codifiche su AWGN
figure; hold on; grid on;
miAWGN = find(modes == "AWGN", 1);
for ci2 = 1:numel(codingModes)
    semilogy(SNRdB, squeeze(BER(ci2,miAWGN,:)), '-o', ...
        'DisplayName', sprintf('AWGN - %s', codingModes(ci2)));
end
xlabel('SNR (dB)'); ylabel('BER');
title('Confronto codifiche (AWGN)');
legend('Location','southwest');

%% ===== 8) Salvataggio =====
saveas(figure(1), "BER_sync_comparison.png");
saveas(figure(2), "BER_coding_comparison.png");
save("results.mat","BER","SNRdB","modes","codingModes","Nframes", ...
    "timingDelay","freqOffsetNorm","phaseOffsetRad","loss_dB");

%% ==================== FUNZIONI ====================

function bits = text_to_bits(txt)
    bytes = uint8(char(txt));
    bits8 = zeros(numel(bytes), 8);
    for i = 1:numel(bytes)
        b = bytes(i);
        for k = 1:8
            bits8(i,k) = bitget(b, 9-k); % MSB -> LSB
        end
    end
    bits = reshape(bits8.', 1, []);
    bits = logical(bits);
end

function txt = bits_to_text(bits)
    bits = bits(:).';
    n = floor(numel(bits)/8)*8;
    bits = bits(1:n);
    bits8 = reshape(bits, 8, []).';
    bytes = zeros(size(bits8,1),1,'uint8');
    for i = 1:size(bits8,1)
        b = uint8(0);
        for k = 1:8
            if bits8(i,k)
                b = bitset(b, 9-k, 1);
            end
        end
        bytes(i) = b;
    end
    txt = char(bytes).';
end

function s = bpsk_mod(bits)
    % BPSK complesso (qui reale, ma complesso per gestire CFO/fase correttamente)
    s = complex(2*double(bits) - 1, 0);
end

function bits = bpsk_demod(r)
    % Decisione su parte reale
    bits = real(r) >= 0;
end

function y = channel_awgn(x, SNRdB)
    Ps = mean(abs(x).^2);
    SNR = 10^(SNRdB/10);
    N0 = Ps/SNR;
    n = sqrt(N0/2) * (randn(size(x)) + 1j*randn(size(x)));
    y = x + n;
end

function y = channel_rayleigh(x, SNRdB)
    h = (randn(size(x)) + 1j*randn(size(x))) / sqrt(2);
    faded = h .* x;

    Ps = mean(abs(faded).^2);
    SNR = 10^(SNRdB/10);
    N0 = Ps/SNR;
    n = sqrt(N0/2) * (randn(size(x)) + 1j*randn(size(x)));

    r = faded + n;
    y = r ./ h; % equalizzazione ideale (assunta)
end

function y = apply_timing_offset(x, delaySamples)
    % Applica un ritardo di "delaySamples" campioni (errore di sincronizzazione)
    if delaySamples <= 0
        y = x; return;
    end
    y = [zeros(1,delaySamples), x];
end

function y = apply_cfo_phase(x, freqOffsetNorm, phaseOffsetRad)
    % Applica offset di frequenza + fase: y[n] = x[n] * exp(j*(2*pi*fo*n + phi))
    n = 0:numel(x)-1;
    y = x .* exp(1j*(2*pi*freqOffsetNorm*n + phaseOffsetRad));
end

function rxCorr = correct_sync_with_preamble(rx, preambleSym, preambleLenHalf)
    % Correzione basata su preamble:
    % 1) Timing: correlazione con preamble noto
    % 2) CFO: stima tramite due metà uguali (fase media tra le metà)
    % 3) Compensazione CFO e fase residua

    Lpre = numel(preambleSym);

    % --- 1) Timing via correlazione ---
    c = xcorr(rx, preambleSym);
    [~, idx] = max(abs(c));
    start = idx - numel(rx);     % stima (approssimata) inizio
    start = max(1, start+1);

    rx2 = rx(start:end);
    if numel(rx2) < Lpre
        rxCorr = rx2; return;
    end

    preRx = rx2(1:Lpre);

    % --- 2) Stima CFO con due metà uguali ---
    A = preRx(1:preambleLenHalf);
    B = preRx(preambleLenHalf+1 : 2*preambleLenHalf);
    phi = angle(sum(conj(A).*B));             % fase media
    cfoEst = phi / (2*pi*preambleLenHalf);    % stima freqOffsetNorm

    % --- 3) Compensazione CFO su rx2 ---
    n = 0:numel(rx2)-1;
    rx3 = rx2 .* exp(-1j*(2*pi*cfoEst*n));

    % --- 4) Compensazione fase residua usando il preamble noto ---
    preRx3 = rx3(1:Lpre);
    phaseEst = angle(sum(conj(preambleSym).*preRx3));
    rxCorr = rx3 .* exp(-1j*phaseEst);
end

%% ===== Coding wrappers =====
function [outBits, meta] = encode_bits(inBits, mode)
    meta = struct();

    if mode == "NONE"
        outBits = logical(inBits);
        meta.pad = 0;

    elseif mode == "REP3"
        outBits = repelem(logical(inBits), 3);
        meta.pad = 0;

    else
        % Hamming(7,4): codifica a blocchi da 4 bit
        b = logical(inBits);
        pad = mod(4 - mod(numel(b),4), 4);
        if pad == 4, pad = 0; end
        b = [b, false(1,pad)];
        meta.pad = pad;

        outBits = false(1, (numel(b)/4)*7);
        p = 1;
        for i = 1:4:numel(b)
            d = b(i:i+3);
            c = hamming74_encode(d);
            outBits(p:p+6) = c;
            p = p + 7;
        end
    end
end

function outBits = decode_bits(inBits, mode, meta)
    if mode == "NONE"
        outBits = logical(inBits);

    elseif mode == "REP3"
        outBits = rep_decode(logical(inBits), 3);

    else
        % Hamming(7,4): decodifica e correzione di 1 errore per codeword
        b = logical(inBits);
        L = floor(numel(b)/7)*7;
        b = b(1:L);

        outBits = false(1, (L/7)*4);
        p = 1;
        for i = 1:7:numel(b)
            c = b(i:i+6);
            d = hamming74_decode(c);
            outBits(p:p+3) = d;
            p = p + 4;
        end

        if isfield(meta,'pad') && meta.pad > 0
            outBits = outBits(1:end-meta.pad);
        end
    end
end

function c = hamming74_encode(d)
    d1=d(1); d2=d(2); d3=d(3); d4=d(4);
    p1 = xor(xor(d1,d2), d4);
    p2 = xor(xor(d1,d3), d4);
    p3 = xor(xor(d2,d3), d4);
    c = [p1 p2 d1 p3 d2 d3 d4];
end

function d = hamming74_decode(c)
    p1=c(1); p2=c(2); d1=c(3); p3=c(4); d2=c(5); d3=c(6); d4=c(7);
    s1 = xor(xor(xor(p1,d1), d2), d4);
    s2 = xor(xor(xor(p2,d1), d3), d4);
    s3 = xor(xor(xor(p3,d2), d3), d4);
    syndrome = s1 + 2*s2 + 4*s3; % 1..7

    if syndrome ~= 0
        c(syndrome) = ~c(syndrome);
    end
    d = [c(3) c(5) c(6) c(7)];
end

function dec = rep_decode(bits, N)
    bits = bits(:).';
    L = floor(numel(bits)/N)*N;
    bits = bits(1:L);
    M = reshape(bits, N, []);
    dec = sum(M,1) >= (N/2);
    dec = logical(dec);
end
