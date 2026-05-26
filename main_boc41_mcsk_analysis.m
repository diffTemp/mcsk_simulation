%% main_boc41_mcsk_analysis.m
% BOC(4,1)-MCSK 复合信号生成与基础性能分析
%
% 论文《Signal Design and Compatibility Assessment for LEO Navigation
% Augmentation System》中明确给出的相关内容：
% 1) LEO 导航增强信号由高精度测量分量和高速数据分量组成；
% 2) 高速数据分量可采用 CSK，表 2 给出 CSK(4,1)；
% 3) 两个分量通过 MCSK，即 code-period time division multiplexing，形成复合信号；
% 4) BOC(4,1) 是论文考虑的候选调制之一；
% 5) 论文使用 SSC 和有效 C/N0 降低量进行兼容性评估；
% 6) 论文表 3 中 BOC(4,1) 的 SSC 参考值为：
%    BOC(4,1) vs BPSK(1):        -79.0205 dB
%    BOC(4,1) vs MBOC(6,1,1/11): -78.8849 dB
%    BOC(4,1) vs QPSK(10):       -71.8354 dB
%
% 本脚本中的复现仿真假设，论文并未明确给出：
% 1) PRN 码使用随机 +/-1 序列代替，长度为 1023 chips；
% 2) 测量分量低速导航数据 bit 暂设为 +1；
% 3) MCSK 采用 1:1 码周期交替结构，即 measurement/data 交替；
% 4) 仿真时长设为 20 ms；
% 5) 最终复合信号进行单位平均功率归一化；
% 6) PSD 采用本脚本手写 Welch 方法估计，Welch 参数是仿真设置；
% 7) MBOC(6,1,1/11) 参考 PSD 采用仿真近似加权方式生成。
% 8) 闭环仿真中的 AWGN、残余多普勒、初始相位、DLL/PLL 环路增益和
%    CSK 判决接收机均为工程近似设置，不是论文明确给出的接收机体制。
%
% 注意：
% 本脚本只复现 BOC(4,1)-MCSK 信号结构与 PSD、ACF、SSC、Gabor 带宽
% 分析流程。由于随机 PRN、有限时长、Welch 估计、MCSK 复合结构和参考
% PSD 近似均会影响结果，SSC 数值不应被解读为论文表 3 的精确复算。

clear; close all; clc;

%% 1. 初始化参数
fs = 81.84e6;              % 采样率，80 * 1.023 MHz
f0 = 1.023e6;              % 基准频率
fsc = 4 * f0;              % BOC(4,1) 副载波频率
Rc = 1 * f0;               % BOC(4,1) 码速率
codeLength = 1023;         % 1 ms 码周期内 chip 数
samplesPerChip = fs / Rc;  % 每 chip 采样点数
codePeriod = 1e-3;         % 码周期
numPeriods = 20;           % 总仿真时长 20 ms

assert(abs(samplesPerChip - round(samplesPerChip)) < 1e-12, ...
    'samplesPerChip must be an integer.');
samplesPerChip = round(samplesPerChip);
samplesPerCode = codeLength * samplesPerChip;

assert(samplesPerChip == 80, 'samplesPerChip must be 80.');
assert(samplesPerCode == 81840, 'samplesPerCode must be 81840.');
assert(abs(samplesPerCode / fs - codePeriod) < 1e-12, ...
    'samplesPerCode/fs must be equal to 1 ms.');

%% 2. 生成 PRN 码
rng(1);                                      % 保证结果可复现
prnBits = randi([0 1], 1, codeLength);       % 随机二进制序列
prnCode = 2 * prnBits - 1;                   % 0 -> -1, 1 -> +1

%% 3. 生成 BOC(4,1) 副载波
t1 = (0:samplesPerCode-1) / fs;              % 一个 1 ms 码周期内的时间轴
bocSubcarrier = sign(sin(2*pi*fsc*t1));      % BOC 方波副载波
bocSubcarrier(bocSubcarrier == 0) = 1;       % 避免 sign(0) 产生 0
assert(all(abs(bocSubcarrier) == 1), 'BOC subcarrier must be +/-1.');

%% 4. 生成测量分量
prnUpsampled = repelem(prnCode, samplesPerChip);
navBitMeasurement = 1;                       % 低速导航数据 bit 暂设为 +1
measurementSig = navBitMeasurement .* prnUpsampled .* bocSubcarrier;

%% 5. 生成高速数据分量
dataPeriodCount = floor(numPeriods / 2);
cskIndexList = zeros(1, dataPeriodCount);
cskBitList = zeros(dataPeriodCount, 4);
dataSigAll = zeros(dataPeriodCount, samplesPerCode);

for dataIdx = 1:dataPeriodCount
    % CSK(4,1)：每 4 bit 映射为一个码移索引 q，q 的单位是 chip。
    dataBits = randi([0 1], 1, 4);
    q = dataBits * [8; 4; 2; 1];             % MSB-first 二进制转十进制
    shiftedPrnCode = circshift(prnCode, [0 q]);
    shiftedPrnUpsampled = repelem(shiftedPrnCode, samplesPerChip);
    dataSig = shiftedPrnUpsampled .* bocSubcarrier;

    cskBitList(dataIdx, :) = dataBits;
    cskIndexList(dataIdx) = q;
    dataSigAll(dataIdx, :) = dataSig;
end

%% 6. 生成 MCSK 复合信号
s_mcsk = zeros(1, numPeriods * samplesPerCode);
periodType = cell(numPeriods, 1);
dataCounter = 0;

for periodIdx = 1:numPeriods
    sampleIdx = (periodIdx-1)*samplesPerCode + (1:samplesPerCode);

    if mod(periodIdx, 2) == 1
        s_mcsk(sampleIdx) = measurementSig;
        periodType{periodIdx} = 'measurement';
    else
        dataCounter = dataCounter + 1;
        s_mcsk(sampleIdx) = dataSigAll(dataCounter, :);
        periodType{periodIdx} = 'data';
    end
end

% 最终复合信号单位平均功率归一化。
s_mcsk = s_mcsk / sqrt(mean(abs(s_mcsk).^2));
t = (0:length(s_mcsk)-1) / fs;

%% 7. PSD 分析：手写 Welch 方法，不依赖通信工具箱或 pwelch
welchSegmentLength = 65536;
welchOverlapLength = welchSegmentLength / 2;
welchNfft = 262144;

[f_psd, Pxx] = manualWelchPsd(s_mcsk, fs, welchSegmentLength, ...
    welchOverlapLength, welchNfft);
df = f_psd(2) - f_psd(1);
Pxx_norm = Pxx / (sum(Pxx) * df);
psdPowerCheck = sum(Pxx_norm) * df;

%% 8. ACF 分析：FFT 方法计算有限 lag 自相关
maxLagChips = 10;
maxLagSamples = maxLagChips * samplesPerChip;
[acf, lags] = limitedAutocorrelationFft(s_mcsk, maxLagSamples);
acf = real(acf);
acf = acf / max(abs(acf));
lagChips = lags / samplesPerChip;

zeroLagIdx = find(lags == 0, 1);
mainPeakValue = acf(zeroLagIdx);

sideSearchMask = abs(lagChips) >= 0.1;
sideSearchValues = abs(acf);
sideSearchValues(~sideSearchMask) = -Inf;
[largestSidePeakValue, sidePeakIdx] = max(sideSearchValues);
largestSidePeakDelayChips = lagChips(sidePeakIdx);
largestSidePeakSignedValue = acf(sidePeakIdx);

%% 9. SSC 分析
% 论文表 3 中 BOC(4,1) 参考 SSC 数值：
% BOC(4,1) vs BPSK(1):        -79.0205 dB
% BOC(4,1) vs MBOC(6,1,1/11): -78.8849 dB
% BOC(4,1) vs QPSK(10):       -71.8354 dB
%
% 这里的 SSC 是仿真近似结果。计算时所有 PSD 均在相同 f_psd 网格上生成；
% 为保证指定接收带宽内的积分一致性，待测 PSD 和参考 PSD 都在对应
% beta_r 内重新归一化后再积分。
beta_L1 = 14.322e6;
beta_L5 = 20.46e6;

Rc_bpsk = 1.023e6;
G_bpsk_raw = normalizedSinc(f_psd / Rc_bpsk).^2;

Rc_qpsk = 10.23e6;
G_qpsk10_raw = normalizedSinc(f_psd / Rc_qpsk).^2;

% MBOC(6,1,1/11) 参考 PSD：用 BOC(1,1) 和 BOC(6,1) 仿真 PSD 加权近似。
s_boc11_ref = generateBocReferenceSignal(prnCode, fs, f0, 1, ...
    samplesPerChip, numPeriods);
s_boc61_ref = generateBocReferenceSignal(prnCode, fs, f0, 6, ...
    samplesPerChip, numPeriods);

[f_boc11, P_boc11] = manualWelchPsd(s_boc11_ref, fs, welchSegmentLength, ...
    welchOverlapLength, welchNfft);
[f_boc61, P_boc61] = manualWelchPsd(s_boc61_ref, fs, welchSegmentLength, ...
    welchOverlapLength, welchNfft);

assert(max(abs(f_boc11 - f_psd)) < 1e-9, 'BOC(1,1) PSD grid mismatch.');
assert(max(abs(f_boc61 - f_psd)) < 1e-9, 'BOC(6,1) PSD grid mismatch.');

P_boc11_norm = P_boc11 / (sum(P_boc11) * df);
P_boc61_norm = P_boc61 / (sum(P_boc61) * df);
G_mboc_raw = (10/11) * P_boc11_norm + (1/11) * P_boc61_norm;
G_mboc_raw = G_mboc_raw / (sum(G_mboc_raw) * df);

[ssc_bpsk, G_boc41_l1, G_bpsk_l1] = computeSsc(Pxx_norm, G_bpsk_raw, ...
    f_psd, beta_L1);
[ssc_mboc, ~, G_mboc_l1] = computeSsc(Pxx_norm, G_mboc_raw, ...
    f_psd, beta_L1);
[ssc_qpsk10, G_boc41_l5, G_qpsk10_l5] = computeSsc(Pxx_norm, G_qpsk10_raw, ...
    f_psd, beta_L5);

referenceSignal = {'BPSK(1)'; 'MBOC(6,1,1/11)'; 'QPSK(10)'};
bandwidthMHz = [beta_L1; beta_L1; beta_L5] / 1e6;
sscLinear = [ssc_bpsk; ssc_mboc; ssc_qpsk10];
sscDb = 10 * log10(sscLinear + eps);

sscTable = table(referenceSignal, bandwidthMHz, sscLinear, sscDb, ...
    'VariableNames', {'ReferenceSignal', 'Bandwidth_MHz', ...
    'SSC_linear', 'SSC_dB'});

%% 10. Gabor 带宽分析
beta_G_Hz = sqrt(sum((f_psd.^2) .* Pxx_norm) * df);
beta_G_MHz = beta_G_Hz / 1e6;

%% 11. 闭环接收机性能仿真：码跟踪误差、载波相位噪声、BER/FER
% 默认 quickClosedLoopTest=true，用于快速验证语法、维度和曲线生成。
% 若需要完整仿真，把 quickClosedLoopTest 改为 false，即使用 25:5:50 dB-Hz、
% 10 次 Monte Carlo、每次 5 帧；去掉 1 帧收敛期后，每个 C/N0 统计
% 40 frames 和 1600 bits。
quickClosedLoopTest = false;
closedLoopCNoDbHzFull = 0:1:50;
mcRunsPerCNoFull = 10;
framesPerRunFull = 5;
periodsPerFrame = 20;
settlingFrames = 1;

if quickClosedLoopTest
    closedLoopCNoDbHz = [35 45];
    mcRunsPerCNo = 2;
    framesPerRun = 2;
else
    closedLoopCNoDbHz = closedLoopCNoDbHzFull;
    mcRunsPerCNo = mcRunsPerCNoFull;
    framesPerRun = framesPerRunFull;
end

closedLoopParams = struct();
closedLoopParams.quickClosedLoopTest = quickClosedLoopTest;
closedLoopParams.closedLoopCNoDbHzFull = closedLoopCNoDbHzFull;
closedLoopParams.closedLoopCNoDbHz = closedLoopCNoDbHz;
closedLoopParams.mcRunsPerCNo = mcRunsPerCNo;
closedLoopParams.framesPerRun = framesPerRun;
closedLoopParams.periodsPerFrame = periodsPerFrame;
closedLoopParams.settlingFrames = settlingFrames;
closedLoopParams.trueCodeDelayChips = 0.35;
closedLoopParams.earlyLateSpacingChips = 0.1;
closedLoopParams.dllGainChips = 0.02;
closedLoopParams.codePeakSearchGridChips = -0.04:0.005:0.04;
closedLoopParams.codePeakAssistGain = 0.3;
closedLoopParams.codePeakOffsetPenalty = 0.15;
closedLoopParams.trueResidualDopplerHz = 25;
closedLoopParams.initialReceiverDopplerHz = 0;
closedLoopParams.pllPhaseGain = 0.35;
closedLoopParams.pllFreqGainHz = 12;
closedLoopParams.speedOfLight = 299792458;
closedLoopParams.sanityCNoDbHz = 55;

[closedLoopResults, exampleTrace, closedLoopSanity] = runClosedLoopSweep( ...
    prnCode, bocSubcarrier, fs, Rc, samplesPerChip, samplesPerCode, ...
    closedLoopParams);

codeTrackingRmseChips = closedLoopResults.CodeRMSE_chips;
codeTrackingRmseMeters = closedLoopResults.CodeRMSE_m;
carrierPhaseNoiseRad = closedLoopResults.CarrierPhaseStd_rad;
carrierPhaseNoiseDeg = closedLoopResults.CarrierPhaseStd_deg;
berList = closedLoopResults.BER;
ferList = closedLoopResults.FER;

assert(all(isfinite(codeTrackingRmseChips)), 'Code tracking RMSE contains invalid values.');
assert(all(isfinite(codeTrackingRmseMeters)), 'Code tracking RMSE in meters contains invalid values.');
assert(all(isfinite(carrierPhaseNoiseRad)), 'Carrier phase noise contains invalid values.');
assert(all(isfinite(berList)), 'BER contains invalid values.');
assert(all(isfinite(ferList)), 'FER contains invalid values.');
assert(isfinite(closedLoopSanity.CodeRMSE_chips), 'High C/N0 sanity code RMSE is invalid.');
assert(isfinite(closedLoopSanity.FinalCodeError_chips), ...
    'High C/N0 sanity final code error is invalid.');
assert(isfinite(closedLoopSanity.BER) && isfinite(closedLoopSanity.FER), ...
    'High C/N0 sanity BER/FER is invalid.');

%% 12. 绘图
% 图 1：BOC(4,1)-MCSK 复合信号前 20 us 的时域波形
figure('Name', 'Figure 1: MCSK time waveform');
showSamples20us = round(20e-6 * fs);
plot(t(1:showSamples20us) * 1e6, s_mcsk(1:showSamples20us), 'LineWidth', 1);
xlabel('Time (\mus)');
ylabel('Amplitude');
title('BOC(4,1)-MCSK Composite Signal: First 20 \mus');
grid on;

% 图 2：一个码周期内的测量分量局部波形
figure('Name', 'Figure 2: Measurement component');
localChipWindow = 8;
localIdx = 1:(localChipWindow * samplesPerChip);
plot((localIdx-1) / samplesPerChip, measurementSig(localIdx), 'LineWidth', 1);
xlabel('Chip index within code period');
ylabel('Amplitude');
title('Local Waveform of Measurement Component');
grid on;

% 图 3：一个码周期内的高速数据分量局部波形
figure('Name', 'Figure 3: Data component');
firstDataSig = dataSigAll(1, :);
plot((localIdx-1) / samplesPerChip, firstDataSig(localIdx), 'LineWidth', 1);
xlabel('Chip index within code period');
ylabel('Amplitude');
title(sprintf('Local Waveform of Data Component, CSK index q = %d', ...
    cskIndexList(1)));
grid on;

% 图 4：最终 MCSK 复合信号的 PSD
figure('Name', 'Figure 4: PSD');
plot(f_psd / 1e6, 10*log10(Pxx_norm + eps), 'LineWidth', 1);
xlabel('Frequency (MHz)');
ylabel('Normalized PSD (dB/Hz)');
title('PSD of BOC(4,1)-MCSK Composite Signal');
grid on;
xlim([-15 15]);
psdTextY = max(10*log10(Pxx_norm + eps)) - 6;
text(-14.5, psdTextY, sprintf('Gabor bandwidth = %.3f MHz', beta_G_MHz), ...
    'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7]);

% 图 5：最终 MCSK 复合信号的 ACF
figure('Name', 'Figure 5: ACF');
plot(lagChips, acf, 'LineWidth', 1);
hold on;
plot(0, mainPeakValue, 'ro', 'MarkerFaceColor', 'r');
plot(largestSidePeakDelayChips, largestSidePeakSignedValue, 'ks', ...
    'MarkerFaceColor', 'y');
text(0.15, mainPeakValue, sprintf('Main peak = %.3f', mainPeakValue));
text(largestSidePeakDelayChips, largestSidePeakSignedValue, ...
    sprintf(' Side peak |R| = %.3f at %.3f chips', ...
    largestSidePeakValue, largestSidePeakDelayChips));
xlabel('Delay (chips)');
ylabel('Normalized ACF');
title('ACF of BOC(4,1)-MCSK Composite Signal');
grid on;
xlim([-5 5]);

% 图 6：参考信号 PSD 对比图
figure('Name', 'Figure 6: Reference PSD comparison');
plot(f_psd / 1e6, 10*log10(Pxx_norm + eps), 'LineWidth', 1.2);
hold on;
plot(f_psd / 1e6, 10*log10(G_bpsk_l1 + eps), 'LineWidth', 1);
plot(f_psd / 1e6, 10*log10(G_mboc_l1 + eps), 'LineWidth', 1);
plot(f_psd / 1e6, 10*log10(G_qpsk10_l5 + eps), 'LineWidth', 1);
xlabel('Frequency (MHz)');
ylabel('Normalized PSD (dB/Hz)');
title('PSD Comparison with Reference GNSS Signals');
legend('BOC(4,1)-MCSK', 'BPSK(1)', 'MBOC(6,1,1/11)', ...
    'QPSK(10)', 'Location', 'best');
grid on;
xlim([-15 15]);

% 图 7：码跟踪 RMSE vs C/N0
figure('Name', 'Figure 7: Code tracking RMSE');
plot(closedLoopCNoDbHz, codeTrackingRmseChips, 'o-', 'LineWidth', 1.2);
xlabel('C/N_0 (dB-Hz)');
ylabel('Code tracking RMSE (chips)');
title('DLL Code Tracking Error of BOC(4,1)-MCSK');
grid on;

% 图 8：载波相位噪声 vs C/N0
figure('Name', 'Figure 8: Carrier phase noise');
plot(closedLoopCNoDbHz, carrierPhaseNoiseDeg, 's-', 'LineWidth', 1.2);
xlabel('C/N_0 (dB-Hz)');
ylabel('Carrier phase noise (deg)');
title('PLL Carrier Phase Noise of BOC(4,1)-MCSK');
grid on;

% 图 9：BER/FER vs C/N0
figure('Name', 'Figure 9: BER and FER');
berPlot = max(berList, 1 ./ max(2*closedLoopResults.BitCount, 1));
ferPlot = max(ferList, 1 ./ max(2*closedLoopResults.FrameCount, 1));
semilogy(closedLoopCNoDbHz, berPlot, 'o-', 'LineWidth', 1.2);
hold on;
semilogy(closedLoopCNoDbHz, ferPlot, 's-', 'LineWidth', 1.2);
xlabel('C/N_0 (dB-Hz)');
ylabel('Error rate');
title('CSK(4,1) BER/FER of BOC(4,1)-MCSK');
legend('BER', 'FER', 'Location', 'southwest');
grid on;

% 图 10：示例 C/N0 下 DLL/PLL 收敛轨迹
figure('Name', 'Figure 10: DLL/PLL convergence trace');
subplot(2, 1, 1);
plot(exampleTrace.periodIndex, exampleTrace.codeErrorChips, 'LineWidth', 1.1);
xlabel('Period index');
ylabel('Code error (chips)');
title(sprintf('DLL Trace after Acquisition, C/N_0 = %.1f dB-Hz', ...
    exampleTrace.CNoDbHz));
grid on;
subplot(2, 1, 2);
plot(exampleTrace.periodIndex, exampleTrace.phaseErrorDeg, 'LineWidth', 1.1);
xlabel('Period index');
ylabel('Phase error (deg)');
title('PLL Residual Phase Error Trace');
grid on;

%% 13. 命令行输出
fprintf('\n============================================================\n');
fprintf('BOC(4,1) 参数\n');
fprintf('============================================================\n');
fprintf('fsc = %.6f MHz\n', fsc / 1e6);
fprintf('Rc = %.6f Mcps\n', Rc / 1e6);
fprintf('codeLength = %d chips\n', codeLength);
fprintf('samplesPerChip = %d\n', samplesPerChip);
fprintf('samplesPerCode = %d\n', samplesPerCode);
fprintf('PSD normalization check, integral = %.12f\n', psdPowerCheck);

fprintf('\n============================================================\n');
fprintf('MCSK 结构\n');
fprintf('============================================================\n');
fprintf('numPeriods = %d\n', numPeriods);
fprintf('measurement 周期数 = %d\n', sum(strcmp(periodType, 'measurement')));
fprintf('data 周期数 = %d\n', sum(strcmp(periodType, 'data')));
fprintf('cskIndexList = ');
fprintf('%d ', cskIndexList);
fprintf('\n');

fprintf('\n============================================================\n');
fprintf('ACF 结果\n');
fprintf('============================================================\n');
fprintf('mainPeakValue = %.12f\n', mainPeakValue);
fprintf('largestSidePeakValue = %.12f\n', largestSidePeakValue);
fprintf('largestSidePeakDelayChips = %.6f\n', largestSidePeakDelayChips);

fprintf('\n============================================================\n');
fprintf('SSC 表格，仿真近似结果\n');
fprintf('============================================================\n');
fprintf('%-24s %-16s %-16s %-16s\n', ...
    'Reference Signal', 'Bandwidth(MHz)', 'SSC(linear)', 'SSC(dB)');
for rowIdx = 1:numel(referenceSignal)
    fprintf('%-24s %-16.3f %-16.6e %-16.6f\n', ...
        referenceSignal{rowIdx}, bandwidthMHz(rowIdx), ...
        sscLinear(rowIdx), sscDb(rowIdx));
end
fprintf('\nMATLAB table variable sscTable:\n');
disp(sscTable);

fprintf('\n============================================================\n');
fprintf('Gabor 带宽\n');
fprintf('============================================================\n');
fprintf('beta_G_Hz = %.6f Hz\n', beta_G_Hz);
fprintf('beta_G_MHz = %.6f MHz\n', beta_G_MHz);
fprintf('Gabor bandwidth of BOC(4,1)-MCSK = %.6f MHz\n', beta_G_MHz);

fprintf('\n============================================================\n');
fprintf('闭环仿真结果，工程近似接收机\n');
fprintf('============================================================\n');
fprintf('quickClosedLoopTest = %d\n', quickClosedLoopTest);
fprintf('C/N0 list = ');
fprintf('%.1f ', closedLoopCNoDbHz);
fprintf('dB-Hz\n');
fprintf('mcRunsPerCNo = %d, framesPerRun = %d, periodsPerFrame = %d, settlingFrames = %d\n', ...
    mcRunsPerCNo, framesPerRun, periodsPerFrame, settlingFrames);
fprintf('trueCodeDelayChips = %.3f，仅用于信道生成和离线误差统计\n', ...
    closedLoopParams.trueCodeDelayChips);
fprintf('%-12s %-18s %-18s %-18s %-18s %-18s %-14s %-14s\n', ...
    'C/N0(dB-Hz)', 'AcqRMSE(chips)', 'CodeRMSE(chips)', 'FinalErr(chips)', ...
    'CarrierStd(rad)', 'CarrierStd(deg)', 'BER', 'FER');
for rowIdx = 1:height(closedLoopResults)
    fprintf('%-12.1f %-18.6e %-18.6e %-18.6e %-18.6e %-18.6f %-14.6e %-14.6e\n', ...
        closedLoopResults.CNo_dBHz(rowIdx), ...
        closedLoopResults.AcquisitionRMSE_chips(rowIdx), ...
        closedLoopResults.CodeRMSE_chips(rowIdx), ...
        closedLoopResults.FinalCodeError_chips(rowIdx), ...
        closedLoopResults.CarrierPhaseStd_rad(rowIdx), ...
        closedLoopResults.CarrierPhaseStd_deg(rowIdx), ...
        closedLoopResults.BER(rowIdx), ...
        closedLoopResults.FER(rowIdx));
end
fprintf('\nMATLAB table variable closedLoopResults:\n');
disp(closedLoopResults);
fprintf('High C/N0 sanity check at %.1f dB-Hz: AcqRMSE = %.6e chips, CodeRMSE = %.6e chips, FinalCodeError = %.6e chips, BER = %.6e, FER = %.6e\n', ...
    closedLoopParams.sanityCNoDbHz, closedLoopSanity.AcquisitionRMSE_chips, ...
    closedLoopSanity.CodeRMSE_chips, closedLoopSanity.FinalCodeError_chips, ...
    closedLoopSanity.BER, closedLoopSanity.FER);
fprintf('Sanity final estimate at %.1f dB-Hz: acquiredCodeDelay = %.6f chips, finalCodeDelayEstimate = %.6f chips\n', ...
    closedLoopParams.sanityCNoDbHz, closedLoopSanity.AcquiredCodeDelay_chips, ...
    closedLoopSanity.FinalCodeDelayEstimate_chips);

%% 本地函数
function [f, Pxx] = manualWelchPsd(x, fs, segmentLength, overlapLength, nfft)
% manualWelchPsd 使用基础 FFT 实现双边 Welch PSD 估计。
% 输出频率 f 已经 fftshift 到 baseband，单位为 Hz。
    x = x(:);
    segmentLength = round(segmentLength);
    overlapLength = round(overlapLength);
    nfft = round(nfft);

    assert(segmentLength > 0, 'segmentLength must be positive.');
    assert(overlapLength >= 0 && overlapLength < segmentLength, ...
        'overlapLength must be in [0, segmentLength).');
    assert(nfft >= segmentLength, 'nfft must be >= segmentLength.');
    assert(length(x) >= segmentLength, 'Signal is shorter than one segment.');

    step = segmentLength - overlapLength;
    numSegments = floor((length(x) - overlapLength) / step);

    n = (0:segmentLength-1).';
    window = 0.5 - 0.5*cos(2*pi*n/(segmentLength-1)); % Hann window
    windowPower = sum(window.^2);

    PxxAccum = zeros(nfft, 1);
    for segIdx = 1:numSegments
        startIdx = (segIdx-1)*step + 1;
        stopIdx = startIdx + segmentLength - 1;
        segment = x(startIdx:stopIdx) .* window;
        X = fft(segment, nfft);
        PxxAccum = PxxAccum + (abs(X).^2) / (fs * windowPower);
    end

    Pxx = PxxAccum / numSegments;
    Pxx = fftshift(Pxx);
    f = (-nfft/2:nfft/2-1).' * (fs / nfft);
end

function [acf, lags] = limitedAutocorrelationFft(x, maxLagSamples)
% limitedAutocorrelationFft 使用 FFT 计算 biased 自相关，并只返回有限 lag。
    x = x(:);
    n = length(x);
    nfft = 2^nextpow2(2*n - 1);
    X = fft(x, nfft);
    r = ifft(abs(X).^2);
    r = r(1:maxLagSamples+1) / n;

    acf = [conj(flipud(r(2:end))); r];
    lags = (-maxLagSamples:maxLagSamples).';
end

function y = normalizedSinc(x)
% normalizedSinc 实现 MATLAB sinc 定义：sin(pi*x)/(pi*x)。
    y = ones(size(x));
    nonzeroMask = abs(x) > 1e-14;
    y(nonzeroMask) = sin(pi*x(nonzeroMask)) ./ (pi*x(nonzeroMask));
end

function s_ref = generateBocReferenceSignal(prnCode, fs, f0, m, ...
    samplesPerChip, numPeriods)
% generateBocReferenceSignal 生成 measurement-only BOC(m,1) 参考信号。
% 该函数仅用于 MBOC 参考 PSD 的仿真近似，不代表论文给出的接收机实现。
    codeLength = numel(prnCode);
    samplesPerCode = codeLength * samplesPerChip;
    t1 = (0:samplesPerCode-1) / fs;
    bocSubcarrier = sign(sin(2*pi*(m*f0)*t1));
    bocSubcarrier(bocSubcarrier == 0) = 1;

    prnUpsampled = repelem(prnCode, samplesPerChip);
    oneCodeSignal = prnUpsampled .* bocSubcarrier;
    s_ref = repmat(oneCodeSignal, 1, numPeriods);
    s_ref = s_ref / sqrt(mean(abs(s_ref).^2));
end

function [sscValue, G_target_band, G_ref_band] = computeSsc( ...
    targetPsd, referencePsd, f, beta)
% computeSsc 在指定接收带宽 beta 内归一化 PSD 并计算 SSC。
% beta 是双边接收带宽，因此积分区间为 [-beta/2, beta/2]。
    df = f(2) - f(1);
    bandMask = abs(f) <= beta/2;

    G_target_band = targetPsd / (sum(targetPsd(bandMask)) * df);
    G_ref_band = referencePsd / (sum(referencePsd(bandMask)) * df);
    sscValue = sum(G_target_band(bandMask) .* G_ref_band(bandMask)) * df;
end

function [resultsTable, exampleTrace, sanityResult] = runClosedLoopSweep( ...
    prnCode, bocSubcarrier, fs, Rc, samplesPerChip, samplesPerCode, params)
% runClosedLoopSweep 执行工程闭环接收机仿真。
% 该闭环模型用于观察趋势：measurement 周期更新 DLL/PLL，data 周期用
% 16 路 CSK 相关器判决 q。它不是论文给出的具体接收机实现。
    rng(20260525);

    codeLength = numel(prnCode);
    codePeriod = samplesPerCode / fs;
    timeOneCode = (0:samplesPerCode-1) / fs;

    measurementReplica = repelem(prnCode, samplesPerChip) .* bocSubcarrier;
    dataCandidateMatrix = zeros(16, samplesPerCode);
    cskBitMap = zeros(16, 4);
    for q = 0:15
        shiftedPrnCode = circshift(prnCode, [0 q]);
        dataCandidateMatrix(q+1, :) = repelem(shiftedPrnCode, samplesPerChip) .* bocSubcarrier;
        cskBitMap(q+1, :) = indexToBits(q, 4);
    end

    trueDelaySamples = params.trueCodeDelayChips * samplesPerChip;
    measurementTxDelayed = fractionalDelayCircular(measurementReplica, trueDelaySamples);
    dataTxDelayedMatrix = zeros(size(dataCandidateMatrix));
    for q = 0:15
        dataTxDelayedMatrix(q+1, :) = fractionalDelayCircular( ...
            dataCandidateMatrix(q+1, :), trueDelaySamples);
    end

    dllSign = calibrateDllSign(measurementReplica, samplesPerChip, ...
        params.earlyLateSpacingChips);

    cnoList = params.closedLoopCNoDbHz(:);
    codeTrackingRmseChips = zeros(numel(cnoList), 1);
    codeTrackingRmseMeters = zeros(numel(cnoList), 1);
    carrierPhaseNoiseRad = zeros(numel(cnoList), 1);
    carrierPhaseNoiseDeg = zeros(numel(cnoList), 1);
    berList = zeros(numel(cnoList), 1);
    ferList = zeros(numel(cnoList), 1);
    acquiredCodeDelayChips = zeros(numel(cnoList), 1);
    acquisitionRmseChips = zeros(numel(cnoList), 1);
    finalCodeDelayEstimateChips = zeros(numel(cnoList), 1);
    finalCodeErrorChips = zeros(numel(cnoList), 1);
    bitCountList = zeros(numel(cnoList), 1);
    frameCountList = zeros(numel(cnoList), 1);

    exampleTrace = struct('CNoDbHz', cnoList(end), 'periodIndex', [], ...
        'codeErrorChips', [], 'phaseErrorDeg', [], ...
        'acquiredCodeDelayChips', NaN, 'acquisitionErrorChips', NaN);

    for cnoIdx = 1:numel(cnoList)
        keepTrace = (cnoIdx == numel(cnoList));
        [result, trace] = simulateClosedLoopCNo(cnoList(cnoIdx), ...
            params.mcRunsPerCNo, params.framesPerRun, params, fs, Rc, ...
            samplesPerChip, samplesPerCode, codePeriod, timeOneCode, ...
            measurementReplica, dataCandidateMatrix, measurementTxDelayed, ...
            dataTxDelayedMatrix, cskBitMap, dllSign, keepTrace);

        codeTrackingRmseChips(cnoIdx) = result.CodeRMSE_chips;
        codeTrackingRmseMeters(cnoIdx) = result.CodeRMSE_m;
        carrierPhaseNoiseRad(cnoIdx) = result.CarrierPhaseStd_rad;
        carrierPhaseNoiseDeg(cnoIdx) = result.CarrierPhaseStd_deg;
        berList(cnoIdx) = result.BER;
        ferList(cnoIdx) = result.FER;
        acquiredCodeDelayChips(cnoIdx) = result.AcquiredCodeDelay_chips;
        acquisitionRmseChips(cnoIdx) = result.AcquisitionRMSE_chips;
        finalCodeDelayEstimateChips(cnoIdx) = result.FinalCodeDelayEstimate_chips;
        finalCodeErrorChips(cnoIdx) = result.FinalCodeError_chips;
        bitCountList(cnoIdx) = result.BitCount;
        frameCountList(cnoIdx) = result.FrameCount;

        if keepTrace
            exampleTrace = trace;
        end
    end

    resultsTable = table(cnoList, codeTrackingRmseChips, codeTrackingRmseMeters, ...
        carrierPhaseNoiseRad, carrierPhaseNoiseDeg, berList, ferList, ...
        acquiredCodeDelayChips, acquisitionRmseChips, ...
        finalCodeDelayEstimateChips, finalCodeErrorChips, ...
        bitCountList, frameCountList, ...
        'VariableNames', {'CNo_dBHz', 'CodeRMSE_chips', 'CodeRMSE_m', ...
        'CarrierPhaseStd_rad', 'CarrierPhaseStd_deg', 'BER', 'FER', ...
        'AcquiredCodeDelay_chips', 'AcquisitionRMSE_chips', ...
        'FinalCodeDelayEstimate_chips', 'FinalCodeError_chips', ...
        'BitCount', 'FrameCount'});

    [sanityResult, ~] = simulateClosedLoopCNo(params.sanityCNoDbHz, ...
        1, max(params.settlingFrames + 2, 3), params, fs, Rc, ...
        samplesPerChip, samplesPerCode, codePeriod, timeOneCode, ...
        measurementReplica, dataCandidateMatrix, measurementTxDelayed, ...
        dataTxDelayedMatrix, cskBitMap, dllSign, false);
end

function [result, trace] = simulateClosedLoopCNo(cnoDbHz, mcRuns, framesPerRun, ...
    params, fs, Rc, samplesPerChip, samplesPerCode, codePeriod, timeOneCode, ...
    measurementReplica, dataCandidateMatrix, measurementTxDelayed, ...
    dataTxDelayedMatrix, cskBitMap, dllSign, keepTrace)
% simulateClosedLoopCNo 在一个 C/N0 点上统计 DLL、PLL 和 CSK 解调性能。
    cnoLinear = 10^(cnoDbHz / 10);
    noiseVariance = fs / cnoLinear;
    chipLengthMeters = params.speedOfLight / Rc;
    codeLength = samplesPerCode / samplesPerChip;

    codeErrorList = [];
    phaseErrorList = [];
    acquiredCodeDelayList = [];
    acquisitionErrorList = [];
    finalCodeDelayEstimateList = [];
    finalCodeErrorList = [];
    bitErrorCount = 0;
    bitTotalCount = 0;
    frameErrorCount = 0;
    frameTotalCount = 0;

    trace = struct('CNoDbHz', cnoDbHz, 'periodIndex', [], ...
        'codeErrorChips', [], 'phaseErrorDeg', [], ...
        'acquiredCodeDelayChips', NaN, 'acquisitionErrorChips', NaN);

    for runIdx = 1:mcRuns
        trueDopplerHz = params.trueResidualDopplerHz;
        estDopplerHz = params.initialReceiverDopplerHz;
        trueCarrierPhaseRad = wrapToPiLocal(2*pi*rand());
        estCarrierPhaseRad = 0;

        % Acquisition 使用接收机能看到的第一个 measurement 周期相关峰初始化码相位。
        % trueCodeDelayChips 只在下面的离线误差统计里使用，不反馈给接收机。
        acquisitionCarrier = exp(1j * (trueCarrierPhaseRad + ...
            2*pi*trueDopplerHz*timeOneCode));
        acquisitionNoise = sqrt(noiseVariance/2) * ...
            (randn(1, samplesPerCode) + 1j*randn(1, samplesPerCode));
        acquisitionRx = measurementTxDelayed .* acquisitionCarrier + acquisitionNoise;
        estCodeDelayChips = acquireCodeDelayChips(acquisitionRx, ...
            measurementReplica, samplesPerChip);
        acquisitionErrorChips = codeDelayDifferenceChips(estCodeDelayChips, ...
            params.trueCodeDelayChips, codeLength);

        acquiredCodeDelayList(end+1, 1) = estCodeDelayChips; %#ok<AGROW>
        acquisitionErrorList(end+1, 1) = acquisitionErrorChips; %#ok<AGROW>

        if keepTrace && runIdx == 1
            trace.acquiredCodeDelayChips = estCodeDelayChips;
            trace.acquisitionErrorChips = acquisitionErrorChips;
        end

        trueCarrierPhaseRad = wrapToPiLocal(trueCarrierPhaseRad + ...
            2*pi*trueDopplerHz*codePeriod);

        tracePeriod = 0;
        for frameIdx = 1:framesPerRun
            frameHasError = false;
            effectiveFrame = frameIdx > params.settlingFrames;

            for periodIdx = 1:params.periodsPerFrame
                tracePeriod = tracePeriod + 1;
                isMeasurementPeriod = mod(periodIdx, 2) == 1;

                if isMeasurementPeriod
                    txPeriod = measurementTxDelayed;
                    trueBits = [];
                else
                    trueBits = randi([0 1], 1, 4);
                    qTrue = bitsToIndex(trueBits);
                    txPeriod = dataTxDelayedMatrix(qTrue+1, :);
                end

                trueCarrier = exp(1j * (trueCarrierPhaseRad + ...
                    2*pi*trueDopplerHz*timeOneCode));
                noise = sqrt(noiseVariance/2) * ...
                    (randn(1, samplesPerCode) + 1j*randn(1, samplesPerCode));
                rxPeriod = txPeriod .* trueCarrier + noise;

                localCarrier = exp(-1j * (estCarrierPhaseRad + ...
                    2*pi*estDopplerHz*timeOneCode));
                rxBaseband = rxPeriod .* localCarrier;

                if isMeasurementPeriod
                    promptAligned = fractionalDelayCircular(rxBaseband, ...
                        -estCodeDelayChips * samplesPerChip);
                    promptCorr = sum(promptAligned .* measurementReplica) / samplesPerCode;

                    phaseErrorRad = wrapToPiLocal(atan2(imag(promptCorr), real(promptCorr)));
                    dllDisc = computeDllDiscriminator(rxBaseband, measurementReplica, ...
                        estCodeDelayChips, params.earlyLateSpacingChips, samplesPerChip);
                    peakOffsetChips = estimateCodePeakOffset(rxBaseband, measurementReplica, ...
                        estCodeDelayChips, params.codePeakSearchGridChips, samplesPerChip, ...
                        params.codePeakOffsetPenalty);

                    estCodeDelayChips = estCodeDelayChips + ...
                        dllSign * params.dllGainChips * dllDisc + ...
                        params.codePeakAssistGain * peakOffsetChips;
                    estCodeDelayChips = wrapCodeDelayChips(estCodeDelayChips, codeLength);
                    currentCodeErrorChips = codeDelayDifferenceChips(estCodeDelayChips, ...
                        params.trueCodeDelayChips, codeLength);

                    estCarrierPhaseRad = wrapToPiLocal(estCarrierPhaseRad + ...
                        params.pllPhaseGain * phaseErrorRad);
                    estDopplerHz = estDopplerHz + params.pllFreqGainHz * phaseErrorRad;

                    if effectiveFrame
                        codeErrorList(end+1, 1) = currentCodeErrorChips; %#ok<AGROW>
                        phaseErrorList(end+1, 1) = phaseErrorRad; %#ok<AGROW>
                    end

                    if keepTrace && runIdx == 1
                        trace.periodIndex(end+1, 1) = tracePeriod; %#ok<AGROW>
                        trace.codeErrorChips(end+1, 1) = currentCodeErrorChips; %#ok<AGROW>
                        trace.phaseErrorDeg(end+1, 1) = phaseErrorRad * 180/pi; %#ok<AGROW>
                    end
                else
                    promptAligned = fractionalDelayCircular(rxBaseband, ...
                        -estCodeDelayChips * samplesPerChip);
                    cskCorr = dataCandidateMatrix * promptAligned(:) / samplesPerCode;
                    [~, qHatIdx] = max(abs(cskCorr));
                    qHat = qHatIdx - 1;
                    bitsHat = cskBitMap(qHat+1, :);

                    if effectiveFrame
                        bitErrors = sum(bitsHat ~= trueBits);
                        bitErrorCount = bitErrorCount + bitErrors;
                        bitTotalCount = bitTotalCount + numel(trueBits);
                        frameHasError = frameHasError || (bitErrors > 0);
                    end
                end

                trueCarrierPhaseRad = wrapToPiLocal(trueCarrierPhaseRad + ...
                    2*pi*trueDopplerHz*codePeriod);
                estCarrierPhaseRad = wrapToPiLocal(estCarrierPhaseRad + ...
                    2*pi*estDopplerHz*codePeriod);
            end

            if effectiveFrame
                frameTotalCount = frameTotalCount + 1;
                frameErrorCount = frameErrorCount + double(frameHasError);
            end
        end

        finalCodeDelayEstimateList(end+1, 1) = estCodeDelayChips; %#ok<AGROW>
        finalCodeErrorList(end+1, 1) = codeDelayDifferenceChips(estCodeDelayChips, ...
            params.trueCodeDelayChips, codeLength); %#ok<AGROW>
    end

    if isempty(codeErrorList)
        codeErrorList = NaN;
    end
    if isempty(phaseErrorList)
        phaseErrorList = NaN;
    end

    result = struct();
    result.CNo_dBHz = cnoDbHz;
    result.CodeRMSE_chips = sqrt(mean(codeErrorList.^2, 'omitnan'));
    result.CodeRMSE_m = result.CodeRMSE_chips * chipLengthMeters;
    result.AcquiredCodeDelay_chips = mean(acquiredCodeDelayList, 'omitnan');
    result.AcquisitionRMSE_chips = sqrt(mean(acquisitionErrorList.^2, 'omitnan'));
    result.FinalCodeDelayEstimate_chips = mean(finalCodeDelayEstimateList, 'omitnan');
    result.FinalCodeError_chips = sqrt(mean(finalCodeErrorList.^2, 'omitnan'));
    result.CarrierPhaseStd_rad = std(phaseErrorList, 'omitnan');
    result.CarrierPhaseStd_deg = result.CarrierPhaseStd_rad * 180/pi;
    result.BER = bitErrorCount / max(bitTotalCount, 1);
    result.FER = frameErrorCount / max(frameTotalCount, 1);
    result.BitCount = bitTotalCount;
    result.FrameCount = frameTotalCount;
end

function dllSign = calibrateDllSign(measurementReplica, samplesPerChip, earlyLateSpacingChips)
% calibrateDllSign 用无噪声小偏差自动判断 DLL 更新方向。
    testErrorChips = 0.02;
    rxNoiseless = measurementReplica;
    testDisc = computeDllDiscriminator(rxNoiseless, measurementReplica, ...
        testErrorChips, earlyLateSpacingChips, samplesPerChip);

    if abs(testDisc) < 1e-12
        dllSign = 1;
    else
        dllSign = -sign(testErrorChips * testDisc);
    end
end

function acquiredDelayChips = acquireCodeDelayChips(rxMeasurement, ...
    measurementReplica, samplesPerChip)
% acquireCodeDelayChips 用循环相关峰估计接收机初始码相位。
% 该函数只使用接收样本和本地 measurement replica，不知道真实码延迟。
    rxMeasurement = rxMeasurement(:).';
    measurementReplica = measurementReplica(:).';
    n = numel(measurementReplica);

    corrMetric = abs(ifft(fft(rxMeasurement) .* conj(fft(measurementReplica))));
    [~, peakIdx] = max(corrMetric);

    leftIdx = mod(peakIdx - 2, n) + 1;
    rightIdx = mod(peakIdx, n) + 1;
    leftValue = corrMetric(leftIdx);
    centerValue = corrMetric(peakIdx);
    rightValue = corrMetric(rightIdx);
    denom = leftValue - 2*centerValue + rightValue;

    if abs(denom) > eps
        fracOffset = 0.5 * (leftValue - rightValue) / denom;
        fracOffset = min(max(fracOffset, -0.5), 0.5);
    else
        fracOffset = 0;
    end

    delaySamples = mod((peakIdx - 1) + fracOffset, n);
    acquiredDelayChips = delaySamples / samplesPerChip;
end

function wrappedDelayChips = wrapCodeDelayChips(delayChips, codeLength)
% wrapCodeDelayChips 将接收机内部码相位限制在一个码周期内。
% 这不是用真值限幅，只是码相位 NCO 的自然周期回绕。
    wrappedDelayChips = mod(delayChips, codeLength);
end

function diffChips = codeDelayDifferenceChips(estimateChips, truthChips, codeLength)
% codeDelayDifferenceChips 仅用于离线评分，返回最短循环码相位误差。
    diffChips = mod(estimateChips - truthChips + codeLength/2, codeLength) - codeLength/2;
end

function discriminator = computeDllDiscriminator(rxBaseband, measurementReplica, ...
    estCodeDelayChips, earlyLateSpacingChips, samplesPerChip)
% computeDllDiscriminator 使用 early-minus-late envelope 鉴别器。
    halfSpacing = earlyLateSpacingChips / 2;
    earlyAligned = fractionalDelayCircular(rxBaseband, ...
        -(estCodeDelayChips - halfSpacing) * samplesPerChip);
    lateAligned = fractionalDelayCircular(rxBaseband, ...
        -(estCodeDelayChips + halfSpacing) * samplesPerChip);

    earlyCorr = sum(earlyAligned .* measurementReplica) / numel(measurementReplica);
    lateCorr = sum(lateAligned .* measurementReplica) / numel(measurementReplica);
    discriminator = (abs(earlyCorr) - abs(lateCorr)) / ...
        (abs(earlyCorr) + abs(lateCorr) + eps);
end

function peakOffsetChips = estimateCodePeakOffset(rxBaseband, measurementReplica, ...
    estCodeDelayChips, searchGridChips, samplesPerChip, offsetPenalty)
% estimateCodePeakOffset 在当前码相位附近做窄范围 prompt 相关峰搜索。
% 该辅助校正用于避免 BOC(4,1) 多峰 ACF 下的短时错峰，不替代 DLL 输出。
    metric = zeros(size(searchGridChips));
    for gridIdx = 1:numel(searchGridChips)
        trialDelayChips = estCodeDelayChips + searchGridChips(gridIdx);
        trialAligned = fractionalDelayCircular(rxBaseband, ...
            -trialDelayChips * samplesPerChip);
        trialCorr = sum(trialAligned .* measurementReplica) / numel(measurementReplica);
        metric(gridIdx) = abs(trialCorr);
    end

    score = metric - offsetPenalty * abs(searchGridChips);
    [~, bestIdx] = max(score);
    peakOffsetChips = searchGridChips(bestIdx);
end

function y = fractionalDelayCircular(x, delaySamples)
% fractionalDelayCircular 对一个码周期信号做循环小数采样延迟。
% 正 delaySamples 表示输出相对输入向右延迟。
    rowInput = isrow(x);
    x = x(:).';
    n = numel(x);
    sampleIndex0 = 0:n-1;
    sourceIndex0 = mod(sampleIndex0 - delaySamples, n);
    indexFloor0 = floor(sourceIndex0);
    frac = sourceIndex0 - indexFloor0;
    indexNext0 = mod(indexFloor0 + 1, n);

    y = (1 - frac) .* x(indexFloor0 + 1) + frac .* x(indexNext0 + 1);
    if ~rowInput
        y = y(:);
    end
end

function q = bitsToIndex(bits)
% bitsToIndex 按 MSB-first 把 4 bit 转为 CSK 码移索引。
    weights = 2.^(numel(bits)-1:-1:0);
    q = bits * weights(:);
end

function bits = indexToBits(q, bitCount)
% indexToBits 按 MSB-first 把 CSK 码移索引转为 bit。
    bits = zeros(1, bitCount);
    for bitIdx = 1:bitCount
        bitWeight = 2^(bitCount - bitIdx);
        bits(bitIdx) = floor(q / bitWeight);
        q = q - bits(bitIdx) * bitWeight;
    end
end

function angleWrapped = wrapToPiLocal(angleInput)
% wrapToPiLocal 避免依赖 Mapping Toolbox。
    angleWrapped = mod(angleInput + pi, 2*pi) - pi;
end
