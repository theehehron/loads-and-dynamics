clear; close all; clc;


acceptance_level = [  20, 0.0053
                     150, 0.04
                     600, 0.04
                    2000,0.0036];

qual_level = [acceptance_level(:, 1), acceptance_level(:, 2).*10^(6/10)];

acceptance_grms = grms(acceptance_level);
qual_grms = grms(qual_level);

loglog(acceptance_level(:, 1), acceptance_level(:,2))
hold on

freqs = 20:.1:2000;

%%
sdof_response_fn100Hz = sdof_PSD_response(acceptance_level, 100, 10, freqs);
sdof_response_fn200Hz = sdof_PSD_response(acceptance_level, 200, 10, freqs);
sdof_response_fn300Hz = sdof_PSD_response(acceptance_level, 300, 10, freqs);

loglog(freqs, sdof_response_fn100Hz, LineStyle='--')
loglog(freqs, sdof_response_fn200Hz, LineStyle='--')
loglog(freqs, sdof_response_fn300Hz, LineStyle='--')

legend("Base Input (6.14 gRMS)", "f_n = 100 Hz SDOF response", "f_n = 200 Hz SDOF response", "f_n = 300 Hz SDOF response", Location="NW")
ylim([0.001, 10])
grid on
title(["Response Power Sepctral Density Curves" ; "SDOF Systems  Q=10, Base Input = MIL-STD-1540C ATP"])
xlabel("Frequency (Hz)")
ylabel("Accel PSD (G^2/Hz)")

SDOF_300Hz_spectrum = [freqs', sdof_response_fn300Hz'];
SDOF_200Hz_spectrum = [freqs', sdof_response_fn200Hz'];
SDOF_100Hz_spectrum = [freqs', sdof_response_fn100Hz'];

grms_100Hz_sdof = grms(SDOF_100Hz_spectrum);
grms_200Hz_sdof = grms(SDOF_200Hz_spectrum);
grms_300Hz_sdof = grms(SDOF_300Hz_spectrum);

%%
figure();
vrs(acceptance_level);
title(["Vibration Response Spectrum"; "SDOF Systems, Q=10, Base Input = MIL-STD-1540C ATP"])