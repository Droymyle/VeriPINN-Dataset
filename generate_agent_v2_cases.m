%% generate_agent_v2_cases.m
% Agent 第二版评测用的特殊工况（单条仿真，流程与 generate_dataset_2.m 完全一致）：
%   训练范围外的参数（Lr 130%、Lm 70%、Ro 5 ohm、Cr 130%）
%   与训练条件不符的工况（Vin 280/340 V、fs 85 kHz；训练数据中 Vin、fs 固定）
% 输出 llc_agent_v2_cases.mat（字段同 llc_dataset_simulink_2.mat，另加 case_name）。
% 运行方式：matlab -batch "run('dataset/generate_agent_v2_cases.m')"

clear; clc;
mdl = 'LLC_simu';
load_system(fullfile(fileparts(mfilename('fullpath')), [mdl '.slx']));
stop_time = 0.02; transient_max_half = 3000; settle_band_frac = 0.02;
set_param(mdl, 'StopTime', num2str(stop_time));

Lr0 = 73e-6; Cr0 = 35e-9; Lm0 = 146e-6; Co0 = 2200e-6; n_nominal = 6.5; D = 45; Vpulse = 1;
% 名称, Vin, fs, Lr, Cr, Lm, Co, Ro
cases = {
  'lr130',  311, 90e3, 1.30*Lr0, 1.02*Cr0, 0.97*Lm0, 0.93*Co0, 2.4;
  'lm070',  311, 90e3, 0.98*Lr0, 0.96*Cr0, 0.70*Lm0, 1.05*Co0, 1.8;
  'ro5',    311, 90e3, 1.03*Lr0, 0.99*Cr0, 1.04*Lm0, 0.88*Co0, 5.0;
  'cr130',  311, 90e3, 0.97*Lr0, 1.30*Cr0, 1.02*Lm0, 0.95*Co0, 2.2;
  'vin280', 280, 90e3, 1.02*Lr0, 0.98*Cr0, 1.03*Lm0, 0.90*Co0, 2.0;
  'vin340', 340, 90e3, 0.96*Lr0, 1.03*Cr0, 0.99*Lm0, 1.04*Co0, 2.6;
  'fs85k',  311, 85e3, 1.01*Lr0, 1.01*Cr0, 0.98*Lm0, 0.86*Co0, 1.7;
};
num_sequences = size(cases, 1);
pts_per_period = 200; record_periods = 10; seq_len = record_periods * pts_per_period;

states_all = zeros(num_sequences, seq_len+1, 4);
v_sq_all = zeros(num_sequences, seq_len); vab_real_all = zeros(num_sequences, seq_len);
dt_half_all = zeros(num_sequences, 1);
vo_transient_all = zeros(num_sequences, transient_max_half+1); settling_time_all = zeros(num_sequences, 1);
vin_s = zeros(num_sequences,1); fs_s = vin_s; Lr_s = vin_s; Cr_s = vin_s; Lm_s = vin_s; Co_s = vin_s; Ro_s = vin_s;

t0 = tic;
for k = 1:num_sequences
    [Vin, fs_k, Lr, Cr, Lm, Co, Ro] = cases{k, 2:8};
    T = 1 / fs_k;
    vin_s(k) = Vin; fs_s(k) = fs_k; Lr_s(k) = Lr; Cr_s(k) = Cr; Lm_s(k) = Lm; Co_s(k) = Co; Ro_s(k) = Ro;

    simOut = sim(mdl);
    ILr = simOut.ILr; ILm = simOut.ILm; Vo = simOut.Vo; VCr = simOut.VCr; Vab = simOut.Vab;

    win_span = record_periods * T; half_period = T / 2;
    t_grid = linspace(stop_time - win_span, stop_time, seq_len + 1)';
    dt_half_all(k) = win_span / seq_len;
    states_all(k, :, :) = [interp1(ILr.Time, ILr.Data, t_grid), interp1(ILm.Time, ILm.Data, t_grid), ...
                           interp1(VCr.Time, VCr.Data, t_grid), interp1(Vo.Time, Vo.Data, t_grid)];
    vab_real_all(k, :) = interp1(Vab.Time, Vab.Data, t_grid(1:end-1));

    t_transient_grid = (0:transient_max_half)' * half_period;
    t_transient_grid = t_transient_grid(t_transient_grid <= Vo.Time(end));
    vo_transient_k = interp1(Vo.Time, Vo.Data, t_transient_grid);
    if numel(vo_transient_k) < transient_max_half + 1
        vo_transient_k(end+1:transient_max_half+1) = vo_transient_k(end);
    end
    tail_n = max(1, floor(numel(vo_transient_k) / 15));
    vo_ss_ref = mean(vo_transient_k(end-tail_n+1:end));
    within = abs(vo_transient_k - vo_ss_ref) <= settle_band_frac * abs(vo_ss_ref);
    settle_idx = numel(vo_transient_k);
    for j = 1:numel(within)
        if all(within(j:end)); settle_idx = j; break; end
    end
    vo_transient_all(k, :) = vo_transient_k';
    settling_time_all(k) = (settle_idx - 1) * half_period;

    n_half_from_end = floor((stop_time - t_grid(1:end-1)) / half_period + 1e-9);
    v_sq_all(k, :) = Vin * (-1).^n_half_from_end;
    fprintf('[%d/%d] %s 完成，累计耗时 %.1fs\n', k, num_sequences, cases{k,1}, toc(t0));
end

states = states_all; v_sq_ideal = v_sq_all; vab_real = vab_real_all;
rload = Ro_s; fs = fs_s; vin = vin_s; duty = D * ones(num_sequences, 1);
Lr = Lr_s; Cr = Cr_s; Lm = Lm_s; Co = Co_s; n = n_nominal * ones(num_sequences, 1);
dt_half = dt_half_all; vo_transient = vo_transient_all; settling_time = settling_time_all;
case_name = char(cases(:, 1));
out_path = fullfile(fileparts(mfilename('fullpath')), 'llc_agent_v2_cases.mat');
save(out_path, 'states', 'v_sq_ideal', 'vab_real', 'rload', 'fs', 'vin', 'duty', 'Lr', 'Cr', 'Lm', 'Co', ...
     'n', 'dt_half', 'vo_transient', 'settling_time', 'case_name', '-v7');
fprintf('已保存至 %s\n', out_path);
