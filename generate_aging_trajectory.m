%% generate_aging_trajectory.m
% 输出电解电容老化轨迹：Lr/Cr/Lm 固定为铭牌值，只让 Co 从 100% 逐步降到 70%，
% 在几个负载下各跑一遍 LLC_simu.slx。仿真/采样/保存流程与 generate_dataset_2.m
% 完全一致（同一 stop_time、窗口、冷启动曲线算法），只是参数不随机抽取。
% Co 低于 80% 的点超出训练范围（±20%），用于检验外推。
%
% 输出 llc_aging_trajectory.mat（字段同 llc_dataset_simulink_2.mat，另加 co_ratio）。
% 运行方式：MATLAB 里直接 Run，或 matlab -batch "run('dataset/generate_aging_trajectory.m')"
% 只跑一条试算耗时：先在工作区设 smoke_test = true 再 run。

if ~exist('smoke_test', 'var'); smoke_test = false; end
clearvars -except smoke_test; clc;
mdl = 'LLC_simu';
load_system(fullfile(fileparts(mfilename('fullpath')), [mdl '.slx']));

stop_time = 0.02;
transient_max_half = 3000;
settle_band_frac = 0.02;
set_param(mdl, 'StopTime', num2str(stop_time));

Vin_nominal = 311; Lr_nominal = 73e-6; Cr_nominal = 35e-9; Lm_nominal = 146e-6;
Co_nominal = 2200e-6; n_nominal = 6.5; D = 45; Vpulse = 1; FS_nominal = 90e3;

co_ratio_list = (1.00:-0.025:0.70)';          % 13 个老化等级
ro_list       = [1.5; 2.0; 3.0];              % 负载 (ohm)
[co_grid, ro_grid] = ndgrid(co_ratio_list, ro_list);
co_ratio      = co_grid(:);
rload_samples = ro_grid(:);
if smoke_test
    co_ratio = co_ratio(1); rload_samples = rload_samples(1);
end
num_sequences = numel(co_ratio);

pts_per_period = 200; record_periods = 10;
seq_len = record_periods * pts_per_period;

states_all = zeros(num_sequences, seq_len+1, 4);
v_sq_all = zeros(num_sequences, seq_len);
vab_real_all = zeros(num_sequences, seq_len);
dt_half_all = zeros(num_sequences, 1);
vo_transient_all = zeros(num_sequences, transient_max_half+1);
settling_time_all = zeros(num_sequences, 1);

t0 = tic;
for k = 1:num_sequences
    fs_k = FS_nominal;
    Vin = Vin_nominal; Lr = Lr_nominal; Cr = Cr_nominal; Lm = Lm_nominal;
    Co  = co_ratio(k) * Co_nominal;
    Ro  = rload_samples(k);
    T   = 1 / fs_k;

    simOut = sim(mdl);
    ILr = simOut.ILr; ILm = simOut.ILm; Vo = simOut.Vo; VCr = simOut.VCr; Vab = simOut.Vab;

    win_span = record_periods * T;
    half_period = T / 2;
    t_grid = linspace(stop_time - win_span, stop_time, seq_len + 1)';
    dt_half_all(k) = win_span / seq_len;

    iLr_k = interp1(ILr.Time, ILr.Data, t_grid);
    iLm_k = interp1(ILm.Time, ILm.Data, t_grid);
    vCr_k = interp1(VCr.Time, VCr.Data, t_grid);
    vo_k  = interp1(Vo.Time,  Vo.Data,  t_grid);
    vab_k = interp1(Vab.Time, Vab.Data, t_grid(1:end-1));
    states_all(k, :, :) = [iLr_k, iLm_k, vCr_k, vo_k];
    vab_real_all(k, :) = vab_k;

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

    n_pts = numel(vo_k);
    drift = abs(mean(vo_k(end-floor(n_pts/4)+1:end)) - mean(vo_k(1:floor(n_pts/4)))) / max(abs(mean(vo_k)), 1e-9);
    if drift > 0.02
        fprintf('  [警告] 第 %d 条 Ro=%.2fohm Co=%.0fuF 记录窗口内 vo 漂移 %.1f%%\n', k, Ro, Co*1e6, drift*100);
    end

    n_half_from_end = floor((stop_time - t_grid(1:end-1)) / half_period + 1e-9);
    v_sq_all(k, :) = Vin * (-1).^n_half_from_end;

    fprintf('[%d/%d] Ro=%.2fohm Co=%.0fuF (%.1f%%) 调节时间%.3fms 累计耗时 %.1fs\n', ...
            k, num_sequences, Ro, Co*1e6, co_ratio(k)*100, settling_time_all(k)*1e3, toc(t0));
end

states = states_all; v_sq_ideal = v_sq_all; vab_real = vab_real_all;
rload = rload_samples;
fs = FS_nominal * ones(num_sequences, 1);
vin = Vin_nominal * ones(num_sequences, 1);
duty = D * ones(num_sequences, 1);
Lr = Lr_nominal * ones(num_sequences, 1);
Cr = Cr_nominal * ones(num_sequences, 1);
Lm = Lm_nominal * ones(num_sequences, 1);
Co = co_ratio * Co_nominal;
n = n_nominal * ones(num_sequences, 1);
dt_half = dt_half_all; vo_transient = vo_transient_all; settling_time = settling_time_all;

if smoke_test
    out_name = 'llc_aging_smoke.mat';
else
    out_name = 'llc_aging_trajectory.mat';
end
out_path = fullfile(fileparts(mfilename('fullpath')), out_name);
save(out_path, 'states', 'v_sq_ideal', 'vab_real', 'rload', 'fs', 'vin', 'duty', ...
     'Lr', 'Cr', 'Lm', 'Co', 'n', 'dt_half', 'vo_transient', 'settling_time', 'co_ratio', '-v7');
fprintf('已保存至 %s\n', out_path);
