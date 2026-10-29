clear;clc;

Vin = 311;
Lr = 73e-6;
Cr = 35e-9;
Lm = 146e-6;
Co = 2200e-6;
Ro = 2;
n  = 6.5;          % 变压器变比 Np:Ns，与 config.py 的 CIRCUIT_PARAMS["n"] 一致

f_1 = 1/(2*pi*sqrt((Lr+Lm)*Cr));
f_2 = 1/(2*pi*sqrt(Lr*Cr));

fprintf("当前二元谐振频率为 %f Hz \n", f_2);
fprintf("当前三元谐振频率为 %f Hz \n", f_1);

Vpulse = 1;
f = 90e3;
T = 1/f;
D = 45;

%% ---------------------------------------------------------------------
% 跑一遍标称参数下的仿真，按跟 generate_dataset.m 一样的密采样格式（同一套
% pts_per_period/record_periods，跟训练数据同分辨率）保存参考波形到
% llc_reference.mat，供 inference.py/inference_inverse.py 用 --reference-npz
% 读取做对照——直接用真实 Simulink 波形当真值，不用现场再跑一遍 Python 数值仿真
% ------------------------------------------------------------------------
mdl = 'LLC_simu';
pts_per_period = 200;   % 需与 config.py DATASET_CONFIG 保持一致
record_periods = 10;
seq_len = record_periods * pts_per_period;

% 冷启动瞬态段配置，需与 generate_dataset.m 保持一致
stop_time = 0.02;
transient_max_half = 3000;
settle_band_frac = 0.02;
set_param(mdl, 'StopTime', num2str(stop_time));

% Data Import/Export 里 "Single simulation output" 开着，To Workspace 输出
% 不会直接进工作区，要接收 sim() 返回值再取字段
simOut = sim(mdl);
ILr = simOut.ILr;
ILm = simOut.ILm;
Vo  = simOut.Vo;
VCr = simOut.VCr;
Vab = simOut.Vab;

win_span = record_periods * T;   % 记录窗口时长（record_periods 个完整周期）
if win_span > stop_time
    error('record_periods=%d 个周期时长超过模型 StopTime(%.3gs)，需加大 StopTime 或减小 record_periods', record_periods, stop_time);
end

% t_grid 取 seq_len+1 个点（含窗口起点），状态在全部点上采样，方波只在前
% seq_len 个点上取值，跟 generate_dataset.m 的窗口/相位约定完全一致
half_period = T / 2;
t_grid = linspace(stop_time - win_span, stop_time, seq_len + 1)';

iLr_k = interp1(ILr.Time, ILr.Data, t_grid);
iLm_k = interp1(ILm.Time, ILm.Data, t_grid);
vCr_k = interp1(VCr.Time, VCr.Data, t_grid);
vo_k  = interp1(Vo.Time,  Vo.Data,  t_grid);
vab_k = interp1(Vab.Time, Vab.Data, t_grid(1:end-1));

states   = reshape([iLr_k, iLm_k, vCr_k, vo_k], [1, seq_len+1, 4]);
vab_real = reshape(vab_k, [1, seq_len]);

% 冷启动瞬态段：跟 generate_dataset.m 同一套算法，共用上面这次 sim 的结果
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
    if all(within(j:end))
        settle_idx = j;
        break;
    end
end
vo_transient  = reshape(vo_transient_k, [1, transient_max_half+1]);
settling_time = (settle_idx - 1) * half_period;
fprintf('冷启动调节时间 %.3f ms\n', settling_time*1e3);

% 理想方波：和 generate_dataset.m 同一套相位约定（以 t=stop_time 为基准，按半周期
% 奇偶交替），忽略死区，纯公式合成
n_half_from_end = floor((stop_time - t_grid(1:end-1)) / half_period + 1e-9);
v_sq_ideal = reshape(Vin * (-1).^n_half_from_end, [1, seq_len]);

% 稳态收敛性提示：记录窗口前1/4 与后1/4 的 vo 均值漂移
n_pts = numel(vo_k);
drift = abs(mean(vo_k(end-floor(n_pts/4)+1:end)) - mean(vo_k(1:floor(n_pts/4)))) / max(abs(mean(vo_k)), 1e-9);
if drift > 0.02
    fprintf('[警告] 记录窗口内 vo 漂移 %.1f%%，可能没收敛到稳态\n', drift*100);
end

dt_half = win_span / seq_len;
rload   = Ro;
fs      = f;
vin     = Vin;
duty    = D;

out_path = fullfile(fileparts(mfilename('fullpath')), 'llc_reference_2.mat');
save(out_path, 'states', 'v_sq_ideal', 'vab_real', 'rload', 'fs', 'vin', 'duty', ...
     'Lr', 'Cr', 'Lm', 'Co', 'n', 'dt_half', 'vo_transient', 'settling_time', '-v7');
fprintf('参考波形已保存至 %s\n', out_path);
