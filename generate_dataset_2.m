%% generate_dataset.m
% 批量运行 LLC_simu.slx，在 (Vin, fs, Lr, Cr, Lm, Co, Ro) 联合空间上采样，提取
% 每条工况记录窗口末尾若干个开关周期内的密采样波形（不再只在开关切换时刻采样），
% 保存成与 mat_to_npz.py 对应的 .mat 文件。
%
% 保存字段含义（每个字段第一维都是样本条数 num_sequences）：
%   states     - 状态波形 [iLr,iLm,vCr,vo]，(num_sequences, seq_len+1, 4)。
%                seq_len+1 个点，含记录窗口起点，单位 [A, A, V, V]
%   v_sq_ideal - 理想原边方波电压（公式合成，忽略死区，按半周期 ±Vin 交替），
%                (num_sequences, seq_len)，单位 V。第 j 列驱动 states 里第 j 步
%                到第 j+1 步的状态转移
%   vab_real   - 真实原边电压（Simulink 测量，含死区等非理想效应），
%                (num_sequences, seq_len)，单位 V，跟 v_sq_ideal 逐列对应同一组
%                时刻（每步区间取中点采样，避开边界死区跳变）
%   dt_half    - (num_sequences, 1)，states/v_sq_ideal 的采样间隔，单位 s
%   rload      - (num_sequences, 1)，负载电阻，单位 ohm，待辨识电路参数
%   fs         - (num_sequences, 1)，开关频率，单位 Hz，本版本固定不扫描
%   vin        - (num_sequences, 1)，输入母线电压，单位 V，本版本固定不扫描
%   duty       - (num_sequences, 1)，占空比，单位 %，固定不扫描
%   Lr/Cr/Lm/Co - 各 (num_sequences, 1)，谐振电感/谐振电容/励磁电感/输出滤波
%                电容，单位 H/F/H/F，待辨识电路参数
%   n          - (num_sequences, 1)，变压器变比 Np:Ns，本版本固定不扫描（每条都
%                是同一个常数），只是跟其余电路参数一样按逐样本形式存，方便以后
%                如果真要扫描 n 时不用改数据结构
%   vo_transient  - 冷启动（vo 从 0 开始）低密度 vo 曲线，(num_sequences,
%                transient_max_half+1)，每半开关周期一个点，单位 V。Co 在稳态
%                观测窗口里几乎没有痕迹（vo 纹波跟真实 Co 相关系数实测只有
%                -0.2），冷启动到稳定的过程主要由 Ro*Co 决定，是目前唯一还有
%                点用的 Co 替代信号，见 physics.py::compute_settling_time
%   settling_time - (num_sequences, 1)，冷启动 vo 进入±2%误差带且不再跑出去的
%                时刻，单位 s，从 vo_transient 直接算出来存一份，训练时不用
%                每次都重新算
%
% fs 按每条样本自己的 Lr/Cr/Lm 重新算 f1/f2 采样。D 固定不扫描。
%
% 前置条件：
%   1) 需先在 Simulink 里保存 LLC_simu.slx（Ctrl+S）
%   2) 五个 To Workspace 块变量名需为 ILr / ILm / Vo / VCr / Vab，Save format = Timeseries
%      （Vab 为原边 H 桥输出到谐振腔的电压，即 v_sq 对应的真实测量点）
%   3) 电路参数固定值与扫描范围需与 config.py 保持一致
%
% 运行方式：MATLAB 里直接 Run，或命令行 run('generate_dataset.m')

clear; clc;
mdl = 'LLC_simu';

%% ---------------------------------------------------------------------
% 冷启动瞬态段配置（跟稳态段共用同一次 sim(mdl) 结果，不用重新仿真）
% ------------------------------------------------------------------------
stop_time = 0.02;              % 原来固定 0.01s，最坏情况 Ro*Co≈10.6ms 太接近上限，加大留余量
transient_max_half = 3000;     % 冷启动瞬态最多记录的半周期数，需覆盖最坏情况的调节时间
settle_band_frac = 0.02;       % 判定"进入稳态且不再跑出去"的相对误差带
set_param(mdl, 'StopTime', num2str(stop_time));

%% ---------------------------------------------------------------------
% 标称电路参数（与 setparams.m 一致，不直接调用它以免其 clear;clc; 清掉下面的变量）
% ------------------------------------------------------------------------
Vin_nominal = 311;      % V
Lr_nominal  = 73e-6;    % H
Cr_nominal  = 35e-9;    % F
Lm_nominal  = 146e-6;   % H
Co_nominal  = 2200e-6;  % F
Ro_nominal  = 2.0;      % ohm
n_nominal   = 6.5;      % 变压器变比 Np:Ns
D           = 45;       % 占空比 (%)
Vpulse      = 1;        % 驱动脉冲幅值

%% ---------------------------------------------------------------------
% fs/Vin/Vpulse/D/n 先固定不扫描，只扫描 Lr/Cr/Lm/Co/Ro（±20%，无实际
% 规格依据，按需调整）与数据集规模
% ------------------------------------------------------------------------
FS_nominal = 90e3;      % Hz，固定开关频率，与 config.py 的 FS_NOMINAL_HZ 保持一致
lr_range  = [0.8*Lr_nominal,  1.2*Lr_nominal];
cr_range  = [0.8*Cr_nominal,  1.2*Cr_nominal];
lm_range  = [0.8*Lm_nominal,  1.2*Lm_nominal];
co_range  = [0.8*Co_nominal,  1.2*Co_nominal];
rload_range = [0.5*Ro_nominal, 2.0*Ro_nominal];

num_sequences = 300;    % 目标可用样本数，留出 mat_to_npz.py 过滤异常样本的缓冲
rng_seed      = 42;

% 每条序列覆盖记录窗口末尾 record_periods 个完整开关周期，每周期内均匀取
% pts_per_period 个点，seq_len = record_periods * pts_per_period 是最终点数
pts_per_period = 200;
record_periods = 10;
seq_len = record_periods * pts_per_period;

rng(rng_seed);
% Lr/Cr/Lm/Co/Ro 五维联合空间用拉丁超立方采样（需要 Statistics and Machine
% Learning Toolbox 的 lhsdesign），比五个独立 rand() 撒点更均匀地覆盖参数空间，
% 同样的样本数下不容易出现某些区域扎堆、某些区域完全没覆盖到的情况。
lhs_unit      = lhsdesign(num_sequences, 5);   % (num_sequences, 5)，每维在(0,1)内分层采样
lr_samples    = lr_range(1)    + (lr_range(2)    - lr_range(1))    * lhs_unit(:, 1);
cr_samples    = cr_range(1)    + (cr_range(2)    - cr_range(1))    * lhs_unit(:, 2);
lm_samples    = lm_range(1)    + (lm_range(2)    - lm_range(1))    * lhs_unit(:, 3);
co_samples    = co_range(1)    + (co_range(2)    - co_range(1))    * lhs_unit(:, 4);
rload_samples = rload_range(1) + (rload_range(2) - rload_range(1)) * lhs_unit(:, 5);
vin_samples   = Vin_nominal * ones(num_sequences, 1);
fs_samples    = FS_nominal  * ones(num_sequences, 1);
n_samples     = n_nominal   * ones(num_sequences, 1);   % 固定不扫描，按逐样本形式存，方便以后扩展

% fs 固定后，个别 Lr/Cr/Lm 组合可能让固定 fs 落到该样本自己的正常工作区
% （f1~1.05*f2）之外（比如深容性区），这里只提示不阻断，工况本身仍然合法
f1_samples = 1 ./ (2*pi*sqrt((lr_samples + lm_samples) .* cr_samples));   % 三元谐振频率
f2_samples = 1 ./ (2*pi*sqrt(lr_samples .* cr_samples));                  % 二元谐振频率

%% ---------------------------------------------------------------------
% 逐工况仿真 + 提取记录窗口末尾密采样波形
% ------------------------------------------------------------------------
states_all   = zeros(num_sequences, seq_len+1, 4);   % 列顺序 [ILr, ILm, VCr, Vo]
v_sq_all     = zeros(num_sequences, seq_len);        % 理想方波（公式合成，正向模型输入）
vab_real_all = zeros(num_sequences, seq_len);        % 真实原边电压（Simulink 测量，供对比/监督用）
dt_half_all  = zeros(num_sequences, 1);
vo_transient_all   = zeros(num_sequences, transient_max_half+1);   % 冷启动低密度vo曲线
settling_time_all  = zeros(num_sequences, 1);                       % 冷启动调节时间

t0 = tic;
for k = 1:num_sequences
    fs_k = fs_samples(k);
    Vin  = vin_samples(k);     % 供模型 Amplitude 参数引用
    Lr   = lr_samples(k);      % 供模型 Inductance 参数引用
    Cr   = cr_samples(k);      % 供模型 Capacitance 参数引用
    Lm   = lm_samples(k);      % 供模型 Inductance 参数引用
    Co   = co_samples(k);      % 供模型 Capacitance 参数引用
    Ro   = rload_samples(k);   % 供模型 Resistance 参数引用
    T    = 1 / fs_k;           % 供模型 Period/PhaseDelay 参数引用

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
        error('第 %d 条工况 fs=%.1fkHz 的 record_periods=%d 个周期时长超过模型 StopTime(%.3gs)，需加大 StopTime 或减小 record_periods', ...
              k, fs_k/1e3, record_periods, stop_time);
    end

    % fs 固定，该样本自己的 Lr/Cr/Lm 组合可能让 fs 落到正常工作区（f1~1.05*f2）之外
    if fs_k < f1_samples(k) || fs_k > 1.05*f2_samples(k)
        fprintf('  [警告] 第 %d 条 Lr=%.1fuH Cr=%.1fnF Lm=%.1fuH 下固定 fs=%.1fkHz 落在正常工作区 [%.1f, %.1f]kHz 之外\n', ...
              k, Lr*1e6, Cr*1e9, Lm*1e6, fs_k/1e3, f1_samples(k)/1e3, 1.05*f2_samples(k)/1e3);
    end

    % t_grid 取 seq_len+1 个点（含窗口起点），状态在全部点上采样，方波只在前
    % seq_len 个点上取值（第 j 个方波值驱动 t_grid(j)->t_grid(j+1) 这一步转移）
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

    % 冷启动瞬态段：跟稳态段共用这次 sim(mdl) 的结果，不用重新仿真。假设模型
    % 状态初始条件是 0（Simulink 积分器默认初始值），从 t=0 开始按半周期采样
    t_transient_grid = (0:transient_max_half)' * half_period;
    t_transient_grid = t_transient_grid(t_transient_grid <= Vo.Time(end));
    vo_transient_k = interp1(Vo.Time, Vo.Data, t_transient_grid);
    if numel(vo_transient_k) < transient_max_half + 1
        % 仿真时长不够覆盖 transient_max_half 个半周期时，末尾补最后一个值，
        % 保证每条样本存下来的长度一致
        vo_transient_k(end+1:transient_max_half+1) = vo_transient_k(end);
    end

    tail_n = max(1, floor(numel(vo_transient_k) / 15));
    vo_ss_ref = mean(vo_transient_k(end-tail_n+1:end));
    band = settle_band_frac * abs(vo_ss_ref);
    within = abs(vo_transient_k - vo_ss_ref) <= band;
    settle_idx = numel(vo_transient_k);   % 找不到就退化成"到仿真结束都没稳定"
    for j = 1:numel(within)
        if all(within(j:end))
            settle_idx = j;
            break;
        end
    end
    settling_time_k = (settle_idx - 1) * half_period;

    vo_transient_all(k, :) = vo_transient_k';
    settling_time_all(k) = settling_time_k;

    % 稳态收敛性提示：记录窗口前1/4 与后1/4 的 vo 均值漂移
    n_pts = numel(vo_k);
    drift = abs(mean(vo_k(end-floor(n_pts/4)+1:end)) - mean(vo_k(1:floor(n_pts/4)))) / max(abs(mean(vo_k)), 1e-9);
    if drift > 0.02
        fprintf('  [警告] 第 %d 条 fs=%.1fkHz Vin=%.1fV Ro=%.2fohm 记录窗口内 vo 漂移 %.1f%%，可能没收敛到稳态，该样本慎用\n', k, fs_k/1e3, Vin, Ro, drift*100);
    end

    % 理想方波：和 v_sq_ideal 同一套相位约定（以 t=stop_time 为基准，按半周期
    % 奇偶交替），忽略死区，纯公式合成，作为正向模型的已知输入
    n_half_from_end = floor((stop_time - t_grid(1:end-1)) / half_period + 1e-9);
    v_sq_all(k, :) = Vin * (-1).^n_half_from_end;
    dt_half_all(k) = win_span / seq_len;

    fprintf('[%d/%d] fs=%.1fkHz Vin=%.1fV Ro=%.2fohm Lr=%.1fuH Cr=%.1fnF Lm=%.1fuH Co=%.0fuF 序列长度%d点 调节时间%.3fms 完成，累计耗时 %.1fs\n', k, num_sequences, fs_k/1e3, Vin, Ro, Lr*1e6, Cr*1e9, Lm*1e6, Co*1e6, seq_len+1, settling_time_k*1e3, toc(t0));
end

%% ---------------------------------------------------------------------
% 保存为 .mat（供 mat_to_npz.py 转换）
% ------------------------------------------------------------------------
states     = states_all;
v_sq_ideal = v_sq_all;
vab_real   = vab_real_all;
rload      = rload_samples;
fs         = fs_samples;
vin        = vin_samples;
duty       = D * ones(num_sequences, 1);
Lr         = lr_samples;
Cr         = cr_samples;
Lm         = lm_samples;
Co         = co_samples;
n          = n_samples;
dt_half    = dt_half_all;
vo_transient  = vo_transient_all;
settling_time = settling_time_all;

out_path = fullfile(fileparts(mfilename('fullpath')), 'llc_dataset_simulink_2.mat');
save(out_path, 'states', 'v_sq_ideal', 'vab_real', 'rload', 'fs', 'vin', 'duty', ...
     'Lr', 'Cr', 'Lm', 'Co', 'n', 'dt_half', 'vo_transient', 'settling_time', '-v7');
fprintf('数据集已保存至 %s\n', out_path);
