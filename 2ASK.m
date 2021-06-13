
clear;                  % 清除工作区中的数据,也就是将之前的数据清除
clc;                    % 清除命令行窗口中的数据
% 单行注释快捷键CTRL + R 取消注释CTRL + T


% init parameter
N = 1000;               % 输入信号长度，发送速率1000Baud
x = randi([0,1],1 ,N)   % 随机生成N个0或者1
s_nrz = x;              % 单极性不归零码型
fb = 1000; %发送端符号速率
fs = 16000; %滤波器采样频率fs = 16000; %滤波器采样频率
alpha = 0.25;           % 滚降系数
delay = 5;              % 时延
snr = -5;               % 信噪比-5， 5
oversamp = fs/fb;       % 过采样率
elv = 0                 % 误码率
c_error = 0;            % 误码计数
f = ((0:N-1)*fs)/N;       % 频域采样


% 使用平方根升余弦定理 可以rcosdesign(beta,span,sps)
h_sqrt = rcosine(1, oversamp, 'fir/sqrt', alpha, delay); 
% 发送端码元进行扩采样
x_oversamp = kron(s_nrz, [1, zeros(1, oversamp-1)]); 

% 计算时域使用卷积，得到脉冲成型的信号
x_shaped = conv(x_oversamp, h_sqrt, 'same');

% 画出时间域波形
figure('name','发送端脉冲成型');
subplot(2, 1, 1);
plot(x_shaped);axis([0 200 -1 1]);
title("——时域波形");
% 经过傅里叶变换后的频域波形
f_x = fft(x_shaped,N);    % f_x = fft(x_shaped, N);N 可设置为《1000
f_x_abs = abs(f_x);     % 经过傅里叶变换后的要取绝对值
% 画出频域波形
% 可以设置横坐标 
% n = 1:N; x_f = n*fs/N; step(x_f,f_x_abs);
subplot(2, 1, 2);
plot(f,f_x_abs);% axis([0 200 -6 6]);
title("——频域波形");
% ==================以上完成平方根升余弦画图================================
% ==============接下来发送端进行2ASK调制后的时域和频域波形===============

fc = 4000;                  % 载波频率

% f_up = 2000                 % 上分支载波频率
% f_down = 4000               % 下分支载波频率

x_len = length(x_shaped);     % 发送信号长度

ln = 0:x_len - 1;

t = ln/fs;                    % 时间t 2ASK调制

cari_x = cos(2*pi*fc*t);    % 载波
m_x_c = x_shaped .* cari_x;  % 模拟相乘法进行调制


% 画出载波后的波形
figure('name','2ASK调制后的波形');
subplot(2, 1, 1);
plot(m_x_c);
title("——时域波形");
axis([0 200 -0.5 0.5]);

% 傅里叶变换操作
f_m_x_c = fft(m_x_c,N);       % m代表模拟，c代表乘法
f_m_x_c_abs = abs(f_m_x_c);
subplot(2, 1, 2);
stem(f, f_m_x_c_abs);
title("——频域波形");
%axis([0 800 -0.5 0.5]);


% ==============接下来接收端进行2ASK调制后的时域和频域波形===============


m_x_c_n = awgn(m_x_c, snr, 'measured', 'db');     % 添加高斯白噪声

% 相干解调
x_c_n = m_x_c_n .* cari_x;                          % 与同频率相乘 .*

figure('name','2ASK调制后接收端的波形');
subplot(2 , 1, 1);
stem(x_c_n);
title("——时域波形");
axis([0 200 -1 1]);
f_x_c_n = fft(x_c_n, N);
f_x_c_n_abs = abs(f_x_c_n);

subplot(2, 1, 2);
stem(f, f_m_x_c_abs);
title("——频域波形");
%axis([0 200 -1 1]);



% ==============接下来接收端匹配滤波时域和频域波形===============

res = conv(x_c_n, h_sqrt);            % 接收端脉冲成型，滤波器还是使用平方根升余弦
figure('name','接收端匹配滤波');
subplot(2, 1, 1);
stem(res);
title("——时域波形");
axis([0 500 -10 20]);
% 傅里叶变换
f_res = fft(res, N);
f_res_abs = abs(f_res);
subplot(2, 1, 2);
plot(f, f_res_abs);                        % stem(x_f,f_res_abs);
title("——频域波形");


% ===============发送端输入信号和接收判决器输出的信号波形对比图============

% 抽样同步
SynPosi = delay * oversamp * 2 + 1;                 % 两个时延*过采样率
SymPosi = SynPosi + (0:oversamp:(N-1) * oversamp);  % 采样点

res_signl = res(SymPosi);                           % 接收端采样信号



% 判决
res_match = zeros(length(res_signl));               % 初始化一个矩阵后续存放0|1

for i = 1:N
    if res_signl(i) > 0.5
        res_match(i) = 1;
        
    elseif res_signl(i) <= 0.5
        res_match(i) = 0;
        
    end
end

% 画图
figure('name', '比较输入和输出信号');
subplot(2, 1, 1);
stem(s_nrz); axis([0 200 -0.5 1.5]);          % 单极性输入信号
title('单极性输入信号');

subplot(2, 1, 2);
stem(res_match); axis([0 200 -0.5 1.5]);       % 单极性输入信号
title('接收端匹配输出信号');        


% 计算误码率
for i = 1:N
    if res_match(i) ~= s_nrz(i)
        c_error = c_error + 1;
        
    end
end

elv = c_error / N;                      % 误码elv

