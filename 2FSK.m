% clear;                  % 清除工作区中的数据,也就是将之前的数据清除
% clc;                    % 清除命令行窗口中的数据
% % 单行注释快捷键CTRL + R 取消注释CTRL + T
% 
% 
% % init parameter
% N = 1000;               % 输入信号长度，发送速率1000Baud
% x = randi([0,1],1 ,N)   % 随机生成N个0或者1
% s_nrz = x;              % 单极性不归零码型
% % 2FSK
% up_x = x;               % 上分支
% down = (1-x);              % 下分支
% fb = 1000;              % 发送速率1000Baud
% fs = 16000;             % 采样频率
% alpha = 0.25;           % 滚降系数
% delay = 5;              % 时延
% snr = -5;               % 信噪比-5， 5
% oversamp = fs/fb;       % 过采样率
% elv = 0                 % 误码率
% c_error = 0;            % 误码计数
% f = ((0:N-1)*fs)/N;       % 频域采样点
% 
% % 使用平方根升余弦定理 可以rcosdesign(beta,span,sps)
% h_sqrt = rcosine(1, oversamp, 'fir/sqrt', alpha, delay); 
% % 发送端码元进行扩采样
% x_oversamp = kron(s_nrz, [1, zeros(1, oversamp-1)]); 
% 
% % 上分支扩采样
% % up_x_oversamp = kron(up_x, [1, zeros(1, oversamp-1)]);
% % 下分支扩采样
% % down_x_oversamp = kron(dowm_x, [1, zeros(1, oversamp-1)]);
% 
% % 计算时域使用卷积，得到脉冲成型的信号
% x_shaped = conv(x_oversamp, h_sqrt, 'same');
% 
% % 计算上，下分支卷积
% %up_x_shaped = conv(up_x_oversamp, h_sqrt, 'same');
% %down_x_shaped = conv(dowm_x_oversamp, h_sqrt, 'same');
% 
% 
% 
% % 画上下分支谱
% %figure('name','2FSK 发送端上分支脉冲成型')
% % subplot(2, 1, 1);stem(up_x_shaped);title("上分支时间域");axis([0 200 -1 1]);
% % f_up_x_shaped = fft(up_x_shaped);f_up_x_shaped_abs = abs(f_up_x_shaped);
% % subplot(2, 1, 2);stem(f, f_up_x_shaped_abs);title('上分支频域');
% 
% %figure('name','2FSK 发送端下分支脉冲成型')
% % subplot(2, 1, 1);stem(down_x_shaped);title("下分支时间域");axis([0 200 -1 1]);
% % f_down_x_shaped = fft(dowm_x_shaped);f_down_x_shaped_abs = abs(f_down_x_shaped);
% % subplot(2, 1, 2);stem(f, f_up_x_shaped_abs);title('下分支频域');
% 
% % 画出时间域波形
% figure('name','发送端脉冲成型');
% subplot(2, 1, 1);
% plot(x_shaped);axis([0 200 -1 1]);
% title("——时域波形");
% % 经过傅里叶变换后的频域波形
% f_x = fft(x_shaped,N);    % f_x = fft(x_shaped, N);N 可设置为《1000
% f_x_abs = abs(f_x);     % 经过傅里叶变换后的要取绝对值
% % 画出频域波形
% % 可以设置横坐标 
% % n = 1:N; x_f = n*fs/N; step(x_f,f_x_abs);
% subplot(2, 1, 2);
% plot(f,f_x_abs);% axis([0 200 -6 6]);
% title("——频域波形");
% 
% 
% 
% 








clear;                  % 清除工作区中的数据,也就是将之前的数据清除
clc;                    % 清除命令行窗口中的数据
% 单行注释快捷键CTRL + R 取消注释CTRL + T


% init parameter
N = 1000;               % 输入信号长度，发送速率1000Baud
x = randi([0,1],1 ,N)   % 随机生成N个0或者1
%s_nrz = x;              % 单极性不归零码型





% 2FSK
up_x = x;               % 上分支
down_x = (1-x);           % 下分支
fb = 1000;              % 发送速率1000Baud
fs = 16000;             % 采样频率
alpha = 0.25;           % 滚降系数
delay = 5;              % 时延
snr = -5;               % 信噪比-5， 5
oversamp = fs/fb;       % 过采样率
elv = 0                 % 误码率
c_error = 0;            % 误码计数
f = ((0:N-1)*fs)/N;       % 频域采样点

% 使用平方根升余弦定理 可以rcosdesign(beta,span,sps)
h_sqrt = rcosine(1, oversamp, 'fir/sqrt', alpha, delay); 
% 发送端码元进行扩采样
%x_oversamp = kron(s_nrz, [1, zeros(1, oversamp-1)]); 

% 上分支扩采样
up_x_oversamp = kron(up_x, [1, zeros(1, oversamp-1)]);
% 下分支扩采样
down_x_oversamp = kron(down_x, [1, zeros(1, oversamp-1)]);

% 计算时域使用卷积，得到脉冲成型的信号
% x_shaped = conv(x_oversamp, h_sqrt, 'same');

% 计算上，下分支卷积
up_x_shaped = conv(up_x_oversamp, h_sqrt);
down_x_shaped = conv(down_x_oversamp, h_sqrt);






% 画上下分支谱
figure('name','2FSK 发送端上分支脉冲成型')
subplot(2, 1, 1);plot(up_x_shaped);title("上分支时间域");axis([0 1000 -0.5 0.5]);
f_up_x_shaped = fft(up_x_shaped,N);f_up_x_shaped_abs = abs(f_up_x_shaped);
subplot(2, 1, 2);plot(f, f_up_x_shaped_abs);title('上分支频域');%axis([0 1000 0 40]);

figure('name','2FSK 发送端下分支脉冲成型')
subplot(2, 1, 1);plot(down_x_shaped);title("下分支时间域");axis([0 1000 -0.5 0.5]);
f_down_x_shaped = fft(down_x_shaped,N);f_down_x_shaped_abs = abs(f_down_x_shaped);
subplot(2, 1, 2);plot(f, f_down_x_shaped_abs);title('下分支频域');%axis([0 1000 0 40]);



f_up = 2000                 % 上分支载波频率
f_down = 4000               % 下分支载波频率


x_up_len = length(up_x_shaped);           % 发送信号长度
x_down_len = length(down_x_shaped);       % 发送信号长度

ln_up = 0:x_up_len - 1;
ln_down = 0:x_down_len - 1;


t_up = ln_up/fs;            % 上分支时间变量 2FSK 调制
t_dowm = ln_down/fs;        % 下分支时间变量 2FSK 调制

% ==================2FSK调制===========================
fsk_x_shape = f_down_x_shaped + f_up_x_shaped;      % 两个分支合成
% 上下分支生成载波
carri_up = cos(2*pi*f_up*ln_up);                    % cos(2*pi*f*t)
carri_down = cos(2*pi*f_down*ln_down);
% 
m_x_up = carri_up .* up_x_shaped;
m_x_down = carri_down .* down_x_shaped;                 % 上下之路相干解调
m_up_down = m_x_up + m_x_down;                          % 调制后合成
% 输出调制后合成的时间和频域波形
figure('name','2FSK 两分支调制合成');
subplot(2, 1, 1);
stem(m_up_down);title("---时间");axis([0 1000 -1 1]);
f_m_up_down = fft(m_up_down,N);
f_m_up_down_abs = abs(f_m_up_down);
subplot(2, 1, 2);stem(f, f_m_up_down_abs);title("---频域");%axis([0 1000 0 1]);



% 信道传输，添加高斯白噪声
up_shaped_n = awgn(up_x_shaped, snr, 'measured', 'db');
down_shaped_n = awgn(down_x_shaped, snr, 'measured', 'db');         % 添加噪声的信号是调制前

% 接受端相干解调
de_up_shaped_n = up_shaped_n .* carri_up;                           % 需要同频同相对上下之路进行解调
de_down_shaped_n = down_shaped_n .* carri_down;

% 由于是2FSK是将两2ASK 相加，所有分开计算就可以，少了一步就是选择合适的频率
% 上支路时间域和频谱
figure('name','接受端上支路解调时域和频域');
subplot(2,1,1),plot(de_up_shaped_n);axis([0 1000 -1 1]);
title('----接收端上支路解调后的信号时域波形');                        % 画出时域波形
f_de_up_shaped_n=fft(de_up_shaped_n,N);                             % 进行傅里叶变换
f_de_up_shaped_n_abs=abs(f_de_up_shaped_n);
subplot(2,1,2),plot(f,f_de_up_shaped_n_abs);
title('----接受端上支路解调信号频谱');                                % 画出频谱

% 下支路
figure('name','接受端下支路解调时域和频域');
subplot(2,1,1),plot(de_down_shaped_n);axis([0 1000 -1 1]);
title('----接收端下支路解调后的信号时域波形');                         % 画出时域波形
f_de_down_shaped_n=fft(de_down_shaped_n,N);                          % 进行傅里叶变换
f_de_down_shaped_n_abs=abs(f_de_down_shaped_n);
subplot(2,1,2),plot(f,f_de_down_shaped_n_abs);
title('----接受端下支路解调信号频谱');                                 % 画出频谱

% 接收端滤波器
%上支路
recv_up=conv(de_up_shaped_n, h_sqrt);     %再次频域乘积，时域卷积
figure('name', '----接收端上支路匹配滤波');
subplot(2,1,1),plot(recv_up);axis([0 800 -1 1]);
title('----接收端上支路时域');  %画出时域波形
f_recv_up=fft(recv_up,N);
f_recv_up_abs=abs(f_recv_up);
subplot(2,1,2),plot(f,f_recv_up_abs);
title('----接收端上支路频谱');  %画出频谱

 
%下支路
recv_down=conv(de_down_shaped_n,h_sqrt);     %解调信号频域乘积，时域卷积
figure('name', '----接收端下支路匹配滤波');
subplot(2,1,1),plot(recv_down);axis([0 800 -1 1]);
title('----接收端上支路时域');  %画出时域波形
f_recv_down=fft(recv_down,N);
f_recv_down_abs=abs(f_recv_down);
subplot(2,1,2),plot(f,f_recv_down_abs);
title('----接收端上支路频谱');  %画出频谱
 

% 同步抽样
%接收端同步、采样
synPosi=delay * oversamp * 2 + 1;
symPosi=synPosi+(0:oversamp:(N-1) * oversamp);
resv_up_signal=recv_up(symPosi);
resv_down_signal=recv_down(symPosi);

% 初始化一个矩阵存放0和1
resv_match = zeros(N);

 
%接收端判决
% 比较判别器
for i=1:N
    if resv_up_signal(i) > resv_down_signal(i)          %  上支路判决值大于下支路判决值 选择上支路频率
        resv_match(i) = 1;
    elseif resv_up_signal(i) < resv_down_signal(i)      %  上支路判决值小于下支路判决值 选择下支路频率
        resv_match(i) = 0;
    end
end


figure('name', '抽样判决后波形');
subplot(2,1,1),stem(resv_match);axis([0 100 -0.5 1.5]);
title('接收端抽抽样判决后的波形'); 
 
subplot(2,1,2),stem(up_x);axis([0 100 -0.5 1.5]); % up_x 也就是单极性码元
title('发送端的原始波形'); 
 
% 计算误码率

for i = 1:N
    if resv_match(i) ~= up_x(i)
        c_error = c_error + 1;
    end
end
elv = c_error / N;
sprintf('抽样判决后误码率：%.5f',elv)
 
 
 

