% ToDo
% NTF(多チャンネル)対応

% betaDでのchord.wavがおかしい?
% できればforなしに改良(ifft)

function [] = nmf(filename)

[y,Fs] = audioread(filename);
% 分離数
k = 3;
beta = 2;
if beta < 1
    phi = 1/(2-beta);
elseif beta >= 1 && beta <= 2
    phi = 1;
elseif beta > 2
    phi = 1/(beta-1);
end

% 元波形プロット
figure(1);
t = 0:1/Fs:(length(y)-1)/Fs;
subplot(k+1,1,1);
plot(t,y);
xlim([0 t(end)]);

% Spectrogram
fftlen = 1024;
noverlap = 512;
Y = myspectrogram(y);
X = abs(Y); % パワースペクトル

figure(2);
subplot(k+1,1,1);
imagesc(0:1/Fs:(length(y)-1)/Fs, 0:1000:Fs/2, X);
set(gca,'YDir','normal');

% 初期値
[n,m] = size(X);
H = rand(n,k);
U = rand(k,m);
max = 100; % 反復回数

% 二乗誤差
% for i = 1 : max
%     U = U.*(H'*X)./(H'*H*U);
%     H = H.*(X*U')./(H*U*(U)');
% end

% beta
for i = 1:max
    U = U.*((H'*(X.*(H*U).^(beta-2))) ./ (H'*(H*U).^(beta-1))).^phi;
    H = H.*(((X.*(H*U).^(beta-2))*U') ./ ((H*U).^(beta-1)*U')).^phi;
end

% X = H*U;
% figure();
% imagesc(0:1/Fs:(length(y)-1)/Fs, 0:1000:Fs/2, X);
% set(gca,'YDir','normal');

% Wiener Filter
W = cell(1,k);
for i = 1:k
    W{i} = Y.*(H(:,i).*U(i,:))./(H*U);
    figure(2);
    subplot(k+1,1,i+1);
    imagesc(0:1/Fs:(length(y)-1)/Fs, 0:1000:Fs/2, abs(W{i}));
    title('Wiener Filter');
    set(gca,'YDir','normal');
end

X = abs(Y)-abs(W{1});
figure();
imagesc(0:1/Fs:(length(y)-1)/Fs, 0:1000:Fs/2, X);
set(gca,'YDir','normal');

% IFFT
a = zeros(1,m*fftlen-(m-1)*noverlap);
for p = 1:k
    X = H(:,p)*U(p,:);
    X = X.*Y./abs(Y);
    [~,m] = size(X);
    z = zeros(1,m*fftlen-(m-1)*noverlap);
    for i = 1:m
        L = X(:,i);
        R = flipud(conj(X(2:end-1,i)));
        j = 1+(fftlen-noverlap)*(i-1);
        z(j:j+fftlen-1) = z(j:j+fftlen-1)+ifft([L;R])';
    end
    
    % 音声ファイル作成
%     soundsc(z,Fs);
    myfilename = sprintf('file%d.wav', p);
    audiowrite(myfilename,z,Fs);
    
    % 波形プロット
    figure(1);
    t = 0:1/Fs:(length(z)-1)/Fs;
    subplot(k+1,1,p+1);
    plot(t,z);
    xlim([0 t(end)]);
    ylim([-1 1]);
    
    % 分離音源を結合
    a = a+z;
end

% 再結合音声
% soundsc(a,Fs);
% figure();
% t = 0:1/Fs:(length(z)-1)/Fs;
% subplot(k+1,1,p+1);
% plot(t,a);
% xlim([0 t(end)]);
end


function spec = myspectrogram(y)
% 初期値
frame_size = 1024; % フレーム幅
frame_shift = 512; % フレームシフト幅
fft_point = 1024; % FFT点数
shift_times = floor((length(y)-frame_shift)/(frame_size-frame_shift)); % シフト回数
window = hann(frame_size); %ハニング窓
% result_fft = zeros(ceil((fft_point+1)/2), shift_times); % 出力データ

% FFT
i = 1:frame_size;
j = 1:shift_times;
I = repmat(i',1,shift_times);
J = repmat(j,frame_size,1);
Y = y(I+frame_shift*(J-1));
WINDOW = repmat(window,1,shift_times);
result_fft = fft(Y.*WINDOW);

spec = result_fft(1:ceil((fft_point+1)/2),:);

end
