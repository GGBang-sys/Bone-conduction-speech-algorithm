%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIR filtering for speech recovery using
%%% bone- and air-conducted speech signals
%%%
%%% Coded on Apr 24, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants
%VOL - Multiple channel (MRSP)
L1(1:3) = [   1       1       1 ] * 128;
L2(1:3) = [   1       1       1 ] * 64;

mu1(1:3) = [ 0.005   0.005   0.005]/3;
mu2(1:3) = [ 0.005   0.005   0.005]/3;

%FIR - Multiple channel (MRSP)
M_fir  = [1        1         1] * 128;
mu_fir = [0.0025    0.0025     0.0025];

%FLA - Multiple channel (MRSP)
M_FLA   = M_fir;       
P_max   = [20      20      20];             
mu_FLA1 = mu_fir;
mu_FLA2 = [3.0E-6  3.0E-6  3.0E-6];

%GFL - Multiple channel (MRSP)
P_gmax  = [6      6      6    ];       
N_d     = [3      3      3    ];
mu_GFL  = [0.001  0.001  0.001]/3; 
 

filename1 ='BvAh_F_bon_whi65_ch1.wav';  % noiseless bone-conducted speech
[x_org, Fs, bits] = wavread(filename1);
x(1, :) = x_org';

filename1 ='BvAh_F_bon_whi65_ch2.wav';  % noiseless bone-conducted speech
[x_org, Fs, bits] = wavread(filename1);
x(2, :) = x_org';

filename1 ='BvAh_F_bon_whi65_ch3.wav';  % noiseless bone-conducted speech
[x_org, Fs, bits] = wavread(filename1);
x(3, :) = x_org';

% Data reading (air-conducted speech)
filename2 ='BvAh_F_air_whi65_ch2.wav';  % noisyless air-conducted speech

[d, Fs, bits] = wavread(filename2);

[q, N] =size(x);   % q = number of ref channels
% q=2;              % ch1,ch2

% - Multiple channel (MRSP) - FIR filtering
disp('FIR ==>');
disp(' ');


y_q = zeros(q, N);
y   = zeros(1,N);
e   = zeros(1,N);
a   = zeros(q, max(M_fir), N);

for n=1:N
    
    % filter output
    for j=1:q
        
        for i=0:M_fir(j)-1
            if n-i > 0
                y_q(j,n) = y_q (j, n) + a(j, i+1, n) * x(j, n-i);
            end
        end
        
    end
    
    % error
    if q==1
        y(n) = y_q(1,n);
    else
         y(n) = sum( y_q(1:q, n) );
    end
    
    e(n) = d(n) - y(n);
    
    % updating
    for j=1:q
        
        for i=0:M_fir(j)-1
            if n-i > 0
                a(j, i+1, n+1) = a(j, i+1, n) + mu_fir(j) * e(n) * x(j, n-i);
            end
        end
        
    end
    
end   % n loop end

% Multiple VOL filtering
disp('VOL ==>');
disp(' ');

y_vq1 = zeros(q,N);
y_vq2 = zeros(q,N);
y_v  = zeros(1,N);
e_v  = zeros(1,N);
h1   = zeros( q, max(L1) );     % history saving
h2   = zeros( q, max(L1), max(L2) );    % no history saving

for n=1:N
    
    % VOL filters output
    for k=1:q

        y_vq1(k, n)=0;
        for i=0:L1(k)-1
            if n-i > 0
                y_vq1(k, n) = y_vq1(k, n) + h1(k, i+1) * x(k, n-i);
            end
        end

        for i=0:L2(k)-1
            for j=i:L2(k)-1
                if n-i > 0 & n-j > 0
                    y_vq2(k, n) = y_vq2(k, n) + h2(k, i+1, j+1) * x(k, n-i)*x(k, n-j);
                end
            end
        end
        
    end

    y_v(n) = sum(y_vq1(1:q, n)) + sum(y_vq2(1:q, n));
    
    % error
    e_v(n) = d(n) - y_v(n);

    % updating (1st-order kernel)
    for k=1:q

        for i=0:L1(k)-1
            if n-i > 0
                h1(k, i+1) = h1(k, i+1) + mu1(k) * e_v(n) * x(k, n-i);
            end
        end

        for i=0:L2(k)-1
            for j=i:L2(k)-1
                if n-i > 0 & n-j > 0
                    h2(k, i+1, j+1) = h2(k, i+1, j+1) ...
                        + mu2(k) * e_v(n) * x(k, n-i)*x(k, n-j);
                end
            end
        end
    end
end   % n loop end


% FLANN filtering
disp('FLANN ==>');
disp(' ');

y_fq    = zeros(q,N);
y_fq1   = zeros(q,N);
y_fq2   = zeros(q,N);
y_fla   = zeros(1,N);
e_fla   = zeros(1,N);
c_fla   = zeros(q,max(M_FLA));
c_fla1  = zeros(q,max(M_FLA),max(P_max));
c_fla2  = zeros(q,max(M_FLA),max(P_max));
save_fla   = zeros(2,N);

for n =1:N
    
    for k=1:q;

        % FIR
        for i=0:M_FLA(k)-1
            if n-i > 0
                y_fq1(k,n) = y_fq1(k,n) + c_fla(k,i+1) * x(k,n-i);
            end
        end

        % FLANN
        for j=1:P_max
            for i=0:M_FLA(k)-1
                if n-i > 0
                    y_fq2(k,n) = y_fq2(k,n) + c_fla1(k,i+1,j) * cos(j*pi*x(k, n-i))...
                        +c_fla2(k,i+1,j) * sin(j*pi*x(k, n-i));
                end
            end
        end
        
    end
    
    y_fla(n) = sum(y_fq1(1:q,n)) + sum(y_fq2(1:q,n));
    
    % error
    e_fla(n) = d(n) - y_fla(n);
    
    % updating
    for k=1:q

        for i=0:M_FLA(k)-1
            if n-i > 0
                c_fla(k,i+1) = c_fla(k,i+1) + mu_FLA1(k) * e_fla(n) * x(k, n-i);
            end
        end

        for j=1:P_max
            for i=0:M_FLA(k)-1
                if n-i > 0
                    c_fla1(k,i+1,j) = c_fla1(k,i+1,j) + mu_FLA2(k) * e_fla(n) * cos(j*pi*x(k, n-i));
                    c_fla2(k,i+1,j) = c_fla2(k,i+1,j) + mu_FLA2(k) * e_fla(n) * sin(j*pi*x(k, n-i));
                end
            end
        end
    end
    
    save_fla(1,n) = c_fla1(1,1);   %save‚·‚é‚±‚Æ‚ÅŽû‘©‹ï‡‚ðŒ©‚é
    save_fla(2,n) = c_fla2(1,1);
    
end  % n loop end


% GFLANN filtering
disp('GFLANN ==>');
disp(' ');
 
c_gfl      = zeros(q,max(M_FLA));   % FIR
c_gfl1     = zeros(q,max(M_FLA),max(P_max));  % FLANN
c_gfl2     = zeros(q,max(M_FLA),max(P_max));
c_cro_cos  = zeros(max(N_d), max(M_FLA), max(P_gmax));   % GFLANN cross-terms
c_cro_sin  = zeros(max(N_d), max(M_FLA), max(P_gmax));

y_gfl_FIRq = zeros(q,N);
y_gfl_FLAq = zeros(q,N);
y_gfl_GFLq = zeros(q,N);
y_gfl      = zeros(1,N);
e_gfl      = zeros(1,N);

for n = 1:N
    
    for k=1:q

        % FIR
        for i=0:M_FLA(k)-1
            if n-i > 0
                y_gfl_FIRq(k,n) = y_gfl_FIRq(k,n) +c_gfl(k,i+1) * x(k,n-i);
            end
        end

        % FLANN
        for j=1:P_max
            for i=0:M_FLA(k)-1
                if n-i > 0
                    y_gfl_FLAq(k,n) = y_gfl_FLAq(k,n) + c_gfl1(k,i+1, j) * cos(j*pi*x(k, n-i))...
                                                      + c_gfl2(k,i+1, j) * sin(j*pi*x(k, n-i));
                end
            end
        end

        %GFLANN - cross-terms
        for pg=1:P_gmax
            for j=1:N_d
                for i=1 : M_FLA(k)-j
                    if n-j-i > 0
                        y_gfl_GFLq(k,n) = y_gfl_GFLq(k,n) ...
                          + c_cro_cos(j,i, pg) * x(k,n-j-i) * cos(pg*pi*x(k, n-i)) ...
                          + c_cro_sin(j,i, pg) * x(k,n-j-i) * sin(pg*pi*x(k, n-i));
                    end
                end
            end
        end
    end
    
    y_gfl(n) = sum(y_gfl_FIRq(1:q,n)) + sum(y_gfl_FLAq(1:q,n)) + sum(y_gfl_GFLq(1:q,n));
    
    % error
    e_gfl(n) = d(n) - y_gfl(n);
    
    % updating : FIR
    for k=1:q

        for i=0:M_FLA(k)-1
            if n-i > 0
                c_gfl(k,i+1) = c_gfl(k,i+1) + mu_GFL(k) * e_gfl(n) * x(k,n-i);
            end
        end
 
        %updating : FLANN
        for j=1:P_max
            for i=0:M_FLA(k)-1
                if n-i > 0
                    c_gfl1(k,i+1,j) = c_gfl1(k,i+1,j) + mu_FLA2(k) * e_gfl(n) * cos(j*pi*x(k, n-i));
                    c_gfl2(k,i+1,j) = c_gfl2(k,i+1,j) + mu_FLA2(k) * e_gfl(n) * sin(j*pi*x(k, n-i));
                end
            end
        end


        % updating : GFLANN - cross-terms
        for pg=1:P_gmax
            for j=1:N_d
                for i=1 : M_FLA(k)-j
                    if n-j-i > 0

                        c_cro_cos(j, i, pg) =  c_cro_cos(j, i, pg) ...
                           + mu_GFL(k) * e_gfl(n) * x(k, n-j-i) * cos(pg*pi*x(k, n-i));

                        c_cro_sin(j, i, pg) =  c_cro_sin(j, i, pg) ...
                           + mu_GFL(k) * e_gfl(n) * x(k, n-j-i) * sin(pg*pi*x(k, n-i));

                    end
                end
            end
        end
    end
    
end   % n loop end


% Plotting
figure(1)
subplot(6,1,1);
plot(1:N, d(1:N), '-');
axis( [1  N  -0.8  0.8] );
ylabel('Air-');

subplot(6,1,2);
plot(1:N, x(2,1:N), '-');
axis( [1  N  -0.8  0.8] );
ylabel('Bone-');

subplot(6,1,3);
plot(1:N, y, '-');
axis( [1  N  -0.8  0.8] );
ylabel('FIR_m');

subplot(6,1,4);
plot(1:N, y_v, '-');
axis( [1  N  -0.8  0.8] );
ylabel('VOL_m');

subplot(6,1,5)
plot(1:N, y_fla, '-');
axis( [1  N  -0.8  0.8] );
ylabel('FLA_m');

subplot(6,1,6)
plot(1:N, y_gfl, '-');
axis( [1  N  -0.8  0.8] );
ylabel('GFL_m');


wavwrite(y, Fs, bits,'FIR_BvAh2_F_w65_q123p1.wav');
wavwrite(y_v, Fs, bits,'VOL_BvAh2_F_w65_q123p1.wav');
wavwrite(y_fla, Fs, bits,'FLA_BvAh2_F_w65_q123p1.wav');
wavwrite(y_gfl, Fs, bits,'GFL_BvAh2_F_w65_q123p1.wav');


%time fre analysis

figure(2)
specgram(d(1:N), 512, Fs, 500);

figure(3)
specgram(x(2,1:N), 512, Fs, 500);

figure(4)
specgram(y(1:N), 512, Fs, 500);

figure(5)
specgram(y_v(1:N), 512, Fs, 500);

figure(6)
specgram(y_fla(1:N), 512, Fs, 500);

figure(7)
specgram(y_gfl(1:N), 512, Fs, 500);

% figure(8)
% Para_no = 10;
% for i=1:q
%     subplot(1,1+q,1+i);
%     wk_a(1, 1:N) = a(i, Para_no, 1:N);
%     plot(1:N, wk_a, '-');
% end
   

% Correlation coefficients
N_eva = floor(N/2);
X1 = y(N_eva:N)-mean(y(N_eva:N));
X2 = d(N_eva:N)'-mean(d(N_eva:N));
Coe_FIR_d = mean( X1 .* X2)/( sqrt(var(X1))*sqrt(var(X2)) )

X1 = y_v(N_eva:N)-mean(y_v(N_eva:N));
X2 = d(N_eva:N)'-mean(d(N_eva:N));
Coe_VOL_d = mean( X1 .* X2)/( sqrt(var(X1))*sqrt(var(X2)) )

X1 = y_fla(N_eva:N)-mean(y_fla(N_eva:N));
X2 = d(N_eva:N)'-mean(d(N_eva:N));
Coe_FLA_d = mean( X1 .* X2)/( sqrt(var(X1))*sqrt(var(X2)) )

X1 = y_gfl(N_eva:N)-mean(y_gfl(N_eva:N));
X2 = d(N_eva:N)'-mean(d(N_eva:N));
Coe_GFL_d = mean( X1 .* X2)/( sqrt(var(X1))*sqrt(var(X2)) )
