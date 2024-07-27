%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% l-form Adaptive Notch Filter Bank
%%%(CPANFB)
%%% Normalized Gradient(NG) algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Constans definitions
omega_0=[0.2 0.3 ]*pi;     %fre of sinusoid
Amp    =[1.0    0.5   0.15];       %Amplitude
Phs    =[pi/6 pi/3 pi/2];      %Phase
std_v  =0.1;             %additive noise std deviation

[wk, p] = size(omega_0);     %p=number of real frequencies

N           = 100000;           %N=adaptation lengh
T           = 1;               %number of independent runs
p_set    = p+1;               %number of active cells

save_mean_NG = zeros(p_set,N);      %mean of NG
save_MSE_NG = zeros(p_set,N);      %MSE(mean squerad error)of NG
save_Pow         = zeros(p_set, N);     % powers of extracted sinusoids

r      = 0.95;           %user para rho  
a_true = -2*cos(omega_0); %ideal filter para 

mu_NG  = 0.0025;           %step size for NG
alphh     = 0.98;            %user para for NG (forgetting factor)
epsil      = 0.01;            %user para for 

alphh_A2 = 0.9995;        % for power estimation

for t=1:T
    
    disp('t==>');disp(t);
    
    %Adaptive IIR notch filterring
    x_clean=zeros(1,N);      %noise-free sinusoid
    x      =zeros(1,N);      %noisy sinusoid
    v      =zeros(1,N);      %additive white noise
    
    e_pre           = zeros(p_set, N);         % pre-filters errors
    x_band_pre = zeros(p_set, N);         % pre-filters bandpass signals
    
    x_NG           = zeros(p_set, N);  % active filters input 
    
    a_NG         = zeros(p_set,N);      %filter coefficient
%     for k=1:p_set
%         a_NG(k, 1) = -2 * cos( 0.15*pi*k );
%         a_NG(k, 2) = -2 * cos( 0.15*pi*k );
%         a_NG(k, 3) = -2 * cos( 0.15*pi*k );
%     end
     
    e                 = zeros(p_set, N);      % active cells errors
    g_NG          = zeros(p_set,N);       %gradient signal
    G_NG         = zeros(p_set,N);       %low-pass filtered g_NG^2       
    x_band       = zeros(p_set, N);     % bandpass signals
    Pow            = zeros(p_set, N);      % sinusoids powers
    
    for n=3:N
        %x_clean(n)
        for i = 1:p
              x_clean(n) = x_clean(n)+...
                                   Amp(i)*sin(omega_0(i)*(n-1)+Phs(i));
        end
        
        %v(n)
        v(n) = randn(1,1)*std_v;        %additive white noise
        
        %x(n)
        x(n) = x_clean(n)+v(n);    
        
        %%%%%%%%%%%
        %NG
        %%%%%%%%%%%
        
        % Pre-filtering
        for k=1:p_set
            if k==1
                e_pre(k,n) = -r*a_NG(k,n)*e_pre(k,n-1)-(r^2)*e_pre(k,n-2)...
                                    +x(n)+a_NG(k,n)*x(n-1)+x(n-2);                
                
                x_band_pre(k, n) = x(n) - e_pre(k,n);
            else
                e_pre(k,n) = -r*a_NG(k,n)*e_pre(k,n-1)-(r^2)*e_pre(k,n-2)...
                                    +e_pre(k-1,n)+a_NG(k,n)*e_pre(k-1,n-1)+e_pre(k-1,n-2);             
                
                x_band_pre(k, n) = e_pre(k-1, n) - e_pre(k,n);
            end
        end
        
        % Preparations for parallel active cells
        for k=1:p_set
            x_NG(k,n) = x(n) - sum( x_band_pre(1:p_set, n) ) + x_band_pre(k, n);
        end
        
        % Parallel active cells
        for k=1:p_set           
            e(k,n) = -r*a_NG(k,n)*e(k,n-1)-(r^2)*e(k,n-2)...          % notch filter
                        + x_NG(k,n)+a_NG(k,n)*x_NG(k,n-1)+x_NG(k,n-2);
            
            g_NG(k,n) = -r*e(k,n-1)+x_NG(k,n-1);   % gradcient signal
            
            x_band(k, n) = x_NG(k, n) - e(k,n);       % bandpass signal
            
            G_NG(k,n) = alphh*G_NG(k,n-1)+(1-alphh)*g_NG(k,n)^2;    %AR(1)model
            
            a_NG(k,n+1) = a_NG(k,n)-mu_NG*e(k,n)*g_NG(k,n)/(epsil+G_NG(k,n));   %NG algorithm
            
            Pow(k, n) = alphh_A2 * Pow(k, n-1) +(1-alphh_A2) * x_band(k, n)^2;   % power estimation
        end % k loop end
        
    end % n loop end 
   
    for k=1:p_set
        if k<p+1
            save_mean_NG(k,1:N) = save_mean_NG(k,1:N) + a_NG(k,1:N)/T;
            
            save_MSE_NG(k,1:N) = save_MSE_NG(k,1:N)+...
                (a_NG(k,1:N)-a_true(k)) .* (a_NG(k,1:N)-a_true(k))/T;
        else
            save_mean_NG(k,1:N) = save_mean_NG(k,1:N) + a_NG(k,1:N)/T;
            
            save_MSE_NG(k,1:N) = save_MSE_NG(k,1:N)+...
                (a_NG(k,1:N)-0) .* (a_NG(k,1:N)-0)/T;
            
        end
    end
    
end  %t loop end

%Ploting
% figure(1);
% subplot(4,1,1);
% plot(1:N,x_clean,'-');
% 
% subplot(4,1,2);
% plot(1:N,v,'-');
% 
% subplot(4,1,3);
% plot(1:N,x,'-');
% 
% subplot(4,1,4);    
% plot(1:N,e_NG(p_set,1:N),'-');


% figure(2);
% for k=1:p_set
%     if k<p+1
%         plot(1:N, acos(-a_NG(k,1:N)/2)/pi,'-',1:N,ones(1,N)*omega_0(k)/pi,'--'); hold on
%     else
%         plot(1:N, acos(-a_NG(k,1:N)/2)/pi,'-'); hold on
%     end
% end
% hold off

figure (3);
for k=1:p_set
    subplot(1, p_set, k);
    plot(1:N, acos(-a_NG(k,1:N)/2)/pi,'-'); hold on
    axis( [1   N  0.05  0.6 ] );
    for k1=1:p
        plot(1:N, ones(1,N)*omega_0(k1)/pi,'--'); hold on
        axis( [1   N  0.05  0.6 ] );
    end  
    hold off
end


figure (4);
for k=1:p_set
    subplot(1, p_set, k);
    plot(1:N, Pow(k, 1:N),'-'); hold on
    for k1=1:p
        plot(1:N, ones(1,N)*(Amp(k1)^2/2),'--'); hold on
    end  
    hold off
end


%%%%%%

figure (5)
for k=1:p_set
    plot(1:N, acos(-a_NG(k,1:N)/2)/pi,'-', 'LineWidth',1); hold on
    axis( [1   N  0.05  0.6 ] );
end
for k1=1:p
    plot(1:N, ones(1,N)*omega_0(k1)/pi,'--', 'LineWidth',1); hold on
    axis( [1   N  0.05  0.6 ] );
end
hold off
set(gca, 'FontSize', 20, 'FontName', 'Arial','LineWidth',1);
%legend('VSS-LMS-A', 'VSS-LMS-F', 'VSS-LMS-G', 'SVSS-FXLMS');
xlabel('Iteration number n', ...
    'fontName', 'Arial', 'fontWeight', 'Bold', 'fontSize', 20);
ylabel('Frequency (* \pi)', ...
    'fontName', 'Arial', 'fontWeight', 'Bold', 'fontSize', 20);

