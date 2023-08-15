
close all;
clear global;
clc;

data = load('salt.txt');
data_backup = data;
n_row = size(data,1);
n_col = size(data,2);

%% calculate the end position of each data
end_index_of_data = zeros(n_col,1);
for n = 2:n_col
    end_index_index = find(data(201:end,n)==0,1,'first');
    if isempty(end_index_index)
        end_index_index = n_row;
    else
        end_index_index = 200+end_index_index-1;
    end
    end_index_of_data(n) = end_index_index;
end

%% detect and replace the outlier points
window = 3; % by moving window
for n = 2:n_col
    for m = 50:end_index_of_data(n)-window
        if abs(data(m,n)-mean(data(m-window:m-1,n)))/data(m,n)>0.05 && abs(data(m,n)-mean(data(m+1:m+window,n)))/data(m,n)>0.05
            data(m,n) = (mean(data(m-window:m-1,n))+mean(data(m+1:m+window,n)))/2;
        end
    end
end
%% detect the max value
time_max = zeros(n_col-1,3);
for n= 2:n_col
    [max_val,max_index] = max(data(:,n));
    time_max(n-1,1) = max_index; % column 1: index
    time_max(n-1,2) = max_val; % colunme 2: max OD600
      if max_index>1
         time_max(n-1,3) = (data(max_index-1,n)+data(max_index,n)+data(max_index+1,n))/3;% mean of 3 consecutive values around the max
      else
         time_max(n-1,3)= max_val;
      end
end
   
    
%%
data_new = data;
% combine the data to make the curves more smooth.
for n = 3:n_row-2
    data_new(n,:) = (data(n-2,:)+data(n-1,:)+data(n,:)+data(n+1,:)+data(n+2,:))/5;% average consecutive 5 values
end
% caculate the gradient.
window = 3;
gradient = zeros(size(data));
for n = window+1:n_row-window
    gradient(n,:) = (sum(data_new(n+1:n+window,:))-sum(data_new(n-window:n-1,:)))/(window+1)/window;
end
% the average of the 11 continous gradient
gradient_new = gradient;
window = 5;
for n = window+1:n_row-window
    gradient_new(n,:) = sum(gradient(n-window:n+window,:))/(2*window+1);
end


%%
start_end_decrease = zeros(n_col,2);
%assign the value to the start and end of the decreasing phase.
for n = 2:n_col
    
    start_index = 0;
    end_index = 0;
    
    [max_val,max_index] = max(gradient_new(:,n));
    start_index = max_index;
    while start_index<n_row
        for m = start_index+1:n_row
            if gradient_new(m,n)<=0
                start_index = m;
                break;
            end 
        end
        if m == n_row
            break;
        end
        for end_index = start_index+1:n_row
            if gradient_new(end_index,n)>-0 % || ((m-start_index)>5&&gradient_new(m,n)>-0.002)
%                 end_index = m;
                break;
            end
        end
        
        if end_index > end_index_of_data(n)
            end_index = end_index_of_data(n);
        end
        if (end_index-start_index)>10
            start_end_decrease(n,1) = start_index;
            start_end_decrease(n,2) = end_index;
            break;
        end
    end
    
end

%%
% locate the start and end of the steady phage.
start_end_steady = zeros(n_col,2);
steady_gradient_threshold = 0.0015;
for n = 2:n_col
    
 
    start_index = 1;
    if start_end_decrease(n,2)>0
        start_index = start_end_decrease(n,2)+1;
    end
    
    end_index = start_index;
    steady_length = 0;
    i = start_index;
    while i<end_index_of_data(n)
        
        if abs(gradient_new(i,n))>steady_gradient_threshold
            i = i+1;
            continue;
        end
        
        for m = i+1:end_index_of_data(n)
            if abs(gradient_new(m,n))>steady_gradient_threshold
                if m-1-i>steady_length
                    start_index = i;
                    end_index = m-1;
                    steady_length = m-1-i;
                end
                i = m+1;
                break;
            end
        end

        if m==end_index_of_data(n)
            if end_index_of_data(n)-i>steady_length
                start_index = i;
                end_index = end_index_of_data(n);
            end
            break;
        end
    end
    
    while end_index-start_index>20
        start_index_temp = start_index;
        end_index_temp = start_index;
        period = 0;
        i=start_index;
        for m = start_index+1:end_index
            if gradient_new(m,n)*gradient_new(i,n)<=0
                if m-i-1>period
                    start_index_temp = i;
                    end_index_temp = m-1;
                    period = m-i-1;
                    i = m;
                else
                    i = m;
                end
            end
        end
        
        if period<80
            break;
        else
            if start_index_temp-start_index<end_index-end_index_temp
                start_index = end_index_temp+1;
            else
                end_index = start_index_temp-1;
            end
        end
    end
    
    if (end_index-start_index)>20 
        data_temp = gradient_new(start_index:end_index,n);
        if abs(sum(data_temp>0)-sum(data_temp<0))/length(data_temp)<0.8;
            start_end_steady(n,1) = start_index;
            start_end_steady(n,2) = end_index;
        end
    end
end

start_end_steady(start_end_steady(:,2)>500,2) = 500;


%%
start_end_increase = zeros(size(data,2),2);
start_end_increase_time = zeros(size(data,2),2);
threshold1 = 0.005;% normal: 0.02; salt: 0.005; hungry: 0.01
threshold2 = 0.005;% normal: 0.1; salt: 0.005;  hungry: 0.01 
for n = 1:size(data,2)
    for m = 2:size(data,1)
        if data(m,n)-data(m-1,n)>threshold1 && data(m+1,n)-data(m,n)>threshold1
            start_end_increase(n,1) = m;
            start_end_increase_time(n,1) = data(m,1);
            break;
        end
    end
end

for n = 1:size(data,2)
    for m = 2:size(data,1)
        if data(m,n)<0.5*max(data(:,n))
           continue; 
        end
        if data(m,n)-data(m-1,n)<threshold2  && data(m+1,n)-data(m,n)<threshold2 && data(m+2,n)-data(m+1,n)<threshold2
            start_end_increase(n,2) = m;
            start_end_increase_time(n,2) = data(m,1);
            break;
        end
    end
end

result_of_increase_stage = zeros(size(data,2)-1,13);
for n = 2:size(data,2)
    if start_end_increase(n,1)==0
        continue;
    end
    start_point = start_end_increase(n,1);
    start_time = start_end_increase_time(n,1);
    start_value = data(start_end_increase(n,1),n);
    end_point = start_end_increase(n,2);
    end_time = start_end_increase_time(n,2);
    end_value = data(start_end_increase(n,2),n);

    x = data(start_point:end_point,1);
    y = data(start_point:end_point,n);
    % exponential??
    fit_exponential  = @(beta,x) beta(1)*(exp(beta(2)*x));
    beta0 = [0 0];  
 
    [beta,resnorm,residual,exitflag,output]= lsqcurvefit(fit_exponential,beta0,x,y);

    sum(residual.^2);
    plot(x,fit_exponential(beta,x), '-s');
    
    a = log2(end_value/start_value)/(end_time-start_time);

    result_of_increase_stage(n-1,1) = data(start_point,1);
    result_of_increase_stage(n-1,2) = start_value;
    result_of_increase_stage(n-1,3) = data(end_point,1);
    result_of_increase_stage(n-1,4) = end_value;
    result_of_increase_stage(n-1,5) = data(end_point+1,1);
    result_of_increase_stage(n-1,6) = data(end_point+1,n);
    result_of_increase_stage(n-1,7) = data(end_point+2,1);
    result_of_increase_stage(n-1,8) = data(end_point+2,n);
    result_of_increase_stage(n-1,9) = a;
    result_of_increase_stage(n-1,10) = mean(data(end_point:end_point+4,n)); % mean of following 4 points
    result_of_increase_stage(n-1,11) = mean(data(end_point:end_point+9,n)); % mean of following 9 points
    result_of_increase_stage(n-1,12) = beta(1);%Y0
    result_of_increase_stage(n-1,13) = beta(2);%r calculated from exponential model
end

%%
% locate the start and end of the increasing phase?
start_end_rise = zeros(n_col,2);
for n = 2:n_col
    
    if start_end_steady(n,1)==0 && start_end_decrease(n,1)==0
        start_end_rise(n,1) = start_end_increase(n,2)+1;
        start_end_rise(n,2) = end_index_of_data(n)-1;
    else%if start_end_decrease(n,1)==0
        end_index_index = find(gradient_new(21:end,n)<0,1,'first');
        end_index_index = 20+end_index_index;
        if end_index_index-start_end_increase(n,2)>30
            start_end_rise(n,1) = start_end_increase(n,2)+1;
            start_end_rise(n,2) = end_index_index;
        end
    end
    
end

N =46;

figure
plot(data(:,N));
hold on;
plot(zeros(n_row,1),'r');
plot(data_new(:,N),'r');
% plot(gradient(:,N)*5,'k');
gradient_new(gradient_new>0.02)=0.02;
gradient_new(gradient_new<-0.02)=-0.02;
plot(gradient_new(:,N)*40,'b');

if start_end_steady(N,1)>0
    plot(start_end_steady(N,1),data(start_end_steady(N,1),N),'rd');
    plot(start_end_steady(N,2),data(start_end_steady(N,2),N),'rd');
end

if start_end_rise(N,1)>0
    plot(start_end_rise(N,1),data(start_end_rise(N,1),N),'ro');
    plot(start_end_rise(N,2),data(start_end_rise(N,2),N),'ro');
end

figure;
plot(data(:,N));
hold on;
start_end_decrease(N,:)
if start_end_decrease(N,1)>0
    plot(start_end_decrease(N,1),data(start_end_decrease(N,1),N),'r*');
    plot(start_end_decrease(N,2),data(start_end_decrease(N,2),N),'r*');
end

if start_end_steady(N,1)>0
    plot(start_end_steady(N,1),data(start_end_steady(N,1),N),'rd');
    plot(start_end_steady(N,2),data(start_end_steady(N,2),N),'rd');
end

if start_end_rise(N,1)>0
    plot(start_end_rise(N,1),data(start_end_rise(N,1),N),'ro');
    plot(start_end_rise(N,2),data(start_end_rise(N,2),N),'ro');
end

plot(start_end_increase(N,1),data(start_end_increase(N,1),N),'r-s');
plot(start_end_increase(N,2),data(start_end_increase(N,2),N),'r-s');


%calculate the decreasing rate of the decrease phase.
result_of_decrease_stage = zeros(size(data,2)-1,2);
for n = 2:n_col
    start_index = start_end_decrease(n,1);
    end_index = start_end_decrease(n,2);
    if start_index==0
        continue;
    end
    x = data(start_index:end_index,1);
    y = data(start_index:end_index,n);
    [a,b]=polyfit(x,y,1);

    fit_exponential  = @(beta,x) beta(1)*(exp(beta(2)*x));
    beta0 = [0 0];  
 
    [beta_decrease,resnorm,residual,exitflag,output]= lsqcurvefit(fit_exponential,beta0,x,y);

    sum(residual.^2);
    plot(x,fit_exponential(beta_decrease,x), '-s');

    result_of_decrease_stage(n-1,1) = a(1)*2; %% original decrease rate calculated by Junjun
    result_of_decrease_stage(n-1,2) = beta_decrease(2); %% r from exponential model

end

gradient_steady = zeros(n_col-1,1);
avg_steady = zeros(n_col-1,1);
sd_steady = zeros(n_col-1,1);
for n = 2:n_col
    start_index = start_end_steady(n,1);
    end_index = start_end_steady(n,2);
    if start_index==0
        continue;
    end
    x = data(start_index:end_index,1);
    y = data(start_index:end_index,n);
    [a,b]=polyfit(x,y,1);
    gradient_steady(n-1) = a(1)*2;
    avg_steady(n-1) = mean(y);
    sd_steady(n-1) = std(y);
end

gradient_rise = zeros(n_col-1,1);
for n = 2:n_col
    start_index = start_end_rise(n,1);
    end_index = start_end_rise(n,2);
    if start_index==0
        continue;
    end
    x = data(start_index:end_index,1);
    y = data(start_index:end_index,n);
    [a,b]=polyfit(x,y,1);
    gradient_rise(n-1) = a(1)*2;
end

time_decrease = zeros(n_col-1,2);
for n = 1:n_col-1
    if start_end_decrease(n+1,1)>0
        time_decrease(n,1) = data(start_end_decrease(n+1,1),1);
        time_decrease(n,2) = data(start_end_decrease(n+1,2),1);
    end
end
time_steady = zeros(n_col-1,2);
for n = 1:n_col-1
    if start_end_steady(n+1,1)>0
        time_steady(n,1) = data(start_end_steady(n+1,1),1);
        time_steady(n,2) = data(start_end_steady(n+1,2),1);
    end
end
time_rise = zeros(n_col-1,2);
for n = 1:n_col-1
    if start_end_rise(n+1,1)>0
        time_rise(n,1) = data(start_end_rise(n+1,1),1);
        time_rise(n,2) = data(start_end_rise(n+1,2),1);
    end
end
result = [result_of_increase_stage,time_max,time_decrease,result_of_decrease_stage,...
    time_steady,gradient_steady,avg_steady,sd_steady,...
    time_rise,gradient_rise];


%%
% to locate the start and the end - logistic growth
start_end_increase = zeros(size(data,2)-1,5);
threshold = 0.005;
for n = 2:size(data,2)
    for m = 2:size(data,1)-1
        if data(m,n)-data(m-1,n)>threshold && data(m+1,n)-data(m,n)>threshold
            start_end_increase(n-1,1) = m;
            a = m;
            break;
        end
    end
    for m = 2:size(data,1)
        if data(m,n)<0.5*max(data(:,n))
           continue; 
        end
       if data(m,n)-data(m-1,n)<threshold  && data(m+1,n)-data(m,n)<threshold && data(m+2,n)-data(m+1,n)<threshold
            start_end_increase(n-1,2) = m;
            b = m;
            break;
       end
    end
    x = data(a:b,1);
%     x = x/max(x);
    y = data(a:b,n);

    plot(x,y,'-r*')
    hold on 
    % Logistic ??
    fit_logistic  = @(beta,x) beta(1)./(1+exp(beta(2)+beta(3)*x));
    beta0 = [0 0 0];  
 
    [beta,resnorm,residual,exitflag,output]= lsqcurvefit(fit_logistic,beta0,x,y);
    
    start_end_increase(n-1,3) = beta(1);
    start_end_increase(n-1,4) = beta(2);
    start_end_increase(n-1,5) = beta(3);
    
    sum(residual.^2);
    plot(x,fit_logistic(beta,x), '-s');

%     fit_exponential  = @(alpha,x) alpha(1)+exp(alpha(2)+alpha(3)*x);
%     alpha0 = [0 0 0];
%     [alpha,resnorm,residual,exitflag,output]= lsqcurvefit(fit_exponential,alpha0,x,y);
%     plot(x,fit_exponential(alpha,x), 'k-^');
%     sum(residual.^2)

    legend('original','logistic');
    
end

result = [result, start_end_increase];

writematrix(result,'/Users/jiesi/Library/CloudStorage/OneDrive-Personal/Onedrive/BioScreen/Bioscreen C-All together-Junjun-210928/3-Calculation using Matlab/results_jiesi/result-salt-Feb072023.csv')

%%

% name of each column: 
% 1 start time of increase stage, 2 start value of increase stage, 
% 3 first end time of increase stage, 4 first end value of increase stage, 
% 5 second end time of increase stage, 6 second end value of increase stage, 
% 7 third end time of increase stage, 8 third end value of increase stage,
% 9 expoentional coefficient (a in y=a^x)
% 10 mean of following 4 points after the end of increase stage
% 11 mean of following 9 points after the end of increase stage
% 12 Y0 from fitted exponential model
% 13 r from fitted exponential model

% 14 time of the max vale,15 max value,16 mean value of the 3 max value

% 17 start time of decrease stage, 18 end value of decrease stage, 19
% gradient of decrease stage (linear, slope), 
% 20 gradient of decrease stage (exponential, r)

% 20 start time of steady stage, 21 end value of steady stage, 
% 22 gradient of steady stage, 23 average value of steady stage, 24
% standard deviation of steady stage,

% 25 start time of rise stage, 26 end value of rise stage, 27
% gradient of rise stage,

% result = [result_of_increase_stage,time_max,time_decrease,gradient_decrease,...
%     time_steady,gradient_steady,avg_steady,sd_steady,...
%     time_rise,gradient_rise];


% fid = fopen('normal_result_V3_6_21.txt','w');
% fprintf(fid,[repmat('%f\t', 1, size(result,2)), '\n'], result');    
% fclose(fid);
