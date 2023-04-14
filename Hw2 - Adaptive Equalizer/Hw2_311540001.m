
clc; clear all; close all;
%% setup
d = randi([-1  1],1,2000); %generate the random number from -1 to 1

for i = 1:length(d)
   if d(i) == 0
      d(i) = 1 ; % replace 0 value by 1
   end
end
%% Try Different W
curve = zeros(3,1);
W_arr = [2.9,3.1,3.3,3.5];
L = 11;
w_plot = zeros(L,3);
figure()
for k = 1:4
    W = W_arr(k);
    h = zeros(3,1);
    delay_t = 7;
    step_size = 0.001; %0.01 is the best step_size
    for n = 1:3
        h(n) = 1/2*(1+cos(2*pi/W*(n-2)));
    end

    u = cconv(d,h,2000);
    u = u';

    % Train and find w 


    w = zeros(L,1);


    d_delay = zeros(1,length(d)+delay_t);
    for i = 1+delay_t:length(d)
        d_delay(i) = d(i-delay_t);
    end


    epoch = round(length(d)/L)-1;
    error_arr = 0;
    for i = 1:epoch
        SGD = 0;
        total_error = 0;
        for n = (i-1)*L+1:i*L
            if i == 1
               u_temp = zeros(L,1);
               u_temp(L-n+1:L) = u(1:n); %check clearly
               err = d_delay(n) - dot(w',u_temp);
               SGD = SGD + u_temp*err;
               disp(n)
               disp(size(u_temp))
               total_error = total_error + err;
            else 

               u_temp = u(n-L+1:n,1);
               disp(n)
               disp(i)
               disp(size(u_temp))
               err = d_delay(n) - dot(transpose(w),u_temp);
               SGD = SGD + u_temp*err;
    %            total_error = cat(1,total_error,err);
               total_error = total_error + err*err;
            end
        end 
        w = w + step_size*SGD;
        ave = total_error / L;
        error_arr = cat(1,error_arr,ave);
    end
    
    w_plot(:,k) = w;
    curve(k) = plot(error_arr);
    hold on;
end
xlabel("Number of adaptation cycles,n")
ylabel("Ensemble-average square error")
legend([curve(1),curve(2),curve(3),curve(4)],["W = 2.9","W = 3.1","W = 3.3","W = 3.5"],"FontSize",20); %[2.9,3.1,3.3,3.5];

figure()
result1 = cconv (w_plot(:,1),u,2000);
plot(result1)
title("Plot of convolving u with best result of learning curve (W=2.9)");


%% Try different t_delay

curve = zeros(3,1);
t_arr = [0,7,14,20];
L = 11;
w_plot = zeros(L,3);
figure()
for k = 1:4
    W = 2.9;
    h = zeros(3,1);
    delay_t = t_arr(k);
    step_size = 0.001; %0.01 is the best step_size
    for n = 1:3
        h(n) = 1/2*(1+cos(2*pi/W*(n-2)));
    end

    u = cconv(d,h,2000);
    u = u';

    % Train and find w 


    w = zeros(L,1);


    d_delay = zeros(1,length(d)+delay_t);
    for i = 1+delay_t:length(d)
        d_delay(i) = d(i-delay_t);
    end


    epoch = round(length(d)/L)-1;
    error_arr = 0;
    for i = 1:epoch
        SGD = 0;
        total_error = 0;
        for n = (i-1)*L+1:i*L
            if i == 1
               u_temp = zeros(L,1);
               u_temp(L-n+1:L) = u(1:n); %check clearly
               err = d_delay(n) - dot(w',u_temp);
               SGD = SGD + u_temp*err;
               disp(n)
               disp(size(u_temp))
               total_error = total_error + err;
            else 

               u_temp = u(n-L+1:n,1);
               disp(n)
               disp(i)
               disp(size(u_temp))
               err = d_delay(n) - dot(transpose(w),u_temp);
               SGD = SGD + u_temp*err;
    %            total_error = cat(1,total_error,err);
               total_error = total_error + err*err;
            end
        end 
        w = w + step_size*SGD;
        ave = total_error / L;
        error_arr = cat(1,error_arr,ave);
    end

    w_plot(:,k) = w;
    curve(k) = plot(error_arr);
    hold on;
end
xlabel("Number of adaptation cycles,n")
ylabel("Ensemble-average square error")
legend([curve(1),curve(2),curve(3),curve(4)],["t delay = 0","t delay = 7","t delay = 14","t delay = 20"],"FontSize",20,"Location","Best"); %[2.9,3.1,3.3,3.5];
figure()
result1 = cconv (w_plot(:,2),u,2000);
plot(result1)
title("Plot of convolving u(n) with best result of learning curve (t=7)","FontSize",15);



%% Try different learning rate

curve = zeros(3,1);
learning_rate = [0.01,0.02,0.025];
L = 11;
w_plot = zeros(L,3);
figure()
for k = 1:3
    W = 2.9;
    h = zeros(3,1);
    delay_t = 7;
    step_size = learning_rate(k); %0.01 is the best step_size
    for n = 1:3
        h(n) = 1/2*(1+cos(2*pi/W*(n-2)));
    end

    u = cconv(d,h,2000);
    u = u';

    %Train and find w 


    w = zeros(L,1);


    d_delay = zeros(1,length(d)+delay_t);
    for i = 1+delay_t:length(d)
        d_delay(i) = d(i-delay_t);
    end


    epoch = round(length(d)/L)-1;
    error_arr = 0;
    for i = 1:epoch
        SGD = 0;
        total_error = 0;
        for n = (i-1)*L+1:i*L
            if i == 1
               u_temp = zeros(L,1);
               u_temp(L-n+1:L) = u(1:n); %check clearly
               err = d_delay(n) - dot(w',u_temp);
               SGD = SGD + u_temp*err;
               disp(n)
               disp(size(u_temp))
               total_error = total_error + err;
            else 

               u_temp = u(n-L+1:n,1);
               disp(n)
               disp(i)
               disp(size(u_temp))
               err = d_delay(n) - dot(transpose(w),u_temp);
               SGD = SGD + u_temp*err;
    %            total_error = cat(1,total_error,err);
               total_error = total_error + err*err;
            end
        end 
        w = w + step_size*SGD;
        ave = total_error / L;
        error_arr = cat(1,error_arr,ave);
    end

    w_plot(:,k) = w;
    curve(k) = plot(error_arr);
    hold on;
end
xlabel("Number of adaptation cycles,n")
ylabel("Ensemble-average square error")
legend([curve(1),curve(2),curve(3)],["step size = 0.01","step size = 0.02","step size = 0.025"],"FontSize",20,"Location","Best"); %[2.9,3.1,3.3,3.5];

figure()
result1 = cconv (w_plot(:,2),u,2000);
plot(result1)
title("Plot of convolving u(n) with system that has best result of learning curve (step size = 0.05)","FontSize",15);

for i = 1:3
    figure()
    stem(w_plot(:,i));
    name = sprintf("Impulse response of system with step size = %f",learning_rate(i));
    title(name);
end 

%% Try with corrupted noise v(n)

curve = zeros(3,1);
learning_rate = [0.01,0.02,0.015];
L = 11;
w_plot = zeros(L,3);
mean = 0;
variance = 0.001;   
figure()
for k = 1:3
    W = 2.9;
    h = zeros(3,1);
    delay_t = 9;
    step_size = learning_rate(k); %0.01 is the best step_size
    for n = 1:3
        h(n) = 1/2*(1+cos(2*pi/W*(n-2)));
    end

    u = cconv(d,h,2000);
    u = u';
    u = imnoise(u,"gaussian",mean,variance);

    %Train and find w 


    w = zeros(L,1);


    d_delay = zeros(1,length(d)+delay_t);
    for i = 1+delay_t:length(d)
        d_delay(i) = d(i-delay_t);
    end


    epoch = round(length(d)/L)-1;
    error_arr = 0;
    for i = 1:epoch-50
        SGD = 0;
        total_error = 0;
        for n = (i-1)*L+1:i*L
            if i == 1
               u_temp = zeros(L,1);
               u_temp(L-n+1:L) = u(1:n); %check clearly
               err = d_delay(n) - dot(w',u_temp);
               SGD = SGD + u_temp*err;
               disp(n)
               disp(size(u_temp))
               total_error = total_error + err;
            else 

               u_temp = u(n-L+1:n,1);
               disp(n)
               disp(i)
               disp(size(u_temp))
               err = d_delay(n) - dot(transpose(w),u_temp);
               SGD = SGD + u_temp*err;
    %            total_error = cat(1,total_error,err);
               total_error = total_error + err*err;
            end
        end 
        w = w + step_size*SGD;
        ave = total_error / L;
        error_arr = cat(1,error_arr,ave);
    end

    w_plot(:,k) = w;
    curve(k) = plot(error_arr);
    hold on;
end
xlabel("Number of adaptation cycles,n")
ylabel("Ensemble-average square error")
string = ["","",""];
for i = 1:3
    string(i) = sprintf("step size = %0.05f",learning_rate(i));
end
title("Experiment with adding noise to u(n)","FontSize",15)
legend([curve(1),curve(2),curve(3)],[string(1),string(2),string(3)],"FontSize",20,"Location","Best"); %[2.9,3.1,3.3,3.5];

figure()
result1 = cconv (w_plot(:,1),u,2000);
plot(result1)
title("Plot of convolving u(n) with system that has best result of learning curve (step size = 0.01)","FontSize",15);

% for i = 1:3
%     figure()
%     stem(w_plot(:,i));
% end

%% 
figure()
stem(h)
title("h(n)")
