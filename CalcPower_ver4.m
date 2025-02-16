%%In this version we select individual segments as moving segments the number is not 2m but just m
%%which does not pertain to a pair moving at the same time, this was not accounted for in the first
%%version.
clear all
close all
clc

    R = (23/2)/100; %Radius of segment in m
    t = 0.375/100; %Thickness of segment in m
    L = 20/100; %Length of segment in m
    nu = 1.4; %Poisson ratio 
    j = 0; %Used for iterating through pairs of moving segments
    w = 1; %Number of waves
    
    min_segments = 4; %Minimum number of segments
    max_segments = 150; %Maximum number of segments
    
    %Iterating number of segments from min to max
    for n = min_segments:max_segments
        min_COT = inf; %Selecting a minimum COT for each iteration
        for w = 1:floor(n/2) %Maximum number of waves possible is n/2 waves
            for b = 0:floor(n/2) %Maximum number of bridged segments possible is half the number of segments
                m = 1;
                %While loop to iterate thorough all possible combinations in integers given
                %constraints are met with, while loop includes the 2 possible constraints set up for
                %this equation
                
                %fmincon was used at first but would not meet all constraints and would not give
                %only integer values due to which this optimization was created
                while ((n-w*(m+b)>0) && (2*m+2*b-(n))<=0)
                    %Calculating COT based on Math_ver9.
                    %COT = (1/nu)*(((n*pi*R^3)/((100)*((n-w*(m+b))*2*L*t^2)))+(((m+b)^4*(L/R)^3)/8))*(n/(w*(m/2)*((m/2)+b)));    
                    COT = (1/nu)*(((n*pi*R^3)/((n-w*(m+b))*2*L*t^2))+(((m+b)^4*(L/R)^3)/8))*(n/(w*(m/2)*((m/2)+b)));                    
                    %If a minimum COT is measured note the w,m,b parameters for that segment number
                    if(COT<min_COT)
                        min_COT = COT;
                        wopt(n) = w;
                        if mod(m,2)~= 0
                            mopt(n) = m-1;
                            bopt(n) = 1;
                        else
                            mopt(n) = m;
                            bopt(n) = b;
                        end
                        COTopt(n) = (1/nu)*((2*pi*n^2*R^3)/((n-wopt(n)*(mopt(n)+bopt(n)))*wopt(n)*mopt(n)*(mopt(n)+2*bopt(n))*(L*t^2))+((n*(mopt(n)+bopt(n))^4*L^3)/(2*wopt(n)*mopt(n)*(mopt(n)+2*bopt(n))*R^3)));
                        %COTopt(n) = (1/nu)*(((n*pi*R^3)/((100)*((n-wopt(n)*(mopt(n)+bopt(n)))*2*L*t^2)))+(((mopt(n)+bopt(n))^4*(L/R)^3)/8))*(n/(wopt(n)*(mopt(n)/2)*((mopt(n)/2)+bopt(n))));
                    end
                    m=m+1;
                end
            end
        end
    end
    
    for n = min_segments:max_segments
        max_vel = 0; %Selecting a maximum velocity for each iteration
        for w = 1:floor(n/2) %Maximum number of waves possible is n/2 waves
            for b = 0:floor(n/2) %Maximum number of bridged segments possible is half the number of segments
                m = 1;
                %While loop to iterate thorough all possible combinations in integers given
                %constraints are met with, while loop includes the 2 possible constraints set up for
                %this equation

                %fmincon was used at first but would not meet all constraints and would not give
                %only integer values due to which this optimization was created
                while ((n-w*(m+b)>0) && (2*m+2*b-(n))<=0)
                    %Calculating ideal velocity based on Math_ver9.
                    opt_vel = w*m*(m+b)/n;

                    %If a ideal velocity is measured note the w,m,b parameters for that segment number
                    if(opt_vel>max_vel)
                        max_vel = opt_vel;
                        vel_wopt(n) = w;
                        if mod(m,2)~= 0
                            vel_mopt(n) = m-1;
                            vel_bopt(n) = 1;
                        else
                            vel_mopt(n) = m;
                            vel_bopt(n) = b;
                        end
                        vel_opt(n) = vel_wopt(n)*vel_mopt(n)*(vel_mopt(n)+vel_bopt(n))/n;
                    end
                    m=m+1;
                end
            end
        end
    end
    
    %Calculate percentage of anchoring segments and ratio of moving to anchoring segments based on optimized velocity
    for i = min_segments:max_segments
        vel_ratio(i) = (vel_wopt(i)*(vel_mopt(i)+vel_bopt(i)))/(i-(vel_wopt(i)*(vel_mopt(i)+vel_bopt(i))));
    end    
        
   %Calculate percentage of anchoring segments and ratio of moving to anchoring segments based on optimized COT
    for i = min_segments:max_segments
        ratio(i) = (wopt(i)*(mopt(i)+bopt(i)))/(i-(wopt(i)*(mopt(i)+bopt(i))));
    end
    
    %Record all COT
    [w_max,m_max,b_max] = size(COT);
    
    %If COT is measured as 0 convert it to NaN
    COT(COT == 0) = NaN;

    for i = min_segments:max_segments
        bending_power_factor(i) = (1/nu)*(i*(mopt(i)+bopt(i))^4*(L/R)^3)/(8*wopt(i)*(mopt(i)/2)*((mopt(i)/2)+bopt(i)));
        compression_power_factor(i) = (1/nu)*(pi*i^2*R^3)/(2*wopt(i)*(mopt(i)/2)*((mopt(i)/2)+bopt(i))*(i-(wopt(i)*(mopt(i)+bopt(i))))*(L*t^2));
        max_velocity(i) = wopt(i)*(mopt(i)/2)*((mopt(i)/2)+bopt(i))/(i*L);
    end
    
    %%
    %Store all w,m,b values for a particular number of segments for both compression and bending factor
    segment_num = 100;
    
    for w = 1:floor(segment_num/2)
        m = 1;
        b = 0;
        move = 1;
        while ((segment_num-w*(m+b)>0) && (2*m+2*b-(segment_num))<=0)
            COT_segment(move,1,w) = (1/nu)*(((segment_num*pi*R^3)/((segment_num-w*(m+b))*2*L*t^2))+(((m+b)^4*(L/R)^3)/8))*(segment_num/(w*(m/2)*((m/2)+b)));
            speed_100seg(move,1,w) = (w*(m/2)*((m/2)+b))/segment_num;           
            move = move+1;
            if mod(move,2) == 0
                m = move;
                b = 0;
            else
                m = move-1;
                b = 1;
            end            
        end
    end
    
    COT_segment(1,:,:) = NaN;
    COT_segment(COT_segment == 0) = NaN;
    speed_100seg(speed_100seg == 0) = NaN;
    speed_100seg(1,:,:) = NaN;
    
    %Number of moving pairs, number of bridged segments, number of waves
    [num_rows,num_cols,num_layers] = size(COT_segment);
    
    for i = 1:num_layers
        m_w0(:,i) = COT_segment(:,1,i); % middle number for number of bridge segments, 1 will give 0 bridges, 2 will give 1 bridge
        v_100(:,i) = speed_100seg(:,1,i);
    end
    
    m_w0 = log(m_w0);
        
    figure
    imagesc(m_w0);
    colormap([1 1 1; parula(256)])
    caxis([floor(min(min(m_w0))) ceil(max(max(m_w0)))])
    colorbar
    
    figure
    imagesc(v_100);
    colormap([1 1 1; parula(256)])
    caxis([-1 ceil(max(max(v_100)))])
    colorbar
%% Plots     

    figure
    plot(min_segments:max_segments,COTopt(min_segments:max_segments));
    xlabel('Number of segments');
    ylabel('COT');
    
    figure
    plot(min_segments:max_segments,ratio(min_segments:max_segments));
    xlabel('Number of segments');
    ylabel('Percentage of anchoring segments');
    
    figure
    plot(min_segments:max_segments,vel_opt(min_segments:max_segments));
    xlabel('Number of segments');
    ylabel('Optimal velocity');
    
    figure
    plot(min_segments:max_segments,vel_ratio(min_segments:max_segments));
    xlabel('Number of segments');
    ylabel('Percentage of anchoring segments for optimal velocity');
    
    figure
    subplot(3,1,1)
    plot(min_segments:max_segments,vel_wopt(min_segments:max_segments));
    xlabel('Number of segments');
    ylabel('Optimized number of waves for velocity');
    subplot(3,1,2)
    plot(min_segments:max_segments,vel_mopt(min_segments:max_segments)+vel_bopt(min_segments:max_segments));
    xlabel('Number of segments');
    ylabel('Optimized number of moving pairs for velocity');    
    subplot(3,1,3)
    plot(min_segments:max_segments,bopt(min_segments:max_segments));
    xlabel('Number of segments');
    ylabel('Optimized number of bridged segments');    
    
    figure
    subplot(3,1,1)
    plot(min_segments:max_segments,wopt(min_segments:max_segments));
    xlabel('Number of segments');
    ylabel('Optimized number of waves');
    subplot(3,1,2)
    plot(min_segments:max_segments,mopt(min_segments:max_segments)+bopt(min_segments:max_segments));
    xlabel('Number of segments');
    ylabel('Optimized number of moving pairs');    
    subplot(3,1,3)
    plot(min_segments:max_segments,bopt(min_segments:max_segments));
    xlabel('Number of segments');
    ylabel('Optimized number of bridged segments');
    
    figure
    plot(min_segments:max_segments,bending_power_factor(min_segments:max_segments),min_segments:max_segments,compression_power_factor(min_segments:max_segments),min_segments:max_segments,COTopt(min_segments:max_segments));
    xlabel('Number of segments');
    ylabel('Power factor');
    legend('Bending power factor','Compression power factor');
    
    figure
    plot(min_segments:max_segments,max_velocity(min_segments:max_segments));
    xlabel('Number of segments');
    ylabel('Maximum velocity');