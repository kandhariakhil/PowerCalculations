%%In this version we select individual segments as moving segments the number is not 2m but just m
%%which does not pertain to a pair moving at the same time, this was not accounted for in the first
%%version.
clear all
close all
clc

    R = (23/2)/100; %Radius of segment in m
    t = 0.375/100; %Thickness of segment in m
    L = 20/100; %Length of segment in m
    nu = 1.3; %Poisson ratio
    %i = 0; 
    j = 0; %Used for iterating through pairs of moving segments
    w = 1; %Number of waves
    
    min_segments = 3;
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
                    COT = (1/nu)*(((n*pi*R^3)/((n-w*(m+b))*2*L*t^2))+(((m+b)^4*(L/R)^3)/8))*(n/(w*m*(m+b)));
                    %(1/nu)*(((n*pi^4*R^3)/((n-w*(2*m+b))*4*L*t^2))+(((2*m+b)^4*(L/R)^3)/(384/5)))*(n/(w*m*(m+b))); ->Simply supported bending
                    %(1/nu)*(((n*pi^4*R^3)/((n-w*(2*m+b))*4*L*t^2)))*(n/(w*m*(m+b))); -> Only comrepssion factor
                    
                    %(1/nu)*(((n*pi^4*R^3)/((n-w*(2*m+b))*4*L*t^2))+(((2*m+b)^4*(L/R)^3)/8))*(n/(w*m*(m+b)));->Cantilevered bending
                    %(1/nu)*(((n*pi*R^3)/((n-w*(2*m+b))*2*L*t^2))+(((2*m+b)^4*(L/R)^3)/8))*(n/(w*m*(m+b)));->Cantilevered bending using curved beam theory
                    %(1/nu)*(((n*pi*R^3)/((n-w*(2*m+b))*2*L*t^2)))*(n/(w*m*(m+b))); -> Only comrepssion factor
                    
                    %If a minimum COT is measured note the w,m,b parameters for that segment number
                    if(COT<min_COT)
                        min_COT = COT;
                        wopt(n) = w;
                        mopt(n) = m;
                        bopt(n) = b;
                        COTopt(n) = min_COT;
                    end
                    m=m+1;
                end
            end
        end
    end
    
    %Record all COT
    [w_max,m_max,b_max] = size(COT);
    
    %If COT is measured as 0 convert it to NaN
    COT(COT == 0) = NaN;
    
    %Calculate percentage of anchoring segments
    for i = min_segments:max_segments
        ratio(i) = (wopt(i)*(mopt(i)+bopt(i)))/(i-wopt(i)*(mopt(i)+bopt(i)));
        perc_anchoring(i) = ((i-wopt(i)*(mopt(i)+bopt(i)))/i)*100;
    end
    
    %calculate bending and compression power factors
    for i = min_segments:max_segments
        bending_power_factor(i) = (1/nu)*(i*(mopt(i)+bopt(i))^4*(L/R)^3)/(8*wopt(i)*mopt(i)*(mopt(i)+bopt(i)));
        compression_power_factor(i) = (1/nu)*(pi*i^2*R^3)/(2*wopt(i)*mopt(i)*(mopt(i)+bopt(i))*(i-(wopt(i)*(mopt(i)+bopt(i))))*(L*t^2));
        max_velocity(i) = wopt(i)*mopt(i)*(mopt(i)+bopt(i));
    end
    
    %%
    %Store all w,m,b values for a particular number of segments for both compression and bending factor
    segment_num = 100;
    
    for w = 1:floor(segment_num/2)
        for b = 0:floor(segment_num/2)
            m = 1;
            while ((segment_num-w*(m+b)>0) && (2*m+2*b-(segment_num))<=0)
                COT_segment(m,b+1,w) = (1/nu)*(((n*pi*R^3)/((n-w*(m+b))*2*L*t^2))+(((m+b)^4*(L/R)^3)/8))*(n/(w*m*(m+b)));
                %(1/nu)*(((n*pi^4*R^3)/((n-w*(2*m+b))*4*L*t^2))+(((2*m+b)^4*(L/R)^3)/(384/5)))*(n/(w*m*(m+b)));
                %(1/nu)*(((n*pi^4*R^3)/((n-w*(2*m+b))*4*L*t^2)))*(n/(w*m*(m+b)));
                %(1/nu)*(((n*pi^4*R^3)/((n-w*(2*m+b))*4*L*t^2))+(((2*m+b)^4*(L/R)^3)/8))*(n/(w*m*(m+b)));
                
                %(1/nu)*(((n*pi*R^3)/((n-w*(2*m+b))*2*L*t^2))+(((2*m+b)^4*(L/R)^3)/8))*(n/(w*m*(m+b)));->Cantilevered bending using curved beam theory
                %(1/nu)*(((n*pi*R^3)/((n-w*(2*m+b))*2*L*t^2)))*(n/(w*m*(m+b))); -> Only comrepssion factor
                m=m+1;
            end
        end
    end
    
    COT_segment(COT_segment == 0) = NaN;
    
    %Number of moving pairs, number of bridged segments, number of waves
    [num_rows,num_cols,num_layers] = size(COT_segment);
    
    for i = 1:num_layers
        m_w0(:,i) = COT_segment(:,1,i); % middle number for number of bridge segments, 1 will give 0 bridges, 2 will give 1 bridge
    end
    
    m_w0 = log(m_w0);
        
    figure
    imagesc(m_w0);
    colormap([1 1 1; parula(256)])
    caxis([floor(min(min(m_w0))) ceil(max(max(m_w0)))])
    colorbar
%%    
 %{   
    figure
    [X,Y] = meshgrid(1:num_cols,1:num_rows);
    for i = 1:num_layers
        surf(X,Y,COT_segment(:,:,i));
        colormap(jet)
        shading interp
        xlabel('Number of bridged segments');
        ylabel('Number of moving pairs');
        zlabel('COT');
        hold on
    end
    
    figure
    for j = 1:12
        subplot(6,2,j)
        surf(X,Y,COT_segment(:,:,j));
        colormap(jet)
        shading interp
        xlabel('Number of bridged segments');
        ylabel('Number of moving pairs');
        zlabel('COT');
    end
%}        
    %%
    
    figure
    plot(min_segments:max_segments,COTopt(min_segments:max_segments));
    xlabel('Number of segments');
    ylabel('COT');
    
    figure
    plot(min_segments:max_segments,ratio(min_segments:max_segments));
    xlabel('Number of segments');
    ylabel('Percentage of anchoring segments');
    
    figure
    subplot(3,1,1)
    plot(min_segments:max_segments,wopt(min_segments:max_segments));
    xlabel('Number of segments');
    ylabel('Optimized number of waves');
    subplot(3,1,2)
    plot(min_segments:max_segments,mopt(min_segments:max_segments));
    xlabel('Number of segments');
    ylabel('Optimized number of moving pairs');    
    subplot(3,1,3)
    plot(min_segments:max_segments,bopt(min_segments:max_segments));
    xlabel('Number of segments');
    ylabel('Optimized number of bridged segments');
    
    figure
    plot(min_segments:max_segments,bending_power_factor(min_segments:max_segments),min_segments:max_segments,compression_power_factor(min_segments:max_segments));
    xlabel('Number of segments');
    ylabel('Power factor');
    legend('Bending power factor','Compression power factor');
    
    figure
    plot(min_segments:max_segments,max_velocity(min_segments:max_segments));
    xlabel('Number of segments');
    ylabel('Maximum velocity');