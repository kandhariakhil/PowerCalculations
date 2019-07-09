clear all
close all
clc

    R = 23/2; %Radius of segment
    t = 0.375; %Thickness of segment
    L = 20; %Length of segment
    nu = 1.3; %Poisson ratio
    %i = 0; 
    j = 0; %Used for iterating through pairs of moving segments
    w = 1; %Number of waves
    
    min_segments = 6;
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
                while ((n-w*(2*m+b)>0) && (4*m+2*b-(n))<=0)
                    %Calculating COT based on Math_ver9.
                    COT = (1/nu)*(((n*pi*R^3)/((n-w*(2*m+b))*2*L*t^2)))*(n/(w*m*(m+b)));
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
    for i = 6:max_segments
        perc_anchoring(i) = ((i-wopt(i)*(2*mopt(i)+bopt(i)))/i)*100;
    end
    
    %calculate bending and compression power factors
    for i = 6:max_segments
        bending_power_factor(i) = (i*(2*mopt(i)+bopt(i))^4*(L/R)^3)/(8*wopt(i)*mopt(i)*(mopt(i)+bopt(i)));
        compression_power_factor(i) = (pi^4*i^2*R^3)/(4*wopt(i)*mopt(i)*(mopt(i)+bopt(i))*(i-(wopt(i)*(2*mopt(i)+bopt(i))))*(L*t^2));
        max_velocity(i) = wopt(i)*mopt(i)*(mopt(i)+bopt(i));
    end
    
    %%
    %Store all w,m,b values for a particular number of segments for both compression and bending factor
    segment_num = 6;
    
    for w = 1:floor(segment_num/2)
        for b = 0:floor(segment_num/2)
            m = 1;
            while ((segment_num-w*(2*m+b)>0) && (4*m+2*b-(segment_num))<=0)
                COT_segment(m,b+1,w) = (1/nu)*(((n*pi*R^3)/((n-w*(2*m+b))*2*L*t^2)))*(n/(w*m*(m+b)));
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
        m_w0(:,i) = COT_segment(:,1,i); %0 bridge segments
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
    plot(6:max_segments,COTopt(6:max_segments));
    xlabel('Number of segments');
    ylabel('COT');
    
    figure
    plot(6:max_segments,perc_anchoring(6:max_segments));
    xlabel('Number of segments');
    ylabel('Percentage of anchoring segments');
    
    figure
    subplot(3,1,1)
    plot(6:max_segments,wopt(6:max_segments));
    xlabel('Number of segments');
    ylabel('Optimized number of waves');
    subplot(3,1,2)
    plot(6:max_segments,mopt(6:max_segments));
    xlabel('Number of segments');
    ylabel('Optimized number of moving pairs');    
    subplot(3,1,3)
    plot(6:max_segments,bopt(6:max_segments));
    xlabel('Number of segments');
    ylabel('Optimized number of bridged segments');
    
    figure
    plot(6:max_segments,bending_power_factor(6:max_segments),6:max_segments,compression_power_factor(6:max_segments));
    xlabel('Number of segments');
    ylabel('Power factor');
    legend('Bending power factor','Compression power factor');
    
    figure
    plot(6:max_segments,max_velocity(6:max_segments));
    xlabel('Number of segments');
    ylabel('Maximum velocity');