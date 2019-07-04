clear all
close all
clc

    R = 23/2; %Radius of segment
    t = 0.3175;
    L = 20; %Length of segment
    nu = 1.3;
    i = 0;
    j = 0;
    w = 1;
    
    max_segments = 150;
    
    for n = 6:max_segments
        min_COT = inf;
        for w = 1:floor(n/2)
            for b = 0:floor(n/2)
                m = 1;
                j = 1;
                while ((n-w*(2*m+b)>0) && (4*m+2*b-(n))<=0)
                    COT = (1/nu)*(((n*pi^4*R^3)/((n-w*(2*m+b))*4*L*t^2))+(((2*m+b)^4*(L/R)^3)/8))*(n/(w*m*(m+b)));
                    %(1/nu)*(((n*pi^4*R^3)/((n-w*(2*m+b))*4*L*t^2))+(((2*m+b)^4*(L/R)^3)/(384/5)))*(n/(w*m*(m+b)));
                    %(1/nu)*(((n*pi^4*R^3)/((n-w*(2*m+b))*4*L*t^2)))*(n/(w*m*(m+b)));
                    %(1/nu)*(((n*pi^4*R^3)/((n-w*(2*m+b))*4*L*t^2))+(((2*m+b)^4*(L/R)^3)/8))*(n/(w*m*(m+b)));
                    
                    if(COT<min_COT)
                        min_COT = COT;
                        wopt(n) = w;
                        mopt(n) = m;
                        bopt(n) = b;
                        COTopt(n) = min_COT;
                    end
                    m=m+1;
                    j=j+1;
                end
            end
        end
    end

    [w_max,m_max,b_max] = size(COT);
    
    COT(COT == 0) = NaN;
    
%     for i = 1:w_max
%         for j = 1:m_max
%             for k = 1:b_max
%                 if (COT(i,j,k) == 0)
%                     COT(i,j,k) = NaN;
%                 end
%             end
%         end
%     end
    
    for i = 6:max_segments
        perc_anchoring(i) = ((i-wopt(i)*(2*mopt(i)+bopt(i)))/i)*100;
    end
    
    for i = 6:max_segments
        bending_power_factor(i) = (i*(2*mopt(i)+bopt(i))^4*(L/R)^3)/(8*wopt(i)*mopt(i)*(mopt(i)+bopt(i)));
        compression_power_factor(i) = (pi^4*i^2*R^3)/(4*wopt(i)*mopt(i)*(mopt(i)+bopt(i))*(i-(wopt(i)*(2*mopt(i)+bopt(i))))*(L*t^2));
        max_velocity(i) = wopt(i)*mopt(i)*(mopt(i)+bopt(i));
    end
    
    %%
    %Store all values for a particular segment
    
    segment_num = 50;
    for w = 1:floor(segment_num/2)
        for b = 0:floor(segment_num/2)
            m = 1;
            j = 1;
            while ((segment_num-w*(2*m+b)>0) && (4*m+2*b-(segment_num))<=0)
                COT_segment(j,b+1,w) = (1/nu)*(((segment_num*pi^4*R^3)/((segment_num-w*(2*m+b))*4*L*t^2))+(((2*m+b)^4*(L/R)^3)/8))*(segment_num/(w*m*(m+b)));
                %(1/nu)*(((n*pi^4*R^3)/((n-w*(2*m+b))*4*L*t^2))+(((2*m+b)^4*(L/R)^3)/(384/5)))*(n/(w*m*(m+b)));
                %(1/nu)*(((n*pi^4*R^3)/((n-w*(2*m+b))*4*L*t^2)))*(n/(w*m*(m+b)));
                %(1/nu)*(((n*pi^4*R^3)/((n-w*(2*m+b))*4*L*t^2))+(((2*m+b)^4*(L/R)^3)/8))*(n/(w*m*(m+b)));
                m=m+1;
                j=j+1;
            end
        end
    end
    
    COT_segment(COT_segment == 0) = NaN;
    
    [num_rows,num_cols,num_layers] = size(COT_segment);
    
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