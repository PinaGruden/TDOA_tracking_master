function plot_results(t_serialdate,lags,Rxy_envelope_ALL,measure,Tracks, parameters)

switch parameters.signal_type
    case 'both'

        Rxy_envelope_ALL_clicks=Rxy_envelope_ALL{1};
        Rxy_envelope_ALL_whistles=Rxy_envelope_ALL{2};

        
        
        %% 1) PLOT Extracted MEASUREMENTS

        f1=figure(1);

        %plot whistle& click cross-correlograms (overlayed)
        ax1 = axes(f1);
        im = imagesc(t_serialdate,lags,Rxy_envelope_ALL_clicks); 
        datetick('x','keeplimits');
        colormap(flipud(gray(256)))
        ylim([-parameters.d/parameters.c,parameters.d/parameters.c])
        xlabel('Local Time (HH:MM:SS)'), ylabel('TDOA (s)'),
        title(['Extracted Measurements from Cross-correlogram based on ', parameters.signal_type])
        set(gca,'FontSize',14)
        caxis([0,10])
        im.AlphaData = 0.5; % change this value to change the background image transparency
        hold all;
        %plot second data
        ax2 = axes(f1);
        im1 = imagesc(t_serialdate,lags,Rxy_envelope_ALL_whistles);
        datetick('x','keeplimits');
        colormap(flipud(gray(256)))
        ylim([-parameters.d/parameters.c,parameters.d/parameters.c])
        xlabel('Local Time (HH:MM:SS)'), ylabel('TDOA (s)'),
        set(gca,'FontSize',14)
        caxis([0,10])
        im1.AlphaData = 0.5; % change this value to change the foreground image transparency
        %link axes
        linkaxes([ax1,ax2])
        %%Hide the top axes
        ax2.Visible = 'off';
        ax2.XTick = [];
        ax2.YTick = [];
        set([ax1,ax2],'Position',[.17 .11 .685 .815]);
%         colorbar(ax1,'Position',[.05 .11 .03 .815]);
%         colorbar(ax2,'Position',[.88 .11 .03 .815]);
       
        %plot measurements
        hold on
        for k=1:measure.T
            if ~isempty(measure.Z{k})
                plot(t_serialdate(k),measure.Z{k}(1,:),'b.'),hold on
            end
        end
        datetick('x','keeplimits');

        h(1) = plot(NaN, NaN,'b.','LineWidth',4);
        legend(h,'Measurements','Location', 'southeast');

        set(gca,'FontSize',14)
        set(findall(gcf,'type','text'),'FontSize',14)

        hold off

        %% 2) PLOT TRACKED TDOAS

        f2=figure(2);

        %plot cross-correlograms (overlayed)
        ax1 = axes(f2);
        im = imagesc(t_serialdate,lags,Rxy_envelope_ALL_clicks); 
        datetick('x','keeplimits');
        colormap(flipud(gray(256)))
        ylim([-parameters.d/parameters.c,parameters.d/parameters.c])
        xlabel('Local Time (HH:MM:SS)'), ylabel('TDOA (s)'),
        title(['Tracked TDOAs from measurements based on ', parameters.signal_type])
        set(gca,'FontSize',14)
        caxis([0,10])
        im.AlphaData = 0.5; % change this value to change the background image transparency
        hold all;
        %plot second data
        ax2 = axes(f2);
        im1 = imagesc(t_serialdate,lags,Rxy_envelope_ALL_whistles); 
        datetick('x','keeplimits');
        colormap(flipud(gray(256)))
        ylim([-parameters.d/parameters.c,parameters.d/parameters.c])
        xlabel('Local Time (HH:MM:SS)'), ylabel('TDOA (s)'),
        set(gca,'FontSize',14)
        caxis([0,10])
        im1.AlphaData = 0.5; % change this value to change the foreground image transparency
        %link axes
        linkaxes([ax1,ax2])
        %%Hide the top axes
        ax2.Visible = 'off';
        ax2.XTick = [];
        ax2.YTick = [];
        set([ax1,ax2],'Position',[.17 .11 .685 .815]);
%         colorbar(ax1,'Position',[.05 .11 .03 .815]);
%         colorbar(ax2,'Position',[.88 .11 .03 .815]);

        %plot tracked TDOAs
        hold on
        for k=1:size(Tracks,2)
            plot(ax2,Tracks(k).time_local, Tracks(k).tdoa,'-','LineWidth',3), hold on
        end
        datetick('x','keeplimits');

        h(1) = plot(NaN, NaN,'b-','LineWidth',3);
        legend(h,'Tracked TDOAs','Location', 'southeast');

        set(gca,'FontSize',14)
        set(findall(gcf,'type','text'),'FontSize',14)

        hold off

%% 3) PLOT CROSS-CORRELOGRAMS:

        figure(3);

        subplot(211) %cross-correlogram for clicks
        imagesc(t_serialdate,lags,Rxy_envelope_ALL_clicks); 
        datetick('x','keeplimits');
        colormap(flipud(gray(256)));
        ylim([-parameters.d/parameters.c,parameters.d/parameters.c]);
        xlabel('Local Time (HH:MM:SS)'); ylabel('TDOA (s)');
        title('Cross-correlogram based on clicks');
        set(gca,'FontSize',14);
        caxis([0,10]);
        colorbar;

        subplot(212) %cross-correlogram for whistles
        imagesc(t_serialdate,lags,Rxy_envelope_ALL_whistles); 
        datetick('x','keeplimits');
        colormap(flipud(gray(256)));
        ylim([-parameters.d/parameters.c,parameters.d/parameters.c]);
        xlabel('Local Time (HH:MM:SS)'); ylabel('TDOA (s)');
        title(['Cross-correlogram based on whistles']);
        set(gca,'FontSize',14);
        caxis([0,10]);
        colorbar;


    otherwise

        %% 1) PLOT CROSS-CORRELOGRAM:

        figure(1);
        imagesc(t_serialdate,lags,Rxy_envelope_ALL), 
        datetick('x','keeplimits');
        colormap(flipud(gray(256)))
        caxis([0,10])
        colorbar
        ylim([-parameters.d/parameters.c,parameters.d/parameters.c])
        xlabel('Local Time (HH:MM:SS)'), ylabel('TDOA (s)'),
        title(['Cross-correlogram based on ', parameters.signal_type])
        set(gca,'FontSize',14)
        set(findall(gcf,'type','text'),'FontSize',14)
        

        %% 2) PLOT EXTRACTED MEASUREMENTS

        figure(2);
        %plot cross-correlogram
        imagesc(t_serialdate,lags,Rxy_envelope_ALL), 
        datetick('x','keeplimits');
        colormap(flipud(gray(256)))
        caxis([0,10])
        colorbar
        ylim([-parameters.d/parameters.c,parameters.d/parameters.c])
        xlabel('Local Time (HH:MM:SS)'), ylabel('TDOA (s)'),
        title(['Extracted Measurements from Cross-correlogram based on ', parameters.signal_type])

        %plot measurements
        hold on
        for k=1:measure.T
            if ~isempty(measure.Z{k})
                plot(t_serialdate(k),measure.Z{k}(1,:),'b.')
            end
        end

        h(1) = plot(NaN, NaN,'b.','LineWidth',4);
        legend(h,'Measurements','Location', 'southeast');

        set(gca,'FontSize',14)
        set(findall(gcf,'type','text'),'FontSize',14)

        hold off

        %% 3) PLOT TRACKED TDOAS

        figure(3);
        %plot cross-correlogram
        imagesc(t_serialdate,lags,Rxy_envelope_ALL), 
        datetick('x','keeplimits');
        colormap(flipud(gray(256)))
        caxis([0,10])
        colorbar
        ylim([-parameters.d/parameters.c,parameters.d/parameters.c])
        xlabel('Local Time (HH:MM:SS)'), ylabel('TDOA (s)'),
        title(['Tracked TDOAs from measurements based on ', parameters.signal_type])

        %plot tracked TDOAs
        hold on
        for k=1:size(Tracks,2)
            plot(Tracks(k).time_local, Tracks(k).tdoa,'-','LineWidth',3)
        end
        datetick('x','keeplimits');
 
        h(1) = plot(NaN, NaN,'b-','LineWidth',3);
        legend(h,'Tracked TDOAs','Location', 'southeast');

        set(gca,'FontSize',14)
        set(findall(gcf,'type','text'),'FontSize',14)

        hold off

end
end