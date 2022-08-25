function plot_results(t_serialdate,lags,Rxy_envelope_ALL,measure,Tracks, parameters)

switch parameters.signal_type
    case 'both'

        Rxy_envelope_ALL_clicks=Rxy_envelope_ALL{1};
        Rxy_envelope_ALL_whistles=Rxy_envelope_ALL{2};

        %% 1) PLOT CROSS-CORRELOGRAM:

        figure;
        ax1 = axes;
        im = imagesc(t_serialdate,lags,Rxy_envelope_ALL_clicks); datetick('x','keeplimits');
        colormap(flipud(gray(256)))
        ylim([-parameters.d/parameters.c,parameters.d/parameters.c])
        xlabel('Local Time (HH:MM:SS)'), ylabel('TDOA (s)'),
        title(['Cross-correlogram based on ', parameters.signal_type])
        set(gca,'FontSize',12)
        caxis([0,10])
        im.AlphaData = 0.5; % change this value to change the background image transparency
%         axis square;
        hold all;
        %plot second data
        ax2 = axes;
        im1 = imagesc(t_serialdate,lags,Rxy_envelope_ALL_whistles); datetick('x','keeplimits');
        colormap(flipud(gray(256)))
        colorbar
        ylim([-parameters.d/parameters.c,parameters.d/parameters.c])
        xlabel('Local Time (HH:MM:SS)'), ylabel('TDOA (s)'),
        set(gca,'FontSize',12)
        caxis([0,10])
        im1.AlphaData = 0.5; % change this value to change the foreground image transparency
%         axis square;
        %link axes
        linkaxes([ax1,ax2])
        %%Hide the top axes
        ax2.Visible = 'off';
        ax2.XTick = [];
        ax2.YTick = [];
        set([ax1,ax2],'Position',[.17 .11 .685 .815]);
        colorbar(ax1,'Position',[.05 .11 .03 .815]);
        colorbar(ax2,'Position',[.88 .11 .03 .815]);
%         set(findall(gcf,'type','text'),'FontSize',14)
        

        %% 2) PLOT TRACKED TDOAS

        figure;

        %plot cross-correlograms (overlayed)
        ax1 = axes;
        im = imagesc(t_serialdate,lags,Rxy_envelope_ALL_clicks); datetick('x','keeplimits');
        colormap(flipud(gray(256)))
        ylim([-parameters.d/parameters.c,parameters.d/parameters.c])
        xlabel('Local Time (HH:MM:SS)'), ylabel('TDOA (s)'),
        title(['Cross-correlogram based on ', parameters.signal_type])
        set(gca,'FontSize',12)
        caxis([0,10])
        im.AlphaData = 0.5; % change this value to change the background image transparency
%         axis square;
        hold all;
        %plot second data
        ax2 = axes;
        im1 = imagesc(t_serialdate,lags,Rxy_envelope_ALL_whistles); datetick('x','keeplimits');
        colormap(flipud(gray(256)))
        colorbar
        ylim([-parameters.d/parameters.c,parameters.d/parameters.c])
        xlabel('Local Time (HH:MM:SS)'), ylabel('TDOA (s)'),
        set(gca,'FontSize',12)
        caxis([0,10])
        im1.AlphaData = 0.5; % change this value to change the foreground image transparency
%         axis square;
        %link axes
        linkaxes([ax1,ax2])
        %%Hide the top axes
        ax2.Visible = 'off';
        ax2.XTick = [];
        ax2.YTick = [];
        set([ax1,ax2],'Position',[.17 .11 .685 .815]);
        colorbar(ax1,'Position',[.05 .11 .03 .815]);
        colorbar(ax2,'Position',[.88 .11 .03 .815]);
       
        
        %plot measurements
%         hold on
%         for k=1:measure.T
%             if ~isempty(measure.Z{k})
%                 plot(t_serialdate(k),measure.Z{k}(1,:),'ro')
%             end
%         end

        %plot tracked TDOAs
        hold on
        for k=1:size(Tracks,2)
            plot(Tracks(k).time_local, Tracks(k).tdoa,'-','LineWidth',3)
        end
        datetick('x','keeplimits');

%         h(1) = plot(NaN, NaN,'ro','LineWidth',2);
%         h(2) = plot(NaN, NaN,'b-','LineWidth',2);
%         legend(h,'Measurements','Tracked TDOAs','Location', 'southeastoutside');

        h(1) = plot(NaN, NaN,'b-','LineWidth',2);
        legend(h,'Tracked TDOAs','Location', 'southeastoutside');

        set(gca,'FontSize',14)
        set(findall(gcf,'type','text'),'FontSize',14)

        hold off



    otherwise

        %% 1) PLOT CROSS-CORRELOGRAM:

        f1= figure;
        imagesc(t_serialdate,lags,Rxy_envelope_ALL), datetick('x','keeplimits');
        colormap(flipud(gray(256)))
        colorbar
        ylim([-parameters.d/parameters.c,parameters.d/parameters.c])
        xlabel('Local Time (HH:MM:SS)'), ylabel('TDOA (s)'),
        title(['Cross-correlogram based on ', parameters.signal_type])
        set(gca,'FontSize',14)
        set(findall(gcf,'type','text'),'FontSize',14)
        f1.PaperUnits='centimeters';
        f1.PaperPosition=[0,0,100,30];%size specified as panorama, 5:1- [left,bottom,width, height]

        %% 2) PLOT EXTRACTED MEASUREMENTS


        figure;
        %plot cross-correlogram
        imagesc(t_serialdate,lags,Rxy_envelope_ALL), datetick('x','keeplimits');
        colormap(flipud(gray(256)))
        colorbar
        ylim([-parameters.d/parameters.c,parameters.d/parameters.c])
        xlabel('Local Time (HH:MM:SS)'), ylabel('TDOA (s)'),
        title(['Cross-correlogram based on ', parameters.signal_type])

        %plot measurements
        hold on
        for k=1:measure.T
            if ~isempty(measure.Z{k})
                plot(t_serialdate(k),measure.Z{k}(1,:),'ro')
            end
        end

        h(1) = plot(NaN, NaN,'ro','LineWidth',2);
        legend(h,'Measurements','Location', 'EastOutside');

        set(gca,'FontSize',14)
        set(findall(gcf,'type','text'),'FontSize',14)

        hold off

        %% 3) PLOT TRACKED TDOAS

        figure;
        %plot cross-correlogram
        imagesc(t_serialdate,lags,Rxy_envelope_ALL), datetick('x','keeplimits');
        colormap(flipud(gray(256)))
        colorbar
        ylim([-parameters.d/parameters.c,parameters.d/parameters.c])
        xlabel('Local Time (HH:MM:SS)'), ylabel('TDOA (s)'),
        title(['Tracked TDOAs based on ', parameters.signal_type])

        %plot measurements
        hold on
        for k=1:measure.T
            if ~isempty(measure.Z{k})
                plot(t_serialdate(k),measure.Z{k}(1,:),'ro')
            end
        end

        %plot tracked TDOAs
        for k=1:size(Tracks,2)
            plot(Tracks(k).time_local, Tracks(k).tdoa,'-','LineWidth',3)
        end
        datetick('x','keeplimits');

        h(1) = plot(NaN, NaN,'ro','LineWidth',2);
        h(2) = plot(NaN, NaN,'b-','LineWidth',2);
        legend(h,'Measurements','Tracked TDOAs','Location', 'EastOutside');

        set(gca,'FontSize',14)
        set(findall(gcf,'type','text'),'FontSize',14)

        hold off

end
end