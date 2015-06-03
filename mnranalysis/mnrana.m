clear all
close all

matlabpool open 4
files = dir('roi*.mat');

c = 1;
for fileval = 1:19
    clear smg pair pair2 pwr pwr2 H gErr tgraph mtx
    c = 1;
    for i = 1:42
        temp = files(i + 42 * (fileval - 1)).name;
        mtx = load(temp); % log10
        tgraph = mtx.roi_data;
        tgraph(isinf(tgraph))=0;
        %     tgraph=full(temp.graph);
        smg(:,:,c) = tgraph;

        c = c+1;
    end
    numrois = size(smg,1);
    roinum(fileval) = numrois;
    ID = [1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 10 10 11 11 12 12 13 13 14 14 15 15 16 16 17 17 18 18 19 19 20 20 21 21];
    % compute the transform of each to obtain frequency spectrum
    L = 210; % length of the time signal
    NFFT = 2^nextpow2(L); 
    Fs = 1;
    T = 1/Fs;
    t = (0:L-1)*T;
    w = rectwin(L);
    pwr = zeros(numrois, NFFT/2+1, 42);
    parfor c = 1:42
        for i = 1:numrois
     % fourier transform of the timeseries
            pwr(i,:,c) = pwelch(smg(i,:,c), w, 0, NFFT);
            pwr(i,:,c) = pwr(i,:,c)/sum(pwr(i,:,c)); % norm by sum    
        end
    end
    pwr(isnan(pwr)) = 0;
    pair = zeros(42, numrois, numrois);
    parfor i = 1:42
        pair(i,:,:) = corr(smg(:,:,i)');
    end
    pair(isnan(pair)) = 0;
    numiter = 1;
    tvalpwr = logspace(-4, 0, 25);
    tvalcorr = linspace(min(min(min(pair))), 1, 25);
    % tvalpwr = 0;
    % tvalcorr = 0;
    for tloc = 1:length(tvalpwr)
        %threshold the scans
        pwr2 = pwr;
        pwr2(pwr2<tvalpwr(tloc)) = 0; 
        %renormalize
        for i = 1:42
            for j = 1:numrois
                sumpwr = sum(pwr2(j,:,i));
                if sumpwr == 0
                    pwr2(j,:,i) = zeros(length(pwr2(j,:,i)),1);
                else
                    pwr2(j,:,i) = pwr2(j,:,i)/sum(pwr2(j,:,i));
                end    
            end
        end
        %use the thresholded power for the scans
        H = zeros(numrois, 42, 42);
        parfor j = 1:numrois
            for a = 1:42
                for b = 1:42
                    %j is the roi number, i is the scan, z is the scan being compared
                    %to
                    H(j,a,b) = norm(sqrt(pwr2(j,:,a)) - sqrt(pwr2(j,:,b)))/sqrt(2);

                end
            end
        end
    %     for a = 1:42
    %         for b = 1:42
    %             g(:,:) = pwr2(:,:,a) - pwr2(:,:,b);
    %             gErr(a,b) = norm(g, 'fro');
    %     
    %         end
    %     end
        gErr = zeros(42, 42);
        parfor i = 1:42
            for j = 1:42
                gErr(i,j) = (1/numrois)*sum(H(:,i,j));
            end
        end

    % %% Compute TRT
    %     figure(1), imagesc(log(gErr)), colorbar
    %     matches = 0;
    %     intra = 0; inter = 0;
    %     intravec = [];
    %     intervec = [];
    %     for i = 1:1:size(gErr,1)
    %         [temp,ind] = sort(gErr(i,:));
    %         q = i-1+2*mod(i,2);
    %         z = 2;
    %         lowval = temp(z);
    %         numvals = 0;
    %         find = 0;
    %         while lowval == temp(2) 
    %             numvals = numvals + 1;
    %             if q == ind(z)
    %                 find = 1;
    %             end
    %             z = z + 1;
    %             if z>length(gErr)
    %                 break
    %             end
    %             lowval = temp(z);
    %         end
    %         matchadd = find/numvals;
    %         matches = matches + matchadd;
    %         intra = intra + gErr(i,q);
    %         temp2 = sum([gErr(i,1:q-1), gErr(i,q+1:end)])/40; %40 bc 1 is intra, 1 is me
    %         inter = inter + temp2;
    %         intravec = [intravec, gErr(i,q)];
    %         if i < q
    %             if i > 1
    %                 intervec = [intervec, gErr(i,q-2), gErr(i,q+1:end)];
    %             else
    %                 intervec = [intervec, gErr(i,q+1:end)];
    %             end
    %         else
    %             if i > 2
    %                 intervec = [intervec, gErr(i,q-1), gErr(i,q+2:end)];
    %             else
    %                 intervec = [intervec, gErr(i,q+2:end)];
    %             end
    %         end
    %     end
    %     intra = intra/42;
    %     inter = inter/42;
        TRTpwr(tloc) = compute_mnr(gErr, ID);
    %     title(strcat('correct matches=', num2str(matches),'/', num2str(size(smg,3))));

        % compute the pairwise correlation and threshold for sake of comparison
        pair2 = pair;
        pair2(pair2<tvalcorr(tloc)) = 0;
        gErr = zeros(42, 42);
        parfor i = 1:42
            for j = 1:42
                gErr(i,j) = sqrt(sum(sum((pair2(i,:,:) - pair2(j,:,:)).^2)));
            end
        end

    %     %% Compute TRT
    %     figure(2), imagesc(log(gErr)), colorbar
    %     matches = 0;
    %     intra = 0; inter = 0;
    %     intravec = [];
    %     intervec = [];
    %     for i = 1:1:size(gErr,1)
    %         [temp,ind] = sort(gErr(i,:));
    %         q = i-1+2*mod(i,2);
    %         z = 2;
    %         lowval = temp(z);
    %         numvals = 0;
    %         find = 0;
    %         while lowval == temp(2) 
    %             numvals = numvals + 1;
    %             if q == ind(z)
    %                 find = 1;
    %             end
    %             z = z + 1;
    %             if z>length(gErr)
    %                 break
    %             end
    %             lowval = temp(z);
    %         end
    %         matchadd = find/numvals;
    %         matches = matches + matchadd;
    %         intra = intra + gErr(i,q);
    %         temp2 = sum([gErr(i,1:q-1), gErr(i,q+1:end)])/40; %40 bc 1 is intra, 1 is me
    %         inter = inter + temp2;
    %         intravec = [intravec, gErr(i,q)];
    %         if i < q
    %             if i > 1
    %                 intervec = [intervec, gErr(i,q-2), gErr(i,q+1:end)];
    %             else
    %                 intervec = [intervec, gErr(i,q+1:end)];
    %             end
    %         else
    %             if i > 2
    %                 intervec = [intervec, gErr(i,q-1), gErr(i,q+2:end)];
    %             else
    %                 intervec = [intervec, gErr(i,q+2:end)];
    %             end
    %         end
    %     end
    %     intra = intra/42;
    %     inter = inter/42;
    %     TRTcorr(tloc) = matches/size(gErr, 1);
        TRTcorr(tloc) = compute_mnr(gErr, ID);
    %     title(strcat('correct matches=', num2str(matches),'/', num2str(size(smg,3))));
    end
    subplot(1,2,1)
    plot(tvalpwr, TRTpwr);
    set(gca, 'xscale', 'log')
    title('Threshold vs TRT for power spectrum')
    xlabel('Threshold 10^x')
    subplot(1,2,2);
    plot(tvalcorr, TRTcorr);
    title('Threshold vs TRT for correlation')
    xlabel('Threshold')
    mincorr(fileval) = min(TRTcorr);
    minpwr(fileval) = min(TRTpwr);
    save(sprintf('v2iteration%d.mat', fileval));
    fprintf(sprintf('Done round %d\n', fileval));
end
matlabpool close
system('shutdown -s')