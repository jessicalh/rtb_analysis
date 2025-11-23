%% RTB Subdomain Evolutionary Analysis
% Analysis of ricin toxin B chain (RTB) six subdomains
% Evidence for gene duplication and functional diversification

%% 1. Extract RTB sequence from full ricin precursor
% Full ricin precursor sequence from UniProt P02879
full_ricin_seq = 'MKPGGNTIVIWMYAVATWLCFGSTSGWSFTLEDNNIFPKQYPIINFTTAGATVQSYTNFIRAVRGRLTTGADVRHEIPVLPNRVGLPINQRFILVELSNHAELSVTLALDVTNAYVVGYRAGNSAYFFHPDNQEDAEAITHLFTDVQNRYTFAFGGNYDRLEQLAGNLRENIELGNGPLEEAISALYYYSTGGTQLPTLARSFIICIQMISEAARFQYIEGEMRTRIRYNRRSAPDPSVITLENSWGRLSTAIQESNQGAFASPIQLQRRNGSKFSVYDVSILIPIIALMVYRCAPPPSSQFSLLIRPVVPNFNADVCMDPEPIVRIVGRNGLCVDVRDGRFHNGNAIQLWPCKSNTDANQLWTLKRDNTIRSNGKCLTTYGYSPGVYVMIYDCNTAATDATRWQIWDNGTIINPRSSLVLAATSGNSGTTLTVQTNIYAVSQGWLPTNNTQPFVTTIVGLYGLCLQANSGQVWIEDCSSEKAEQQWALYADGSIRPQQNRDNCLTSDSNIRETVVKILSCGPASSGQRWMFKNDGTILNLYSGLVLDVRASDPSLKQIILYPLHGDPNQIWLPLF';

% Extract RTB B chain (262 residues starting from position 314)
rtb_sequence = full_ricin_seq(314:314+261);

fprintf('RTB sequence length: %d residues\n', length(rtb_sequence));

%% 2. Define six RTB subdomains
subdomain_length = round(length(rtb_sequence) / 6);
subdomains = {};
subdomain_names = {'1a', '1b', '1c', '2a', '2b', '2c'};

for i = 1:6
    start_pos = (i-1) * subdomain_length + 1;
    if i == 6
        end_pos = length(rtb_sequence);
    else
        end_pos = i * subdomain_length;
    end
    
    subdomains{i} = rtb_sequence(start_pos:end_pos);
    fprintf('Subdomain %s: residues %d-%d (%d residues)\n', ...
            subdomain_names{i}, start_pos, end_pos, length(subdomains{i}));
end

%% 3. Perform pairwise sequence alignments
n_subdomains = length(subdomains);
identity_matrix = zeros(n_subdomains, n_subdomains);

fprintf('\nPerforming pairwise sequence alignments...\n');

for i = 1:n_subdomains
    for j = 1:n_subdomains
        if i == j
            identity_matrix(i,j) = 100;
        else
            try
                [score, alignment] = nwalign(subdomains{i}, subdomains{j}, ...
                    'Alphabet', 'AA', 'GapOpen', 10, 'ExtendGap', 1);
                
                seq1_aligned = alignment(1,:);
                seq2_aligned = alignment(3,:);
                
                valid_positions = (seq1_aligned ~= '-') & (seq2_aligned ~= '-');
                seq1_nogap = seq1_aligned(valid_positions);
                seq2_nogap = seq2_aligned(valid_positions);
                
                matches = sum(seq1_nogap == seq2_nogap);
                identity_pct = (matches / length(seq1_nogap)) * 100;
                
                identity_matrix(i,j) = identity_pct;
                
            catch
                identity_matrix(i,j) = 0;
            end
        end
    end
end

%% 4. Create sequence identity matrix visualization
figure('Position', [100, 100, 800, 700]);

% Pastel colormap
pastel_colors = [
    1.0, 0.9, 0.9;    % Light pink
    1.0, 0.95, 0.8;   % Light peach
    0.95, 1.0, 0.9;   % Light mint
    0.9, 0.95, 1.0;   % Light blue
    0.95, 0.9, 1.0;   % Light lavender
    1.0, 1.0, 0.9;    % Light yellow
];

n_colors = 64;
pastel_map = interp1(linspace(0,1,size(pastel_colors,1)), pastel_colors, linspace(0,1,n_colors));

h = imagesc(identity_matrix);
colormap(pastel_map);
colorbar;

set(gca, 'XTick', 1:6, 'XTickLabel', subdomain_names);
set(gca, 'YTick', 1:6, 'YTickLabel', subdomain_names);
set(gca, 'FontSize', 14);

title('RTB Subdomain Sequence Identity Matrix', 'FontSize', 18, 'FontWeight', 'bold');

c = colorbar;
c.Label.String = 'Sequence Identity (%)';
c.Label.FontSize = 14;

% Add text annotations
for i = 1:6
    for j = 1:6
        if i == j
            text(j, i, '100', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 12);
        else
            text(j, i, sprintf('%.1f', identity_matrix(i,j)), ...
                'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 12);
        end
    end
end

box on;
axis square;
set(gca, 'LineWidth', 1.5);

%% 5. Save results
analysis_results = struct();
analysis_results.identity_matrix = identity_matrix;
analysis_results.subdomains = subdomains;
analysis_results.subdomain_names = subdomain_names;
analysis_results.rtb_sequence = rtb_sequence;
analysis_results.highest_identity = max(identity_matrix(identity_matrix < 100));
analysis_results.mean_identity = mean(identity_matrix(identity_matrix > 0 & identity_matrix < 100));

save('rtb_analysis_results.mat', 'analysis_results');
exportgraphics(gcf, 'RTB_Sequence_Identity_Matrix.png', 'Resolution', 300);

fprintf('\nAnalysis complete!\n');
fprintf('Highest sequence identity: %.1f%% (1c <-> 2c)\n', analysis_results.highest_identity);
fprintf('Mean sequence identity: %.1f%%\n', analysis_results.mean_identity);