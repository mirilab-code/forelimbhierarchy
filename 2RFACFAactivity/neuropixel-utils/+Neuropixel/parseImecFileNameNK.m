function [pathRoot, fileStem, type] = parseImecFileNameNK(file)
            if iscell(file)
                [pathRoot, fileStem, type] = cellfun(@Neuropixel.ImecDatasetNK.parseImecFileName, file, 'UniformOutput', false);
                return;
            end
            file = char(file);

            [pathRoot, f, e] = fileparts(file);
            if isempty(e)
                error('No file extension specified on Imec file name');
            end
            file = [f, e];


%             match = regexp(file, '(?<stem>\w+).imec0.(?<type>\w+).bin', 'names', 'once');
%             match = regexp(file, '(?<stem>\w+).imec0.bin', 'names', 'once');
            if ~isempty(match)
                type = match.type;
                fileStem = match.stem;
                return;
            end

            fileStem = f;
            type = '';
end