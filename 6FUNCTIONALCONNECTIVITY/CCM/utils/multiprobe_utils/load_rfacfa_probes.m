function [rfa_probe, cfa_probe, analogin] = load_rfacfa_probes(path)

analogin = load(strcat(path, '\analogin.mat')).analogin;

if isfile(strcat(path, '\neurons_probe3_curated'))
    cfa_probe = load(strcat(path, '\neurons_probe3_curated')).neurons;
else
    cfa_probe = load(strcat(path, '\neurons_probe3')).neurons;
end

if isfile(strcat(path, '\neurons_probe4_curated'))
    rfa_probe = load(strcat(path, '\neurons_probe4_curated')).neurons;
else
    rfa_probe = load(strcat(path, '\neurons_probe4')).neurons;
end
