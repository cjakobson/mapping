%% Import VCF to get genotype data

%Ex.
%  filename='/Users/richardshe/Documents/F6cross_VCFs/Plate1.output.vcf';
%  af_upper = 0.6; af_lower = 0.4; ac_ratio = 90;
function [vcf,header,af,filler] = read_vcf(filename)

fid = fopen(filename);
filler = [];
temp = 'start';
while strcmp(temp(1:2),'#C') == 0
    temp = fgets(fid);
    if strcmp(temp(1:2),'#C') == 0
        filler = [filler,temp];
    end
    if isempty(temp)
        temp = 'empty';
    end
end

header = strsplit(temp,char(9));
vcf = textscan(fid,repmat('%s',1,length(header)),'delimiter',char(9));

fclose(fid);

info = vcf{8};
for i = 1:length(info)
    left = strfind(info{i},'AF=');
    right = strfind(info{i}(left(1):length(info{i})),';');
    afstring = info{i}(left(1)+3:left(1)+right(1)-2);
    if ~isempty(strfind(afstring,','));
        afstring = strsplit(afstring,',');
        if ~ischar(afstring)
            for j = 1:length(afstring)
                af{i}(j) = afstring(j);
            end
        else
            af{i} = str2num(afstring);
        end
    else
        af{i} = str2num(afstring);
    end
    left = strfind(info{i},'AC=');
    right = strfind(info{i}(left(1):length(info{i})),';');
    acstring = info{i}(left(1)+3:left(1)+right(1)-2);
    if ~isempty(strfind(acstring,','));
        acstring = strsplit(acstring,',');
        if ~ischar(acstring)
            for j = 1:length(acstring)
                ac{i}(j) = acstring(j);
            end
        else
            ac{i} = str2num(acstring);
        end
    else
        ac{i} = str2num(acstring);
    end
end

% for i = 1:length(af)
%     afmax(i) = max(af{i});
%     acmax(i) = max(ac{i});
% end




