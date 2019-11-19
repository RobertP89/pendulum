function mytf_ip = imperfectTF( mytf, fcut)
mytf_ip = mytf;
[num_a, den_a] = tfdata(mytf, 'v');
order_num_a = 0;
order_den_a = 0;

% If length of num is 3 then max order is 2
% but the actual order is determined from the lowest index nonzero entry

for i = 1:length(num_a)
    if num_a(i) ~= 0
        order_num_a = length(num_a)-i;
        break
    end
end

for i = 1:length(den_a)
    if den_a(i) ~= 0
        order_den_a = length(den_a)-i;
        break
    end
end

pzDiff_a = order_num_a - order_den_a;


fcut_a = fcut;
wcut_a = 2*pi*fcut_a;
LPF_a = tf(wcut_a, [1 wcut_a]);

if pzDiff_a > 0
    for i = 1:pzDiff_a
        mytf_ip = LPF_a * mytf_ip;
    end
end
end