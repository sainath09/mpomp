
function x = mp(d,k)

    C = dctmtx(d);
    I = eye(d);
    D = cat(2,I,C);
    x = randperm(d);
    x = x';
    r = x
    i = 0;
    x = zeros(d,1);
    while i<k
        ip=[];
        for j =  1:size(D,2)
            ip(j) = dot(r,D(:,j));
        end
        [d_product , idnx] = max(abs(ip));
        x = x + ip(idnx)*D(:,idnx);
        r = r - ip(idnx)*D(:,idnx);
        dot(r,D(:,idnx));
        i = i+1;
    end
    x;
end

        
        