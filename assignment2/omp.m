
function x = omp(d,k)

    C = dctmtx(d);
    I = eye(d);
    D = cat(2,I,C);
    x = randperm(d);
    x = x';
    r = x
    i = 1;
    x = zeros(d,1);
    u = zeros(k);
    while i<k
        ip=[];
        for j =  1:size(D,2)
            ip(j) = dot(r,D(:,j));
        end
        [d_product , idnx] = max(abs(ip));
        if i==1
            u(:,1) = D(:,idnx);
        else 
            negval = zeros(d,1);
            for j = 1:i-1
                negval =negval +  (dot(D(:,idnx),u(:,j)) / dot(u(:,j),u(:,j)))*u(:,j);
            u(:,i) = D(:,idnx) - negval;
            end
        end
        error_val = dot(r,u(:,i)) * u(:,i) /dot(u(:,i),u(:,i));
        x = x + error_val;
        r = r - error_val;
        dot(r,D(:,idnx));
        i = i+1;
    end
    x;
end

        
        