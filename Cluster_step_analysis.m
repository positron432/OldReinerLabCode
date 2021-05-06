%MAJOR DIFFERENCES BETWEEN YOUR VERSION AND MINE
%You no longer have to set jmax. The code will automatically figure out the
%correct number for this (you will get an error in the process but ignore
%it, the code will continue to run). Delta has been renamed exptStepSize,
%and a new variable called hmulti is now defined on line 23, which replaces
%the number that was changed before for each run. Edit it like you would
%the number from before.


clearvars -except exptNumStep exptStepSize hmulti cur;


fileID = fopen('C:\Users\PHYS_reiner\Documents\MBA data analysis\ANALYSIS\4-MBA\H2O\pump_off_capture\6_11_08\2019_06_11_0008_filtered_100Hz.txt');
cur = fscanf(fileID,'%f');
w = size(cur);
length = w(1);

locationswitch = 100;
locationup(1) = 1;
locationdown(1) = 1;

%This is the hmulti that you edit now to change the sensitivity
hmulti = 600;
jmax = 10000;
exptStepSize = 2.3;
for(m=1:2)
try
for(j=1:jmax)
    
    clear Sdown;
    clear Sup;
    start = 1;

if(locationswitch == 1)
    start = locationup(j)+1; 
end
if(locationswitch == 0)
    start = locationdown(j)+1;
end



subcur = cur(start:start+1500);

mu(j) = mean(subcur);
sig = std(subcur);


delta = exptStepSize;

A = delta/sig^2;
%h = 700*delta^2/sig^2;    %This works for longer steps
h = hmulti*delta/sig;


%if(h<1500)
%    h=2000;
%end

Gup(start) = 0;
Sup(start) = 0;
Gdown(start) = 0;
Sdown(start) = 0;
i=start;

%Steps locator
while(Gup(i) <= h && Gdown(i) >= -1*h && i<length)
    smallsup(i+1) = A*(cur(i) - (mu(j)+delta/2)); % running current for upsteps
    Sup(i+1) = Sup(i) + smallsup(i+1);
    Gup(i+1) = max(0,Gup(i) + smallsup(i+1));
    if(Gup(i+1) > h)
        locationup(j+1) = find(Sup==min(Sup));
        locationswitch = 1;
        
    end
    
    smallsdown(i+1) = A*(cur(i) - (mu(j)-delta/2));  % running current for downsteps
    Sdown(i+1) = Sdown(i) + smallsdown(i+1);
    Gdown(i+1) = min(0,Gdown(i) + smallsdown(i+1));
    
    if(Gdown(i+1) < -1*h)
        locationdown(j+1) = find(Sdown==max(Sdown));
        locationswitch = 0;
        
    end
    
   i=i+1;
   
end


end

sizemu = size(mu);
sizemul = sizemu(2);

sizelocdown = size(locationdown);
sizelocup = size(locationup);

sizelocdownlength = sizelocdown(2);
sizelocuplength = sizelocup(2);


if(sizelocdownlength > sizelocuplength)
    locationup(sizelocdownlength) = 0;
else
    locationdown(sizelocuplength) = 0;
end

location = locationup+locationdown;
location(1) = [];


for(k=1:2*(sizemul-1))
    test = mod(k,2);
    
    if(test == 1)
    locationdouble(k) = location((k+1)/2);    
    end
    
    if(test == 0)
    locationdouble(k) = location(k/2);
    end
    
    
end




munew(1) = mean(cur(1:location(1)));

for(k=1:sizemul-1)
    
lower = location(k);
upper = location(k+1);

munew(k+1) = mean(cur(lower:upper));
    
end



for(k=1:2*sizemul)
    test = mod(k,2);
    
    if(test == 1)
    mudouble(k) = munew((k+1)/2);    
    end
    
    if(test == 0)
    mudouble(k) = munew(k/2);
    end
    
    
end




locationdouble = [0 locationdouble(1:end)];
mudouble(end) =[]; 



for(k=1:jmax - 1)
locationdiff(k) = location(k+1) - location (k);
mudiff(k) = munew(k+1) - munew(k);
end

clear jmax;

locationdiff = [0 locationdiff(1:end)];
mudiff = [0 mudiff(1:end)];

mut = transpose(munew);
locationt = 20e-6*transpose(location);
mudifft = transpose(mudiff);
locationdifft = 20e-6*transpose(locationdiff);


out = [mut locationt mudifft locationdifft];

dlmwrite('C:\Users\PHYS_reiner\Documents\MBA data analysis\ANALYSIS\4-MBA\H2O\pump_off_capture\CUMSUM_2019_06_11_0008_filtered_100Hz.txt',out);

%dlmwrite('C:\Users\jereiner\Documents\My Lab\Data\2018\May\May30\CUMESUMresults_mu_file31.txt',mut);
%dlmwrite('C:\Users\jereiner\Documents\My Lab\Data\2018\May\May30\CUMESUMresults_location_file31.txt',locationdoublet);

catch
jmax = j - 1;
end
end
t = linspace(0,length,length);
figure()

plot(t,cur,'r',locationdouble,mudouble,'b')
