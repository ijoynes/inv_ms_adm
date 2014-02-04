close
clear
clc
node_index = [9789,27923];

nt = 901;
U_mag = nan(nt,2);
U_dir = nan(nt,2);
for i = 1 : nt
    i
    load(['E:\ijoynes\thesis_data_backup\run_023\Flow\Flow_' num2str(i-1) '.mat'],'uv');
    U_mag(i,:) = sqrt(sum(uv(node_index,:).^2,2))';
    U_dir(i,:) = 270-atan2(uv(node_index,2),uv(node_index,1))'*180/pi;
end

plot(U_dir(:,2))
for i = 1 : nt
    for j = 1 : 2
        while U_dir(i,j) >= 360
            U_dir(i,j) = U_dir(i,j) - 360;
        end
        while U_dir(i,j) < 0
            U_dir(i,j) = U_dir(i,j) + 360;
        end
    end
end
plot(U_dir(:,2))

for i = 2 : nt
    for j = 1 : 2
        if U_dir(i,j) > U_dir(i-1,j) + 180
            U_dir(i,j) = U_dir(i,j) - 360;
        elseif U_dir(i,j) < U_dir(i-1,j) - 180
            U_dir(i,j) = U_dir(i,j) + 360;
        end
            
            
    end
end
plot(U_dir(:,1))


plot(
% plot(U_dir(:,2))
% for i = 2 : nt-1
%     for j = 1 : 2
%         if abs(U_dir(i,j)-U_dir(i-1,j)) > 180
%             U_dir(i,j) = U_dir(i,j) + 360;
%         end
%         
%     end
% end
% for i = 2 : nt-1
%     for j = 1 : 2
%         if abs(U_dir(i,j)-U_dir(i-1,j)) > 180
%             U_dir(i,j) = U_dir(i,j) - 360;
%         end
%         
%     end
% end
plot(U_dir(1:end-1,2)-U_dir(2:end,2))


t = (0:840)/60;
time_delay = [1,2,5,10,20,30,60];
for i = 1 : 1
figure(i)
%plotyy(t,U_mag(1:841,i)-U_mag(1+time_delay(1):841+time_delay(1),i),t,U_dir(1:841,i)-U_dir(1+time_delay(1):841+time_delay(1),i))
plot(t,U_dir(1:841,i)-U_dir(1+time_delay(1):841+time_delay(1),i))
axis([0,14,-90,90])
end

U_dir_rms = nan(61,1);
U_vel_rms = nan(61,1);

for i = 0 : 60
U_dir_rms(i+1) = sqrt(mean((U_dir(1:841,1)-U_dir(1+i:841+i,1)).^2));
U_vel_rms(i+1) = sqrt(mean((U_mag(1:841,1)-U_mag(1+i:841+i,1)).^2));
end

plot(0:60,U_dir_rms,0:60,U_vel_rms);
grid

