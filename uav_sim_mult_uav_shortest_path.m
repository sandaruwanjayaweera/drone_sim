close all
clear all

N_uav               = 3;
altitude            = 8;
d_in_range          = 10;
Init_coord          = [0 0];
model               = multirotor;

t_full          = 1*60;
t_req           = 0.5*60;
t_down          = 0.25*60;
t_refill        = 0.25*60;
charge_rate     = 40;
speed           = 10;

buildingsgeo    = readgeotable("university.osm",CoordinateSystemType="geographic"); %readgeotable("linnanmaa.osm",Layer="buildings");
buildingsLayer  = readgeotable("university.osm",Layer="buildings");

% Suppose buildings is a geotable from readgeotable(.osm)
meanLat = mean(buildingsgeo.Shape.Latitude, 'omitnan'); %mean(buildings.Centroid.Latitude, 'omitnan');
meanLon = mean(buildingsgeo.Shape.Longitude, 'omitnan'); % mean(buildings.Centroid.Longitude, 'omitnan');
origin = [meanLat, meanLon, 0];

Scenario = uavScenario(ReferenceLocation=origin, UpdateRate=1000, MaxNumFrames=N_uav*2+2);

% xTerrainLimits = [-300 400];
% yTerrainLimits = [-1300 500];
% color = [0.6 0.6 0.6];
% addMesh(scene,"terrain",{"gmted2010",xTerrainLimits,yTerrainLimits},color)
xBuildingLimits = [-500 300]; % longitudes horizontal axis
yBuildingLimits = [-400 400]; %[-1300 500]; % lattitudes vertical axis
color = [0.6431 0.8706 0.6275];
addMesh(Scenario,"buildings",{"university.osm",xBuildingLimits,yBuildingLimits,"auto"},color)

[geolimlat, geolimlon, geolimalt] = local2latlon(xBuildingLimits,yBuildingLimits,[0 0],origin);
fig1 = figure;
fig1.Position = [0 0 abs(xBuildingLimits(2) - xBuildingLimits(1)) abs(yBuildingLimits(2) - yBuildingLimits(1))];
geolimits([geolimlat(1) geolimlat(2)], [geolimlon(1) geolimlon(2)])
geoplot(buildingsLayer,ColorVariable="MaxHeight", FaceColor="black")
% print(gcf,'building_height.png','-dpng','-r600');

% drawnow
% % truesize(gcf, [abs(yBuildingLimits(2) - yBuildingLimits(1)) abs(xBuildingLimits(2) - xBuildingLimits(1))])
% % save the image using getframe() and imwrite()
% screenshot = frame2im(getframe(gca));
% imwrite(screenshot,'building_height.png')

strScreenshot = getframe(gca); % Get screenshot
rgbImage = strScreenshot.cdata;
figure
imshow(rgbImage)
axis off;
impixelinfo;
imwrite(rgbImage, 'building_height.png');

image = imread('building_height.png');
testpoint1   = [0 0];
testpoint2  = [161.174 162.851];
testpoint3  = [20 10]; % the simulation considers y axis as 1st coordinate. Consider this when creating waypoints
% RGB = insertMarker(image,[0.9724*testpoint1(1)+488.2744 -0.981*testpoint1(2)+396.7568], "x", MarkerColor="red", Size=10); % interpolation to match with pixels
% RGB = insertMarker(RGB,[0.9724*testpoint2(1)+488.2744 -0.981*testpoint2(2)+396.7568], "x", MarkerColor="green", Size=10);
% RGB = insertMarker(RGB,[0.9724*testpoint3(1)+488.2744 -0.981*testpoint3(2)+396.7568], "x", MarkerColor="blue", Size=10);
% figure
% imshow(RGB)

grayimage = rgb2gray(image);
bwimage = grayimage < 50;
% grayimage(round)
grid = binaryOccupancyMap(bwimage);
figure
show(grid)

% inflate(grid, 4) % inflate to avoid being too close to buildings

% create grid/graph to find the shortest path
xlim1   = xBuildingLimits(1);
xlim2   = xBuildingLimits(2);
ylim1   = yBuildingLimits(1);
ylim2   = yBuildingLimits(2);
% [I,J]   = ndgrid(1:(yBuildingLimits(2) - yBuildingLimits(1))/10 + 1, 1:(xBuildingLimits(2) - xBuildingLimits(1))/10 + 1);
[I,J]   = ndgrid(xlim1:speed:xlim2, ylim1:speed:ylim2);
IJ      = [I(:),J(:)];
D       = pdist2(IJ,IJ);
A1      = (D>=0.001 & D<=sqrt(200)*1.001).*D; %adjacency matrix
A2      = A1;

% create a list of unocccupied coordinates
freepoints = [];
coord_idx  = 0;
for y = ylim1:speed:ylim2
    for x = xlim1:speed:xlim2
        coord_idx = coord_idx + 1;
        if getOccupancy(grid, getOccMapCoordinates(x, y)) % if occupied
            A1(coord_idx, :) = 0;
            A1(:, coord_idx) = 0;
        else
            freepoints = [freepoints; x, y]; % Store the unoccupied coordinates
        end
        % if ~getOccupancy(grid, getOccMapCoordinates(x, y))
            % if sqrt(x^2+y^2) < 20 % for now, limit destination to nearby
                % freepoints = [freepoints; x, y]; % Store the unoccupied coordinates
                % setOccupancy(grid, getOccMapCoordinates(x, y), 1);
            % end
        % end
    end
end

G1       = graph(A1);
G2       = graph(A2);
obj     = spatialgraph2D(G1,IJ(:,1),IJ(:,2));
figure
p = plot(obj);
% set(gca, 'YDir','reverse')
% set(gca, 'XDir','reverse')

reachpoints = [];
pos_init = sub2ind([(xlim2 - xlim1)/speed + 1, (ylim2 - ylim1)/speed + 1], ...
                                    (Init_coord(1) - xlim1)/speed + 1, ...
                                    (Init_coord(2) - ylim1)/speed + 1);

for i = 1:length(freepoints)
    des_point   = freepoints(i,:);
    pos_des     = sub2ind([(xlim2 - xlim1)/speed + 1, (ylim2 - ylim1)/speed + 1], ...
                                    (des_point(1) - xlim1)/speed + 1, ...
                                    (des_point(2) - ylim1)/speed + 1);
    dis_init = distances(G1, pos_init, pos_des);
    if ~isnan(dis_init) && ~isinf(dis_init) && dis_init ~= 0
        reachpoints = [reachpoints; des_point];
    end
end

reachpoint_index = randperm(length(reachpoints));
% getOccupancy(grid, getOccMapCoordinates(testpoint1(1), testpoint1(2)))
% getOccupancy(grid, getOccMapCoordinates(testpoint2(1), testpoint2(2)))
% getOccupancy(grid, getOccMapCoordinates(testpoint3(1), testpoint3(2)))
figure
show(grid)

disp('Section 1')
%%

InitialPosition     = zeros(N_uav, 3);
InitialOrientation  = zeros(N_uav, 3);
destination         = zeros(N_uav, 3);
uavPosition         = zeros(3, N_uav);
Waypoints           = cell(N_uav,1);
plat                = cell(N_uav,1);

skip_list           = zeros(N_uav, 1);
for uav_id = 1:N_uav
    if uav_id > length(reachpoints)
        disp(string("too many UAVs for the environment. Only " + length(reachpoints) + " are considered for simulation."))
        break
    end
    InitialPosition(uav_id,:)     = [Init_coord -altitude];
    InitialOrientation(uav_id,:)  = [pi/2 0 0];
    while(true)
        destination(uav_id,:)         = [reachpoints(reachpoint_index(uav_id), :) -(altitude+0)];
    
        pos1                          = sub2ind([(xlim2 - xlim1)/speed + 1, (ylim2 - ylim1)/speed + 1], ...
                                        (InitialPosition(uav_id,1) - xlim1)/speed + 1, ...
                                        (InitialPosition(uav_id,2) - ylim1)/speed + 1);
        
        pos2                          = sub2ind([(xlim2 - xlim1)/speed + 1, (ylim2 - ylim1)/speed + 1], ...
                                        (destination(uav_id,1) - xlim1)/speed + 1, ...
                                        (destination(uav_id,2) - ylim1)/speed + 1);
    
        [posidx, len_path]            = shortestpath(G1, pos1, pos2);
        if len_path == 0 || len_path == inf
            % [posidx, len_path]            = shortestpath(G2, pos1, pos2);
            % path_points                   = [I(posidx)' J(posidx)'];
            % Waypoints{uav_id}             = [path_points zeros(length(posidx),1)-(altitude+10)];        
            disp(string("No possible path for UAV " + uav_id))
            skip_list(uav_id)         = 1;
            % continue
        else
            path_points                   = [I(posidx)' J(posidx)'];
            Waypoints{uav_id}             = [path_points zeros(length(posidx),1)-altitude];
            break
        end
    end

    uavPosition(:,uav_id)         = InitialPosition(uav_id,:)';
    plat{uav_id}                  = uavPlatform("UAV"+uav_id, Scenario, ...
                                                "ReferenceFrame","NED", ...
                                                "InitialPosition",InitialPosition(uav_id,:), ...
                                                "InitialOrientation",eul2quat(InitialOrientation(uav_id,:)));
    updateMesh(plat{uav_id},"quadrotor",{1.2},[0 0 1],eul2tform([0 0 pi]));
end  

figure
[ax,plotFrames] = show3D(Scenario);
xlim(xBuildingLimits);
ylim(yBuildingLimits);
% ax.Clipping = "off";
ax.ClippingStyle = "rectangle";
axis(ax,"equal")
hold(ax,"on")

setup(Scenario)

disp('Section 2')
%% Channel modeling (https://se.mathworks.com/help/comm/ug/mobility-modeling-with-ray-tracing-channel.html)

sv = siteviewer( ...
    Buildings="university.osm", ...
    Basemap="satellite");

fc = 5.9e9;
c = physconst("LightSpeed");
lambda = c/fc;

ueElement = phased.IsotropicAntennaElement;
ueArray = phased.NRRectangularPanelArray( ...
    Size=[1 1 1 1], ...
    ElementSet={ueElement}, ...
    Spacing=lambda*[.5 .5 1 1]);

% gNBElement = phased.NRAntennaElement( ...
%     PolarizationAngle=-45, ...
%     Beamwidth=120);
gNBElement = phased.IsotropicAntennaElement;

bsArray = phased.NRRectangularPanelArray( ...
    Size=[2 2 1 1], ...
    ElementSet={gNBElement}, ...
    Spacing=lambda*[.5 .5 1 1]);

% rxs = rxsite(Name="NR5G", ...
%     Latitude=65.059560, ...
%     Longitude=25.470443, ...
%     Antenna=bsArray, ...
%     AntennaHeight=17);
[bssite_lat, bssite_lon, bssite_height] = local2latlon(Init_coord(1), Init_coord(2), 17, origin);
rxstation = rxsite(Name="NR5G", ...
    Latitude=bssite_lat, ...
    Longitude=bssite_lon, ...
    Antenna=bsArray, ...
    AntennaHeight=bssite_height);
show(rxstation);

txstation = txsite(Latitude=bssite_lat, ...
    Longitude=bssite_lon, ...
    Antenna=ueArray, ...
    TransmitterPower=1, ...
    TransmitterFrequency=fc, ...
    AntennaHeight=bssite_height);
show(txstation);

numreflecs = 2;
numdiffracs = 1;
maxrelpathloss = 20;
rtpropmdl = propagationModel("raytracing", ...
    MaxNumReflections=numreflecs, ...
    MaxNumDiffractions=numdiffracs, ...
    MaxRelativePathLoss=maxrelpathloss);

longricepropmdl = propagationModel("longley-rice");


disp('Section 3')
%%
% Simulate 
rrays           = cell(N_uav,1);
rlos            = cell(N_uav,1);
rssrt           = cell(N_uav,1);
rsslongrice     = cell(N_uav,1);
goal_reach      = zeros(N_uav, 1);
uav_trav_dist   = zeros(N_uav, 1);
rss_idx         = zeros(N_uav, 1);
waypoint_idx    = zeros(N_uav, 1);
uav_speed       = ones(1, N_uav)*speed;

uav_denm        = zeros(N_uav, 3); % fuel amount, req_id, arrival_time
uav_denm(:,1)   = t_full*speed;
uav_denm(:,2)   = inf;
uav_denm(:,3)   = 0;
req_id          = 1;
rss_idx_rep     = zeros(N_uav, 1);
rrays_rep       = cell(N_uav,1);
rlos_rep        = cell(N_uav,1);
rssrt_rep       = cell(N_uav,1);
rsslongrice_rep = cell(N_uav,1);

uav_state       = zeros(N_uav,1); % 0 = idle, 1 = req_sent, 2 = accepted, 3 = rejected, 4 = moving to BS
bs_state        = 0; % 0 = unreserved, 1 = reserved, 2 = drone_charging
bs_req_count    = 0;

charge_counter      = zeros(N_uav, 1);
dis_list            = zeros(N_uav, 1) + t_full*speed;

reject_counter      = zeros(N_uav, 1);
speed_idx           = zeros(N_uav, 1);
prev_destination    = zeros(N_uav, 3);

cm = parula(N_uav);

idx = 0;
req_count       = 5; % the most probable DENM message count during simulation
while idx < (t_full-t_req)*req_count+1          
    idx             = idx + 1;
    waypoint_idx    = waypoint_idx + 1;
    speed_idx       = speed_idx + 1;
    for uav_id = 1:N_uav
        
        if(uav_state(uav_id) == 2)
            uav_state(uav_id)               = 4 % moving to BS

            prev_destination(uav_id,:)      = destination(uav_id,:);
            InitialPosition(uav_id,:)       = uavPosition(:,uav_id)';
            destination(uav_id,:)           = [Init_coord -altitude];
            pos1                            = sub2ind([(xlim2 - xlim1)/speed + 1, (ylim2 - ylim1)/speed + 1], ...
                                                (InitialPosition(uav_id,1) - xlim1)/speed + 1, ...
                                                (InitialPosition(uav_id,2) - ylim1)/speed + 1);
            
            pos2                            = sub2ind([(xlim2 - xlim1)/speed + 1, (ylim2 - ylim1)/speed + 1], ...
                                                (destination(uav_id,1) - xlim1)/speed + 1, ...
                                                (destination(uav_id,2) - ylim1)/speed + 1);

            [posidx, len_path]              = shortestpath(G1, pos1, pos2);
            if len_path == 0 || len_path == inf
                % [posidx, len_path]        = shortestpath(G2, pos1, pos2);
                % path_points               = [I(posidx)' J(posidx)'];
                % Waypoints{uav_id}         = [path_points zeros(length(posidx),1)-(altitude+10)];        
                disp(string("No possible path to reach BS " + uav_id))
                skip_list(uav_id)           = 1;
            else
                path_points                 = [I(posidx)' J(posidx)'];
                Waypoints{uav_id}           = [path_points zeros(length(posidx),1)-altitude];
            end
            waypoint_idx(uav_id)            = 1;
        end

        if(uav_state(uav_id) == 3) % rejected
            reject_counter(uav_id)          = reject_counter(uav_id) + 1
            if(reject_counter(uav_id) >= t_refill) % wait for refill to finish and make another request
                reject_counter(uav_id)      = 0;
                uav_state(uav_id)           = 0
            end
        end

        proximity = destination(uav_id,:)' - uavPosition(:,uav_id);
        if (sqrt((proximity'*proximity)) <= d_in_range) || (skip_list(uav_id) == 1)
            goal_reach(uav_id) = 1;
            % continue;

            if(uav_state(uav_id) == 4) % moved to BS
                charge_counter(uav_id)      = charge_counter(uav_id) + 1
                bs_state                    = 2;
                % if(charge_counter(uav_id) >= t_refill) % fully charged
                if(charge_counter(uav_id) >= (t_full*speed - dis_list(uav_id))/charge_rate) % fully charged
                    charge_counter(uav_id)  = 0;
                    uav_state(uav_id)       = 0
                    bs_state                = 0;
                    uav_trav_dist(uav_id)   = 0;
                    dis_list(uav_id)       = t_full*speed;
                    speed_idx(uav_id)       = 2;
                    % waypoint_idx(uav_id)    = length(Waypoints{uav_id, 1});

                    InitialPosition(uav_id,:)       = uavPosition(:,uav_id)';
                    destination(uav_id,:)           = prev_destination(uav_id,:);
                    pos1                            = sub2ind([(xlim2 - xlim1)/speed + 1, (ylim2 - ylim1)/speed + 1], ...
                                                        (InitialPosition(uav_id,1) - xlim1)/speed + 1, ...
                                                        (InitialPosition(uav_id,2) - ylim1)/speed + 1);
                    
                    pos2                            = sub2ind([(xlim2 - xlim1)/speed + 1, (ylim2 - ylim1)/speed + 1], ...
                                                        (destination(uav_id,1) - xlim1)/speed + 1, ...
                                                        (destination(uav_id,2) - ylim1)/speed + 1);
        
                    [posidx, len_path]              = shortestpath(G1, pos1, pos2);
                    if len_path == 0 || len_path == inf
                        % [posidx, len_path]        = shortestpath(G2, pos1, pos2);
                        % path_points               = [I(posidx)' J(posidx)'];
                        % Waypoints{uav_id}         = [path_points zeros(length(posidx),1)-(altitude+10)];        
                        disp(string("No possible path to reach prev_destination " + uav_id))
                        skip_list(uav_id)           = 1;
                    else
                        path_points                 = [I(posidx)' J(posidx)'];
                        Waypoints{uav_id}           = [path_points zeros(length(posidx),1)-altitude];
                    end
                    waypoint_idx(uav_id)            = 2;
                else
                    waypoint_idx(uav_id)    = length(Waypoints{uav_id, 1}) - 1;
                    % waypoint_idx(uav_id)    = length(Waypoints{uav_id, 1});
                end
            else
                % When destination is reached, update destination
                InitialPosition(uav_id,:)     = uavPosition(:,uav_id)';
                while(true)
                    destination(uav_id,:)         = [reachpoints(randi(length(reachpoints)), :) -(altitude+0)];
                    pos1                          = sub2ind([(xlim2 - xlim1)/speed + 1, (ylim2 - ylim1)/speed + 1], ...
                                                    (InitialPosition(uav_id,1) - xlim1)/speed + 1, ...
                                                    (InitialPosition(uav_id,2) - ylim1)/speed + 1);
                    
                    pos2                          = sub2ind([(xlim2 - xlim1)/speed + 1, (ylim2 - ylim1)/speed + 1], ...
                                                    (destination(uav_id,1) - xlim1)/speed + 1, ...
                                                    (destination(uav_id,2) - ylim1)/speed + 1);
                
                    [posidx, len_path]            = shortestpath(G1, pos1, pos2);
                    if len_path == 0 || len_path == inf
                        % [posidx, len_path]            = shortestpath(G2, pos1, pos2);
                        % path_points                   = [I(posidx)' J(posidx)'];
                        % Waypoints{uav_id}             = [path_points zeros(length(posidx),1)-(altitude+10)];        
                        disp(string("No possible path for UAV " + uav_id))
                        skip_list(uav_id)         = 1;
                    else
                        path_points                   = [I(posidx)' J(posidx)'];
                        Waypoints{uav_id}             = [path_points zeros(length(posidx),1)-altitude];
                        break
                    end
                end
                waypoint_idx(uav_id)              = 1;
            end
        end
    
        if(waypoint_idx(uav_id)==1)
            uav_trav_dist(uav_id) = uav_trav_dist(uav_id) + norm(Waypoints{uav_id, 1}(waypoint_idx(uav_id),:) - InitialPosition(uav_id,:)); % dist(Waypoints{uav_id, 1}(idx,:), InitialPosition(uav_id,:));
            dis_list(uav_id)      = dis_list(uav_id) - (norm(Waypoints{uav_id, 1}(waypoint_idx(uav_id),:) - InitialPosition(uav_id,:)));
        else
            uav_trav_dist(uav_id) = uav_trav_dist(uav_id) + norm(Waypoints{uav_id, 1}(waypoint_idx(uav_id),:) - Waypoints{uav_id, 1}(waypoint_idx(uav_id)-1,:));
            dis_list(uav_id)      = dis_list(uav_id) - (norm(Waypoints{uav_id, 1}(waypoint_idx(uav_id),:) - Waypoints{uav_id, 1}(waypoint_idx(uav_id)-1,:)));
        end

        % Update states.
        uavPosition(:,uav_id)       = Waypoints{uav_id, 1}(waypoint_idx(uav_id),:)';
      
        % Plot the path.
        desPosition = plot3(uavPosition(1, uav_id), uavPosition(2, uav_id), -uavPosition(3, uav_id), ".", 'Color', cm(uav_id,:), Parent=ax);
        goalPosition = plot3(destination(uav_id, 1), destination(uav_id, 2), -destination(uav_id, 3), "ob", Parent=ax);
        legend(ax,[desPosition goalPosition],["Desired Position", "Goal"])
    
        % Advance scene simulation time.
        advance(Scenario);

        % logic to make the DENM request
        % uav_speed(uav_id) = uav_trav_dist(uav_id)/speed_idx(uav_id);
        % if(uav_trav_dist(uav_id)/uav_speed(uav_id) >= (t_full-t_req) && uav_state(uav_id) == 0) % 750/5 >= 150
        if(dis_list(uav_id) <= t_req*speed && uav_state(uav_id) == 0) % 750/5 >= 150
            uav_state(uav_id)   = 1 % req sent
            uav_trav_dist
            speed_idx
            uav_speed

            [txsite_lat, txsite_lon, txsite_height] = local2latlon(uavPosition(1,uav_id), uavPosition(2,uav_id), uavPosition(3,uav_id), origin);
        
            tx = txsite(Latitude=txsite_lat, Longitude=txsite_lon, Antenna=ueArray, ...
                        TransmitterPower=1, ...
                        TransmitterFrequency=fc, ...
                        AntennaHeight=-txsite_height);
            show(tx, ShowAntennaHeight=true)
    
            tx.AntennaAngle = angle(tx,rxstation);
            raytrace_t                     = raytrace(tx,rxstation,rtpropmdl);
            if(~isempty(raytrace_t{1,1})) % check if it is possible to send signal
                rss_idx(uav_id)     = rss_idx(uav_id) + 1;
                rrays{uav_id}(rss_idx(uav_id)) = raytrace_t; 
                rlos{uav_id}(rss_idx(uav_id))  = los(tx,rxstation);
        
                rssrt{uav_id}(rss_idx(uav_id))          = computeSigStrengthFromRays(rrays{uav_id}(rss_idx(uav_id)), tx, rxstation, fc);
                rsslongrice{uav_id}(rss_idx(uav_id))    = sigstrength(rxstation, tx, longricepropmdl);

                % UAV-DENM information
                uav_denm(uav_id, 1)    = dis_list(uav_id)
                uav_denm(uav_id, 2)    = req_id
                req_id                 = req_id + 1;
                bs_req_count           = bs_req_count + 1;
            else
                disp('Cannot send req from UAV');
            end
        end
    end

    if(bs_req_count > 0) % BS replying to UAVs
        if(bs_state == 0)
            immediate_uav                                       = find((uav_denm(:,1) == min(uav_denm(:,1))) & (uav_state == 1)) % find UAV with (minimum fuel) & (UAV which has sent request)
            immediate_uav                                       = min(immediate_uav) % find uav with first request
            [immediate_lat, immediate_lon, immediate_height]    = local2latlon(uavPosition(1,immediate_uav), uavPosition(2,immediate_uav), uavPosition(3,immediate_uav), origin);
            immediate_rx = rxsite(Name="NR5G", ...
                Latitude=immediate_lat, ...
                Longitude=immediate_lon, ...
                Antenna=bsArray, ...
                AntennaHeight=-immediate_height);
            show(immediate_rx);
            raytrace_t                                          = raytrace(txstation, immediate_rx, rtpropmdl);   
            if(~isempty(raytrace_t{1,1})) % check if it is possible to send signal
    
                % UAV-DENM information
                % uav_denm(:,1)   = t_full;
                % uav_denm(:,2)   = inf;
                uav_denm(immediate_uav,1)   = t_full*speed;
                uav_denm(immediate_uav,2)   = inf;
                % req_id                      = req_id + 1;
                if(bs_state == 0) % unreserved (accept)
                    bs_req_count                                                    = 0;
                    % bs_req_count                                                  = bs_req_count - 1;
                    rss_idx_rep(immediate_uav)                                      = rss_idx_rep(immediate_uav) + 1;
                    rrays_rep{immediate_uav}(rss_idx_rep(immediate_uav))            = raytrace_t; 
                    rlos_rep{immediate_uav}(rss_idx_rep(immediate_uav))             = los(txstation,immediate_rx);
                    rssrt_rep{immediate_uav}(rss_idx_rep(immediate_uav))            = computeSigStrengthFromRays(rrays_rep{immediate_uav}(rss_idx_rep(immediate_uav)), txstation, immediate_rx, fc);
                    rsslongrice_rep{immediate_uav}(rss_idx_rep(immediate_uav))      = sigstrength(immediate_rx, txstation, longricepropmdl);
                    uav_state(immediate_uav)                                        = 2 % accepted
    
                    % bs_req_count                                                  = bs_req_count - sum(uav_state == 1);
                    if(sum(uav_state == 1) > 0)
                        for(rej_uav_id = find(uav_state == 1)')
                            [rej_uav_lat, rej_uav_lon, rej_uav_height]                  = local2latlon(uavPosition(1,rej_uav_id), uavPosition(2,rej_uav_id), uavPosition(3,rej_uav_id), origin);
                            rej_uav_rx = rxsite(Name="NR5G", ...
                                Latitude=rej_uav_lat, ...
                                Longitude=rej_uav_lon, ...
                                Antenna=bsArray, ...
                                AntennaHeight=-rej_uav_height);
                            show(rej_uav_rx);
                            raytrace_t                                                  = raytrace(txstation, rej_uav_rx, rtpropmdl);
                            if(~isempty(raytrace_t{1,1})) % check if it is possible to send signal
                                rss_idx_rep(rej_uav_id)                                 = rss_idx_rep(rej_uav_id) + 1;
                                rrays_rep{rej_uav_id}(rss_idx_rep(rej_uav_id))          = raytrace_t; 
                                rlos_rep{rej_uav_id}(rss_idx_rep(rej_uav_id))           = los(txstation,rej_uav_rx);
                                rssrt_rep{rej_uav_id}(rss_idx_rep(rej_uav_id))          = computeSigStrengthFromRays(rrays_rep{rej_uav_id}(rss_idx_rep(rej_uav_id)), txstation, immediate_rx, fc);
                                rsslongrice_rep{rej_uav_id}(rss_idx_rep(rej_uav_id))    = sigstrength(rej_uav_rx, txstation, longricepropmdl);
                            else
                                disp('Cannot send rej from BS 1');
                            end
                        end
                    end
                    uav_state(uav_state == 1)                                       = 3 % reject everyone else
    
                    bs_state                                                        = 1;
                else % reserved (rej)
                    bs_req_count                                                    = 0;
                    % bs_req_count                                                  = bs_req_count - sum(uav_state == 1);
                    for(rej_uav_id = find(uav_state == 1)')
                        [rej_uav_lat, rej_uav_lon, rej_uav_height]                  = local2latlon(uavPosition(1,rej_uav_id), uavPosition(2,rej_uav_id), uavPosition(3,rej_uav_id), origin);
                        rej_uav_rx = rxsite(Name="NR5G", ...
                            Latitude=rej_uav_lat, ...
                            Longitude=rej_uav_lon, ...
                            Antenna=bsArray, ...
                            AntennaHeight=-rej_uav_height);
                        show(rej_uav_rx);
                        raytrace_t                                                  = raytrace(txstation, rej_uav_rx, rtpropmdl);
                        if(~isempty(raytrace_t{1,1})) % check if it is possible to send signal
                            rss_idx_rep(rej_uav_id)                                 = rss_idx_rep(rej_uav_id) + 1;
                            rrays_rep{rej_uav_id}(rss_idx_rep(rej_uav_id))          = raytrace_t; 
                            rlos_rep{rej_uav_id}(rss_idx_rep(rej_uav_id))           = los(txstation,rej_uav_rx);
                            rssrt_rep{rej_uav_id}(rss_idx_rep(rej_uav_id))          = computeSigStrengthFromRays(rrays_rep{rej_uav_id}(rss_idx_rep(rej_uav_id)), txstation, immediate_rx, fc);
                            rsslongrice_rep{rej_uav_id}(rss_idx_rep(rej_uav_id))    = sigstrength(rej_uav_rx, txstation, longricepropmdl);
                        else
                            disp('Cannot send rej from BS 2');
                        end
                    end
                    uav_state(uav_state == 1)                                       = 3 % reject everyone
                end
            else
                disp('Cannot send accept from BS');
            end
        else
            bs_req_count                                                    = 0;
            % bs_req_count                                                  = bs_req_count - sum(uav_state == 1);
            for(rej_uav_id = find(uav_state == 1)')
                [rej_uav_lat, rej_uav_lon, rej_uav_height]                  = local2latlon(uavPosition(1,rej_uav_id), uavPosition(2,rej_uav_id), uavPosition(3,rej_uav_id), origin);
                rej_uav_rx = rxsite(Name="NR5G", ...
                    Latitude=rej_uav_lat, ...
                    Longitude=rej_uav_lon, ...
                    Antenna=bsArray, ...
                    AntennaHeight=-rej_uav_height);
                show(rej_uav_rx);
                raytrace_t                                                  = raytrace(txstation, rej_uav_rx, rtpropmdl);
                if(~isempty(raytrace_t{1,1})) % check if it is possible to send signal
                    rss_idx_rep(rej_uav_id)                                 = rss_idx_rep(rej_uav_id) + 1;
                    rrays_rep{rej_uav_id}(rss_idx_rep(rej_uav_id))          = raytrace_t; 
                    rlos_rep{rej_uav_id}(rss_idx_rep(rej_uav_id))           = los(txstation,rej_uav_rx);
                    rssrt_rep{rej_uav_id}(rss_idx_rep(rej_uav_id))          = computeSigStrengthFromRays(rrays_rep{rej_uav_id}(rss_idx_rep(rej_uav_id)), txstation, immediate_rx, fc);
                    rsslongrice_rep{rej_uav_id}(rss_idx_rep(rej_uav_id))    = sigstrength(rej_uav_rx, txstation, longricepropmdl);
                else
                    disp('Cannot send rej from BS 3');
                end
            end
            uav_state(uav_state == 1)                                       = 3 % reject everyone            
        end
    end

    % if sum(goal_reach) == length(goal_reach)
    %     break
    % end
end

disp('Section 4')
%%
if sum(goal_reach) >= sum(skip_list)
    figure
    hold on
    for uav_id = 1:N_uav
        if (skip_list(uav_id) == 1)
            continue
        else
            % use RSS and receiver sensitivities to find the PER here
            plot([1:length(rsslongrice{uav_id})], rsslongrice{uav_id}, ':', 'Color', cm(uav_id,:), 'DisplayName', string("RSS longrice uav " +uav_id));
            plot([1:length(rssrt{uav_id})], rssrt{uav_id}, '-', 'Color', cm(uav_id,:), 'DisplayName', string("RSS RT uav " +uav_id));
        end
    end
    hold off
    legend()
    % xlim([0 360]); ylim([-80 -20])
    title("Signal Strength Along Route")
else
    disp('No reachable destination was generated.')
end
disp('Section 5')

function coord = getOccMapCoordinates(show3dX, show3Dy)
    coord = [0.979*show3dX + 487.2107 0.981*show3Dy + 419.2432];
    % coord = [1.2899*show3dX + 1125.1004 1.2971*show3Dy + 514.7704];
end

function ss = computeSigStrengthFromRays(rays,tx,rx,fc)
N11 = length(rays{1});
linRxPWRtmp = 0;
dbP1Tpwr = 10 * log10(tx.TransmitterPower);
azTxToRx = [];
elTxToRx = [];
azRxToTx = [];
elRxToTx = [];
for i = 1:N11
    azTxToRx(i) = ...
        rem(rays{1}(i).AngleOfDeparture(1) - tx.AntennaAngle,180);
    elTxToRx(i) = rays{1}(i).AngleOfDeparture(2);
    azRxToTx(i) = ...
        rem(rays{1}(i).AngleOfArrival(1) - rx.AntennaAngle,180);
    elRxToTx(i) = rays{1}(i).AngleOfArrival(2);
end
txPatternGain = directivity(tx.Antenna,fc,[azTxToRx; elTxToRx]);
rxPatternGain = directivity(rx.Antenna,fc,[azRxToTx; elRxToTx]);
for i = 1:N11
    % +30 to convert to dBm
    dbRxPWR_i = dbP1Tpwr + 30 + txPatternGain(i) + ...
        rxPatternGain(i) - rays{1}(i).PathLoss;
    linRxPWRtmp = linRxPWRtmp + ...
        exp(1i*rays{1}(i).PhaseShift) * 10^(dbRxPWR_i/20);
end
estTotPwr11 = 20 * log10(abs(linRxPWRtmp));
ss = estTotPwr11;
end