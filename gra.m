%%%% Data Uplaod %%%%
dist_data=readtable("dist_cepii.xls");
GDP_data=readtable("GDP.csv");
Exp_data=readtable("Export.csv");

%%%% Data Cleaning %%%%
dist_data.Properties.VariableNames
dist_data=removevars(dist_data,["col45","contig","comlang_off","contig","comlang_ethno" ...
    ,"colony","comcol","curcol","col45","smctry","curcol","distw","distwces","distcap"]);
unique(Exp_data.LOCATION)
%ctry =Exp_data('AUS','AUT','BEL','BRA','CAN','CHE','CHL')
%Exp_data=Exp_data(Exp_data.LOCATION('ctry'),:)
%Exp_data=Exp_data.{('AUS','AUT','BEL','BRA','CAN','CHE','CHL',),:}
%Exp_data=Exp_data.LOCATION([1,2,3,4,5,6,7,],:)
%a = magic(5);
%a(1,3);
Selected_countries = {'AUS','AUT','BEL','BRA','CAN','CHE','CHL','CHN','COL','CRI','CZE','DEU','DNK'};
Exp_data=Exp_data(ismember(Exp_data.LOCATION,Selected_countries),:);
Exp_data.Properties.VariableNames;
%Exp_data1=removevars(Exp_data,["Time","UnitCode","PowerCodeCode","PowerCode","ReferencePeriodCode" ...
%   ,"ReferencePeriod","FlagCodes","Flags","Unit"]);
%joinedData1 = innerjoin(Exp_data1,GDP_data,"Keys","LOCATION");
GDP_data.Properties.VariableNames
GDP_data1=removevars(GDP_data,["Country","TRANSACT","Transaction","MEASURE","Measure" ...
    ,"TIME","UnitCode","Unit","PowerCodeCode","PowerCode","ReferencePeriodCode","ReferencePeriod","FlagCodes", ...,
    "Flags"]);
%joinedData1 = innerjoin(Exp_data1,GDP_data1,"Keys","LOCATION");
%Final_Data2 = join(dist_data,joinedData1,"LeftKeys","iso_o","RightKeys",...
%    "LOCATION")

% Join tables
%joinedData3 = outerjoin(joinedData1,dist_data,"Type","left","LeftKeys",...
%   "LOCATION","RightKeys","iso_o",true,)
%joinedData = outerjoin(GDP_data1,Exp_data1,"Keys","LOCATION");
%Final_Data2 = outerjoin(dist_data,joinedData,"LeftKeys","iso_o","RightKeys",...
%   "LOCATION_GDP_data1");
%joinedData1 = innerjoin(GDP_data1,Exp_data1,"Keys","LOCATION")
%joinedData1 = join(GDP_data1,Exp_data1,"Keys","LOCATION")
%Selected_countries = {'AUS','AUT','BEL','BRA','CAN','CHE','CHL','CHN','COL','CRI','CZE','DEU','DNK'};
%Exp_data1=Exp_data(ismember(Exp_data.LOCATION,Selected_countries),:);
%joinedData1 = innerjoin(Exp_data1,GDP_data1,"Keys","LOCATION");
%Final_Data2 = innerjoin(dist_data,joinedData1,"LeftKeys","iso_o","RightKeys",...
%   "LOCATION");

%%% Data Merging %%%
Selected_countries = {'AUS','AUT','BEL','BRA','CAN','CHE','CHL','CHN','COL','CRI','CZE','DEU','DNK'};
Exp_data=Exp_data(ismember(Exp_data.LOCATION,Selected_countries),:);
Exp_data.Properties.VariableNames;
Exp_data1=Exp_data1(ismember(Exp_data1.LOCATION,Selected_countries),:);
joinedData3 = innerjoin(Exp_data1,GDP_data1,"Keys","LOCATION");
Final_Data2 = innerjoin(dist_data,joinedData3,"LeftKeys",["iso_o","iso_d"],"RightKeys",...
   ["LOCATION","PARTNER"]);
All_Merge_Data=removevars(Final_Data2,"Year");