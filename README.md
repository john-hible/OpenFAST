# OpenFAST
Folder with OpenFAST software, IEA 15MW model and 2.6.0 ROSCO controller (SETUP FOR WINDOWS)
(These instructions are for the IEA 15MW model but the same can be said for the NREL 5MW model - which was predominantly used for the project)

For use with 'Shutting Down Wind Turbines Safely' 4YP. 'IEA-15-250-RWT' contains the generic folders and controller files whilst the 'IEA-15-250-RWT-Monopile' folder has the folders needed for that specific variant of wind turbine model
The generic folder contains libdiscon.dll (in ServoData) which is the ROSCO controller file, and 'Monopile' contains the DISCON.IN file (in ServoData) which lets you change the controller settings, 'Monopile' is also where the '.fst' file is (to run the simulation through), and where the output file is generated to (as default anyway)

TO RUN THE SOFTWARE - Right click on empty space within 'monopile' folder – then click on ‘open in terminal’ – then in the terminal write ‘’OpenFAST path’ ‘model_fst_file’ e.g. C:\Users....\openfast_x64.exe IEA-15-240-RWT-Monopile.fst

You can change the wind input in 'InflowFile', initial conditions in 'ElastoDyn' and controller settings in DISCON.IN as well as many other parameters, I'd suggest just looking through all of the files to understand yourself
