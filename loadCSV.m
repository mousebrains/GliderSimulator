% Load the CSV file generated by GliderSimulator

myDir = fileparts(mfilename('fullpath'));
fn = fullfile(myDir, 'tpw', 'sim.csv');

a = readtable(fn);
