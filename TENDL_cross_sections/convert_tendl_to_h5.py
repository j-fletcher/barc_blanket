import openmc

ba137m = openmc.data.IncidentNeutron.from_njoy('n-Ba137m.tendl',njoy_exec='/home/jwfletch/NJOY2016/build/njoy')
ba137m.export_to_hdf5('Ba137_m1.h5')

c14 = openmc.data.IncidentNeutron.from_njoy('n-C014.tendl',njoy_exec='/home/jwfletch/NJOY2016/build/njoy')
c14.export_to_hdf5('C14.h5')

cd113m = openmc.data.IncidentNeutron.from_njoy('n-Cd113m.tendl',njoy_exec='/home/jwfletch/NJOY2016/build/njoy')
cd113m.export_to_hdf5('Cd113_m1.h5')

nb93m = openmc.data.IncidentNeutron.from_njoy('n-Nb093m.tendl',njoy_exec='/home/jwfletch/NJOY2016/build/njoy')
nb93m.export_to_hdf5('Nb93_m1.h5')

ra228 = openmc.data.IncidentNeutron.from_njoy('n-Ra228.tendl',njoy_exec='/home/jwfletch/NJOY2016/build/njoy')
ra228.export_to_hdf5('Ra228.h5')
