nEvents=300
#To activate the continuos readout, assign a positive value to the rate 
rate=50.e3

#o2sim -n $nEvents -m PIPE ITS TPC >& sim.log
o2sim -n $nEvents -m PIPE ITS TPC FRAME >& sim.log
root -b -q run_digi.C\($rate\) >& dig.log
root -b -q run_clus.C >& clus.log
root -b -q run_trac_its.C\($rate\) >& rec_its.log

#root -b -q runCATracking.C >& rec_tpc.log

root -l -b -q ~/alice/O2/Detectors/TPC/reconstruction/macro/convertClusterToClusterHardware.C'++g("o2clus.root","o2clusTPC_HW.root")'  #convert old clusters to ClusterHardware
root -l -b -q ~/alice/O2/Detectors/TPC/reconstruction/macro/runHardwareClusterDecoderRoot.C'++g("o2clusTPC_HW.root","o2clusTPC_Native.root")'   #convert ClusterHardware to ClusterNative
root -l -b -q ~/alice/O2/Detectors/TPC/reconstruction/macro/runCATrackingClusterNative.C'++g("o2clusTPC_Native.root", "tracksFromNative.root", "cont refX=80 bz=-5.0068597793")'     #Run tracking on ClusterNative type
