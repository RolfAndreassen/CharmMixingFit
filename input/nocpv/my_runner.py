#!usr/bin/python
import os,sys
the_cpv_types=['nocpv']
the_graph_types=['graph']
the_meas_types = ['','_babar','_belle_babar','_belle_babar_cdf']
print 'commands to run:'
for cpv_type in the_cpv_types:
    for graph_type in the_graph_types:
        for meas_type in the_meas_types:
            
            print('./mixingfit input/nocpv/%s%s_inputs %s | tee output/nocpv/%s_%s.out'%(cpv_type,meas_type,graph_type,cpv_type,meas_type))
            print('cp testing/finalplot.pdf output/nocpv/finalplot_%s_%s_%s.pdf'%(cpv_type,meas_type,graph_type))

                  
print 'all done!'
