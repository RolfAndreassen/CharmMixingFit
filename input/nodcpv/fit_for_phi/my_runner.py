#!usr/bin/python
import os,sys
the_cpv_types=['nodcpv']
the_meas_types=['','_nobabar','_nobelle_babar','_nocdf']
the_graph_types=['graph']
the_agamma_types = ['','_hfag_agamma']
print 'commands to run:'
for cpv_type in the_cpv_types:
    for meas_type in the_meas_types:
        for graph_type in the_graph_types:
            for agamma_type in the_agamma_types:
                print('./mixingfit input/nodcpv/fit_for_phi/%s%s%s_inputs %s | tee output/nodcpv/fit_for_phi/%s_%s_%s_phifit.out'%(cpv_type,meas_type,agamma_type,graph_type,cpv_type,meas_type,agamma_type))
#                print('cp testing/finalplot.pdf output/nodcpv/finalplot_%s_%s_%s.pdf'%(cpv_type,meas_type,agamma_type))

print 'all done!'
