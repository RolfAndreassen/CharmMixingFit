#!usr/bin/python
import os,sys
the_cpv_types=['nodcpv']
the_meas_types=['','_nobabar','_nobelle_babar','_nocdf']
the_graph_types=['graph_12','graph_12_xphi','graph_12_yphi']
the_agamma_types = ['','_hfag_agamma']
print 'commands to run:'
for cpv_type in the_cpv_types:
    for meas_type in the_meas_types:
        for graph_type in the_graph_types:
            for agamma_type in the_agamma_types:
                print('./mixingfit input/nodcpv/fit_for_12/%s%s%s_inputs %s | tee output/nodcpv/fit_for_12/%s_%s_%s_12_fit.out'%(cpv_type,meas_type,agamma_type,graph_type,cpv_type,meas_type,agamma_type))
                print('cp testing/finalplot.pdf output/nodcpv/finalplot_%s_%s_%s_%s.pdf'%(cpv_type,meas_type,agamma_type,graph_type))

print 'all done!'
