#!usr/bin/python
import os,sys
the_cpv_types=['nodcpv']
the_graph_types=['graph']
the_agamma_types = ['','_agamma','_nocdf_agamma']
print 'commands to run:'
for cpv_type in the_cpv_types:
    for graph_type in the_graph_types:

        for agamma_type in the_agamma_types:
            print('./mixingfit input/nodcpv/%s%s_inputs %s >> output/nodcpv/%s_%s_.out'%(cpv_type,agamma_type,graph_type,cpv_type,agamma_type))
                #os.system('ls testing/')
                #os.system('./mixingfit %s_%s_%s %s >> output/%s_%s_%s.out'%(cpv_type,meas_type,agamma_type,graph_type,cpv_type,meas_type,agamma_type))
                #os.system('ls testing/')
                #if not 'alex' in cpv_type:
            print('cp testing/finalplot.pdf output/nodcpv/finalplot_%s_%s_%s.pdf'%(cpv_type,graph_type,agamma_type))

print 'all done!'
