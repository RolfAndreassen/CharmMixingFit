#!usr/bin/python
import os,sys
the_cpv_types=['allcpv','allcpv_alex']
the_graph_types=['graph','graph_qop_phi']
the_meas_types = ['no_belle_babar','no_belle_babar_cdf']
the_agamma_types = ['hfag_agamma','lhcb_agamma']
print 'commands to run:'
for cpv_type in the_cpv_types:
    for graph_type in the_graph_types:
        for meas_type in the_meas_types:
            for agamma_type in the_agamma_types:
                print('./mixingfit input/%s_%s_%s %s >> output/%s_%s_%s.out'%(cpv_type,meas_type,agamma_type,graph_type,cpv_type,meas_type,agamma_type))
                #os.system('ls testing/')
                #os.system('./mixingfit %s_%s_%s %s >> output/%s_%s_%s.out'%(cpv_type,meas_type,agamma_type,graph_type,cpv_type,meas_type,agamma_type))
                #os.system('ls testing/')
                if not 'alex' in cpv_type:
                    print('cp testing/finalplot.pdf output/finalplot_%s_%s_%s_%s.pdf'%(cpv_type,meas_type,graph_type,agamma_type))

print 'all done!'
