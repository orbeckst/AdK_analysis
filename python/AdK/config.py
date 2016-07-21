# $Id: config.py.m4 2892 2009-02-08 23:36:01Z oliver $

basedir = '/Users/oliver/Biop/Projects/DIMS/AdK'

# defined positions in NMP,LID angle space
# included directly from data/new_angles/xray_angles.pickle
# (see txt/x-ray_angles.txt)
angles = {'2ori_b': (42.8143, 111.053),
          '2bbw_b': (55.5245, 116.835),
          '1aky_a': (41.8905, 111.924),
          '4ake_b': (72.0281, 142.311),
          '1e4v_b': (43.4444, 105.227),
          '2ak3_b': (51.8663, 165.609),
          '1ake_b': (43.4063, 105.905),
          '1ak2_a': (69.8738, 96.8588),
          '1ake_a': (43.26, 105.843),
          '1zak_b': (44.4289, 116.726),
          '2rh5_a': (66.6401, 141.545),
          '1e4y_b': (43.2864, 115.43),
          '2osb_a': (42.3541, 112.735),
          '2eu8_a': (42.8848, 111.889),
          '2rh5_b': (66.6188, 139.794),
          '2ar7_b': (54.2351, 131.718),
          '3aky_a': (42.2889, 110.256),
          '2ar7_a': (61.8025, 119.266),
          '1dvr_b': (68.7141, 109.298),
          '1zak_a': (44.2819, 121.837),
          '2ak2_a': (67.7287, 100.416),
          '2rh5_c': (67.3607, 142.29),
          '2eck_b': (43.6509, 106.318),
          '1zio_a': (43.5657, 125.091),
          '1ank_b': (43.8986, 106.586),
          '1dvr_a': (66.361, 112.552),
          '2bbw_a': (55.3401, 128.213),
          '1e4v_a': (43.7801, 104.706),
          '1zip_a': (44.2243, 122.767),
          '2osb_b': (41.9574, 110.917),
          '2oo7_a': (42.7474, 112.87),
          '1s3g_a': (43.9408, 107.945),
          '4ake_a': (72.9281, 146.644),
          '1p3j_a': (42.3782, 108.732),
          '1zin_a': (44.3503, 123.293),
          '2eck_a': (43.3888, 106.308),
          '2ori_a': (43.133, 112.676),
          '2aky_a': (41.3324, 108.728),
          '2rgx_a': (43.3607, 122.505),
          '2p3s_a': (43.1418, 109.212),
          '2c9y_a': (63.7843, 88.126),
          '1ank_a': (43.6262, 107.02),
          '2eu8_b': (42.808, 110.969),
          '2ak3_a': (48.7946, 156.329),
          '1e4y_a': (42.7578, 111.76),
         }
angles['1AKE'] = angles['1ake_a']
angles['4AKE'] = angles['4ake_a']
angles['reference_state'] = angles['1AKE']

# use to set umbrella.DCD_DEFAULT_PATTERNS
# specific for our AdK trajectory repository on greenwulf
dcd_default_patterns = [
    '/nfs/greenwulf/scratchg/oliver/Projects/DIMS/AdK_external/animal/pmf/*.dcd',      # umbrella
    '/nfs/greenwulf/scratchg/oliver/Projects/DIMS/AdK_external/deathspud/pmf/*.dcd',   # umbrella
    '/nfs/greenwulf/scratchg/oliver/Projects/DIMS/AdK_external/darthtater/pmf/*.dcd',  # umbrella
    '/nfs/greenwulf/scratchg/oliver/Projects/DIMS/AdK_external/timberwulf/pmf/*.dcd',  # umbrella
    '/nfs/greenwulf/xenon/denniej0/oli/test/pmf/*.dcd',                # umbrella
    '/nfs/greenwulf/xenon/denniej0/oli/test/long/*.dcd',               # equilibrium             
    '/nfs/greenwulf/xenon/denniej0/oli/project/oc[0-9][0-9][0-9].dcd', # DIMS
    '/nfs/greenwulf/xenon/denniej0/oli/project/co[0-9][0-9][0-9].dcd', # DIMS                 
    ]

dcd_testing_patterns = [
    '/home/oliver/Projects/DIMS/AdK_external_fake/animal/pmf/*.dcd',      # umbrella
    '/home/oliver/Projects/DIMS/AdK_external_fake/deathspud/pmf/*.dcd',   # umbrella
    '/home/oliver/Projects/DIMS/AdK_external_fake/darthtater/pmf/*.dcd',  # umbrella                        
    '/home/oliver/Projects/DIMS/AdK_external_fake/greenwulf/test/pmf/*.dcd',                # umbrella
    '/home/oliver/Projects/DIMS/AdK_external_fake/greenwulf/test/long/*.dcd',               # equilibrium             
    '/home/oliver/Projects/DIMS/AdK_external_fake/greenwulf/project/oc[0-9][0-9][0-9].dcd', # DIMS
    '/home/oliver/Projects/DIMS/AdK_external_fake/greenwulf/project/co[0-9][0-9][0-9].dcd', # DIMS
    '/home/oliver/Projects/DIMS/AdK_external/greenwulf/project/oc[0-9][0-9][0-9].dcd', # real for testing    
    ]
