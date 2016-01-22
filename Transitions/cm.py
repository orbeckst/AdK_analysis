# $Id: cm.py 2303 2008-09-20 11:42:59Z oliver $
"""Extensions of matplotlib colormaps."""

import matplotlib
import matplotlib.colors
#from matplotlib.cbook import reversed
LUTSIZE = matplotlib.rcParams['image.lut']

# bottom is black, not blue;
# top is bright red, not dark red (to use contrast)
_jet2_data =   {'red':   ((0.0, 0, 0), (0.35, 0, 0), (0.66, 1, 1), (0.89,0.8, 0.8),
                          (1.0, 1, 1)),
                'green': ((0.0, 0, 0), (0.125,0, 0), (0.375,1, 1), (0.64,1, 1),
                          (0.8, 0.3, 0.3), (1.0,0,0)),
                'blue':  ((0.0, 0, 0),
                          (0.0, 0, 0.5), (0.11, 1, 1), (0.34, 1, 1), (0.65,0, 0),
                          (1.0, 0, 0))
                }

jet2 = matplotlib.colors.LinearSegmentedColormap('jet2',    _jet2_data, LUTSIZE)
