# $Id: __init__.py 3193 2009-03-28 14:42:20Z oliver $
"""Module to faciliate PMF calculations and analysis for AdK runs.

== Basic usage ==

Loading all files to find starting positions.

1) Load module

   >>> import PMF.angles

2) Create database from existing files (see angles.setup()):

   >>> db = PMF.angles.setup()

3) Query for specific frames:

Find all frames for which the LID angle is 90.3 +/- 1.5 and the NMP
angle is 50.5 +/- 1.5:

   >>> db.qaround(LID=90.3, NMP=50.5, delta=1.5)

Find all frames for 100 <= LID angle <=111 and NMP >= 80:

   >>> db.qrange(LID=(100,111), NMP=(80, None))

4) Reloading the database: either just use setup() again or do

   >>> db = PMF.angles.AngleDB(filename)

5) Adding more data:

   >>> db.insertmany('data/angles/lig_*_angles.pickle')

would add all files that match this pattern.


== WHAM ==

Load the default umbrella-sampled windows only:

   >>> import PMF
   >>> db = PMF.angles.setup(pmfonly=True)

Add information about reference values and force constants for each
window (this is all stored in the table 'trajectories' in the database
and is the equivalent of the 'metadata' file):

  >>> db.populate_trajectories()

Now set up a project for WHAM:

  >>> project = PMF.wham.Project(db)

This uses the default working directory 'pmf/wham'.

Prepare the metadata file:

  >>> project.write_meta_file()

Finally, write all data files that have been matched:

  >>> project.write_data_files()

Then run wham-2d:

 cd pmf/wham
 ../../bin/wham-2d Px=0 30 90 60   Py=0 80 150 70  0.001 300 0 wham2d.meta free.dat

Note:
* Need to increase LINESIZE to 200 in file_read.c.


See docs in the submodules for more; the above is not exhaustive.

"""

__all__ = ["angles","wham", "util", "umbrella", "observables", "analysis"]

import angles,wham,util,umbrella,observables,analysis
