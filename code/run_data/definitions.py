import os
import h5py as h5
from astropy.table import Table


sim_flags_dict = {
      # Fiducial 
    "OldWinds_RemFryer2012": "--mass-loss-prescription BELCZYNSKI2010 --OB-mass-loss VINK2001 --VMS-mass-loss VINK2011 --RSG-mass-loss NJ90 --WR-mass-loss BELCZYNSKI2010 \
 --remnant-mass-prescription FRYER2012 --kick-magnitude-distribution MAXWELLIAN ",
      # Old winds with Muller mandel rem
    "OldWinds_RemMullerMandel": "--mass-loss-prescription BELCZYNSKI2010 --OB-mass-loss VINK2001 --VMS-mass-loss VINK2011 --RSG-mass-loss NJ90 --WR-mass-loss BELCZYNSKI2010 \
 --remnant-mass-prescription MULLERMANDEL --kick-magnitude-distribution MULLERMANDEL ",
      # no BH kick
    "OldWinds_RemFryer2012_noBHkick": "--mass-loss-prescription BELCZYNSKI2010 --OB-mass-loss VINK2001 --VMS-mass-loss VINK2011 --RSG-mass-loss NJ90 --WR-mass-loss BELCZYNSKI2010 \
 --remnant-mass-prescription FRYER2012 --kick-magnitude-distribution MAXWELLIAN \
 --black-hole-kicks ZERO ",
      # no NS or BH kick
    "OldWinds_RemFryer2012_noNSBHkick": "--mass-loss-prescription BELCZYNSKI2010 --OB-mass-loss VINK2001 --VMS-mass-loss VINK2011 --RSG-mass-loss NJ90 --WR-mass-loss BELCZYNSKI2010 \
 --remnant-mass-prescription FRYER2012 --kick-magnitude-distribution MAXWELLIAN \
 --black-hole-kicks ZERO --kick-magnitude-distribution ZERO ",
      # NO main sequence mass loss
    "OldWinds_RemFryer2012_noMSwinds": "--mass-loss-prescription BELCZYNSKI2010 --OB-mass-loss VINK2001 --VMS-mass-loss VINK2011 --RSG-mass-loss NJ90 --WR-mass-loss BELCZYNSKI2010 \
 --remnant-mass-prescription FRYER2012 --kick-magnitude-distribution MAXWELLIAN \
 --OB-mass-loss NONE --VMS-mass-loss VINK2011 ",
      # NO WR winds
    "OldWinds_RemFryer2012_noWRwinds": "--mass-loss-prescription BELCZYNSKI2010 --OB-mass-loss VINK2001 --VMS-mass-loss VINK2011 --RSG-mass-loss NJ90 --WR-mass-loss BELCZYNSKI2010 \
 --remnant-mass-prescription FRYER2012 --kick-magnitude-distribution MAXWELLIAN \
 --wolf-rayet-multiplier 0 ",
      # NO winds
    "OldWinds_RemFryer2012_NOwinds": "--mass-loss-prescription BELCZYNSKI2010 --overall-wind-mass-loss-multiplier 0 \
 --remnant-mass-prescription FRYER2012 --kick-magnitude-distribution MAXWELLIAN ",
      # New winds with Fryer remnants
    "NewWinds_RemFryer2012": "--OB-mass-loss VINK2021 --VMS-mass-loss SABHAHIT2023 --RSG-mass-loss DECIN2023  --WR-mass-loss SANDERVINK2023 \
 --remnant-mass-prescription FRYER2012 --kick-magnitude-distribution MAXWELLIAN ",
      # New winds with Muller Mandel remnants
    "NewWinds_RemMullerMandel": "--OB-mass-loss VINK2021 --VMS-mass-loss SABHAHIT2023 --RSG-mass-loss DECIN2023  --WR-mass-loss SANDERVINK2023 \
 --remnant-mass-prescription MULLERMANDEL --kick-magnitude-distribution MULLERMANDEL ",
}