# Surya Saha and Cathy Gresham @ AgBase
# Dec 2020
#######################################


from goatools.base import download_go_basic_obo
obo_fname = download_go_basic_obo()

import sys
from goatools.obo_parser import GODag

godag = GODag("go-basic.obo", optional_attrs={'part_of'})

# shortest path is level and longest is depth and includes the root
# GO:0000001      level-06        depth-07        mitochondrion inheritance [biological_process]
# https://www.ebi.ac.uk/QuickGO/term/GO:0000001
# GO:0000002      level-06        depth-06        mitochondrial genome maintenance [biological_process]
# https://www.ebi.ac.uk/QuickGO/term/GO:0000002

with open('write_hier_all.txt', 'w') as f:
        for rec in godag:
                go_term=godag[rec]
                print(go_term,file=f)
f.close()