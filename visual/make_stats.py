#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd
import os.path as osp

data = sys.argv[1]
tag = sys.argv[2]
out_dir = osp.dirname(data)

tbl = pd.read_csv(data,
                  sep = '\t',
                  index_col = 0,
                  header = 0,
                 )

new_lbl = tbl['bio_celltype'].values

name,counts = np.unique(new_lbl,
                        return_counts = True)

stats = pd.DataFrame(counts,
                     index = pd.Index(name),
                     columns = pd.Index(['Number of Cells']),
                     )

opth_sts = osp.join(out_dir,
                    '.'.join([tag,
                             'stats.tsv',
                            ]
                           )
                   )

stats.to_csv(opth_sts,
             sep = '\t',
             header = True,
             index = True,
             index_label = 'Cell Type',
            )

