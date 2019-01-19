# Copyright(c) 2019, Yuanhua Huang
# Licensed under the Apache License 2.0 at
# http://www.apache.org/licenses/LICENSE-2.0

from .version import __version__

import pyximport; pyximport.install()
from .utils.gtf import *
from .utils.fasta import *


