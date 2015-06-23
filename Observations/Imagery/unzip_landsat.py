# This file unzips all directories in the current file path.
#
# LMK, UW, 6/20/2015

import os

files = os.listdir('.')

for file in files:
  if file.endswith('tar.gz'):
    filename = file[0:21]
    if not(os.path.exists(filename)):
      print filename
      os.mkdir(filename)
      os.system('tar xf '+file+' -C '+filename)
