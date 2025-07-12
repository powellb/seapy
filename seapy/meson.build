python_sources = [
  '__init__.py',
  'cdl_parser.py',
  'debug.py',
  'environ.py',
  'filt.py',
  'lib.py',
  'mapping.py',
  'oa.py',
  'plot.py',
  'qserver.py',
  'river.py',
  'tidal_energy.py',
  'tide.py',
  'timeout.py',
  'constituents.npz',
]

py3.install_sources(
  python_sources,
  subdir: 'seapy'
)

subdir('model')
subdir('roms')
subdir('seawater')
subdir('external')
subdir('hawaii_coast')
