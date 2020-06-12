import sys
import re

'''
'python time_range_out.py text_to_bufr_log_file'
Constructs a set of multiline regexs that match the time field outputs in the
bufr converter log file. But would only match if the time is out of range - not
in the 06Z file time window.
'''

#regex = (
#  'YEAR[ \t]+2018\.0+\n'
#  'MNTH[ \t]+12\.0+\n'
#  'DAYS[ \t]+8\.0+\n'
#  'HOUR[ \t]+3\.0+\n'
#  'MINU[ \t]+27\.0+\n'
#  'SECO[ \t]+23\.0+\n'
#)

regexs = []
regexs.append((
  'YEAR[ \t]+2018\.0+\n'
  'MNTH[ \t]+12\.0+\n'
  'DAYS[ \t]+8\.0+\n'
  'HOUR[ \t]+3\.0+\n'
  'MINU[ \t]+0\.0+\n'
  'SECO[ \t]+0\.0+'
))

regexs.append((
  'YEAR[ \t]+2018\.0+\n'
  'MNTH[ \t]+12\.0+\n'
  'DAYS[ \t]+8\.0+\n'
  'HOUR[ \t]+[0-2]\.0+\n'
  'MINU[ \t]+([+-]?(\d*[.])?\d+)\n'
  'SECO[ \t]+([+-]?(\d*[.])?\d+)'
))

regexs.append((
  'YEAR[ \t]+2018\.0+\n'
  'MNTH[ \t]+12\.0+\n'
  'DAYS[ \t]+8\.0+\n'
  'HOUR[ \t]+[1-2][0-9]\.0+\n'
  'MINU[ \t]+([+-]?(\d*[.])?\d+)\n'
  'SECO[ \t]+([+-]?(\d*[.])?\d+)'
))

regexs.append((
  'YEAR[ \t]+2018\.0+\n'
  'MNTH[ \t]+12\.0+\n'
  'DAYS[ \t]+8\.0+\n'
  'HOUR[ \t]+9\.0+\n'
  'MINU[ \t]+([+-]?(\d*[.])?\d+)\n'
  'SECO[ \t]+([+-]?(\d*[.])?\d+)'
))

lines = open(sys.argv[1], 'r').read()
for regex in regexs:
    matchres = re.findall(regex, lines, re.MULTILINE)
    print('regex:', regex)
    if matchres:
        print('regex found')
        print(matchres[0], '\n')
    else:
        print('Regex was not matched\n')
