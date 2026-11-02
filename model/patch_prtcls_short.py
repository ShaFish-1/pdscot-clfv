# -*- coding: utf-8 -*-
"""短名版 prtcls1.mdl 修片（按名字匹配，覆盖 9 个奇数粒子）。
用法: python3 patch_prtcls_short.py [prtcls1.mdl]"""
import sys

path = sys.argv[1] if len(sys.argv) > 1 else 'prtcls1.mdl'
lines = open(path).readlines()

map_ = {
 '~N1': "N1  |~N1|~N1|9000012   |1    |MN1 |!wN1 |1    |   |(N1)  |(N1)  |\n",
 '~N2': "N2  |~N2|~N2|9000014   |1    |MN1p|!wN1p|1    |   |(N2)  |(N2)  |\n",
 '~N3': "N3  |~N3|~N3|9000016   |1    |MN2 |!wN2 |1    |   |(N3)  |(N3)  |\n",
 '~N4': "N4  |~N4|~N4|9000018   |1    |MN2p|!wN2p|1    |   |(N4)  |(N4)  |\n",
 '~N5': "N5  |~N5|~N5|9000020   |1    |MN3 |!wN3 |1    |   |(N5)  |(N5)  |\n",
 '~N6': "N6  |~N6|~N6|9000022   |1    |MN3p|!wN3p|1    |   |(N6)  |(N6)  |\n",
 '~ep': "eta+|~ep|~em|9000005   |0    |mEp |!wEp |1    |   |(ep)  |(em)  |\n",
 '~eR': "etaR|~eR|~eR|9000006   |0    |mR  |!wR  |1    |   |(eR)  |(eR)  |\n",
 '~eI': "etaI|~eI|~eI|9000007   |0    |mI  |!wI  |1    |   |(eI)  |(eI)  |\n",
}

cnt = 0
for i, l in enumerate(lines):
    parts = l.split('|')
    if len(parts) > 1 and parts[1].strip() in map_:
        lines[i] = map_[parts[1].strip()]
        cnt += 1

open(path, 'w').writelines(lines)
print('patched', cnt, 'lines')
