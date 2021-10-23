"""

ss.sdat <- input.rdを取得
center <- centerを取得

rc : 隣り合うマスク同士の中心間距離

Ptのmaskを一番近いセンターのmaskに合わせる
PtG : それぞれのPtの重心を出す．

rs : shorter distance from center
rs = sqrt(PtG^2 - center[mask])
rl : longer distance from center

mask 2~5
rl1 =sqrt(PtG^2 - center[mask+1])
rl2 =sqrt(PtG^2 - center[mask-1])
rl = bigger(rl1 or rl2)

rp : prove radius
rp = (rs * rc^2)/(-rl^2+rs^2+rc^2) -rs

"""
