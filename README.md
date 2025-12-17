# Preparation

ROOTをインストール

# Momentum vs Rate

1. Plot.pyの中のCSVファイルの名前を変更
2. python Plot.py

# ビームプロファイルのプロット

1. インプットのrootファイルはKEKCCにあるので、ログイン
2. ROOTの環境設定
3. z軸をHzにしたい場合は、T0のレートを計算し、SFの値を変更する
4. キャリブレーションを更新しておく(Jet chamberのtime->position, LGのADC->momentum)
5. "root -l -b -q tracking.C'("RUNNUMBER")'"
6. output.rootとPDFファイルができる

# Z軸をMomentumにしたビームプロファイル

キャリブレーションを正しく行う、Z軸のビン幅の調整を正しくすることを行わないと、いいプロットはできない
場合によっては、LGにビームが当たっていないイベントが多い可能性があるので、あらかじめ、低エネルギーのカットは入れておく
やることは"python PlotLG.py"



