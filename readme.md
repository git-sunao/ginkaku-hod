# GINKAKU HOD
GINKAKU simulationのためのHODを求める。

## 目的
やることはS16A/S19Aの時と同じで、BOSSで測定したwp, ngをHOD modelでfittingして、HOD parameterを決める。
S19Aの時は、Takahashi full-sky simulation (WMAP cosmology)を使ったので、HODの測定の時もWMAP cosmologyを仮定した。

今回のGINKAKUではPlanck 2018を仮定しているので、HOD測定の時もconsistentにPlanck 2018を使う。

## data
Surhudのコードで測定したwpを使う。ngは宮武さんが昔測定したものを使う。
銀河サンプルはLOWZ, CMASS1, CMASS2 (HSC Y1/Y3の2x2ptの設定と同じ)で、
これはBOSSのflux-limited sampleからluminosity cutを適用して作成したもの。

## model 
最新のmodeling codeを使おう(大里さんによると、宮武さんが使っているコードにバグがあったらしいので...)。
hscs19a3x2pt-likellihoodから必要な部品を抽出。

