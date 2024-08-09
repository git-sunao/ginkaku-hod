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

宮武さんの[notebook](https://github.com/HSC-S19A-cosmology-analysis/S19A-mock-hod/blob/main/plotSignalModelWithAbundance_allscales.ipynb)を見ると
`gw:/work/hironao.miyatake/hsc-cmass/real_data/chains/bossdr11-hsc-fid-b0/data/`にデータがある。
必要なデータは
- `wp_sig_z*.dat`
- `ng_sig_z*.dat`
- `wp_cov_z*_z*.dat`
- `ng_cov_10p_z*_z*.dat` (<- 同じdirectoryに`ng_cov_50p_z*_z*dat`と`ng_cov_5p_z*_z*.dat`があるがこれは...?)

copy_data.shでコピーできる。

## model 
最新のmodeling codeを使おう(大里さんによると、宮武さんが使っているコードにバグがあったらしいので...(要確認))。
hscs19a3x2pt-likellihoodから必要なmodulesを抽出。

