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

[!Note]
測定はWMAP9を仮定している(see [More+ Section II.D.2](https://arxiv.org/pdf/2304.00703))ので、fittingの時にmeasurement correctionを入れる必要あり。

## model 
darkemulator + measurement correction module(hsc y3 likeより)を使う。

これに加えてabundanceにもmeasurement correctionを入れる必要あり。
視線方向と視線と垂直方向の距離を$\pi$, $r_{\rm p}$とする。これはcosmology dependent。観測結果固定で
$$
n_{\rm g}r_{\rm p}{\rm d}r_{\rm p}{\rm d}\pi \propto n_{\rm g}\frac{\chi^2}{E} = {\rm const}
$$
なので宇宙論を変えると
$$
n_{\rm g}(\mathcal{C}_{\rm meas}) 
=
n_{\rm g}(\mathcal{C})
\frac{\chi^2(\mathcal{C})}{\chi^2(\mathcal{C}_{\rm meas})} 
\frac{E(\mathcal{C}_{\rm meas})}{E(\mathcal{C})}
$$