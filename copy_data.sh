#!/bin/bash
# This works only on gw

# copy bossdr11 clustering datafrom Hironao's directory
dirname=/work/hironao.miyatake/hsc-cmass/real_data/chains/bossdr11-hsc-fid-b0/data/
cp ${dirname}/wp_sig_z*.dat ./data/
cp ${dirname}/wp_cov_z*_z*.dat ./data/
cp ${dirname}/ng_sig_z*.dat ./data/
cp ${dirname}/ng_cov_10p_z*_z*.dat ./data/