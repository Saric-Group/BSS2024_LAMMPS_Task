#conda create --name BSS24 freud pandas ovito --channel conda-forge   # old

conda create --name BSS24 freud pandas ovito=3.10.6 conda-forge::matplotlib scipy  --strict-channel-priority -c https://conda.ovito.org -c conda-forge
# 'import matplotlib' has to be before 'import ovito'


