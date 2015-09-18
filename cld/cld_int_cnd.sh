#Initial conditions clds figures
#control 0 time steps
clds -D 55 -E -e "err" -G 10000 -g 8000 -H 15000 -h 2000 -k 3. -l 103 -M 1 -m 3.e-6 -N 1.e6 -n 0 -o "cdf" -p 20 -r 5 -S 1.2 -s .4 -X -x 50
mv err /cgd/data/zender/data/H.T.final.init.50.err
mv cdf /cgd/data/zender/data/H.T.final.init.50.cdf
#
#control 1 time step
clds -D 55 -E -e "err" -G 10000 -g 8000 -H 15000 -h 2000 -k 3. -l 103 -M 1 -m 3.e-6 -N 1.e6 -n 1 -o "cdf" -p 20 -r 5 -S 1.2 -s .4 -X -x 50
mv err /cgd/data/zender/data/H.T.final.init1.50.err
mv cdf /cgd/data/zender/data/H.T.final.init1.50.cdf
#
#truncated 0 time steps
clds -D 55 -E -e "err" -G 10000 -g 8000 -H 15000 -h 2000 -k 3. -l 103 -M 1 -m 20.e-6 -N 1.e6 -n 0 -o "cdf" -p 20 -r 5 -S 1.2 -s .4 -X -x 35
mv err /cgd/data/zender/data/H.T.final.init.35.err
mv cdf /cgd/data/zender/data/H.T.final.init.35.cdf
#
#truncated 1 time step
clds -D 55 -E -e "err" -G 10000 -g 8000 -H 15000 -h 2000 -k 3. -l 103 -M 1 -m 20.e-6 -N 1.e6 -n 1 -o "cdf" -p 20 -r 5 -S 1.2 -s .4 -X -x 35
mv err /cgd/data/zender/data/H.T.final.init1.35.err
mv cdf /cgd/data/zender/data/H.T.final.init1.35.cdf
#
#spherical 1 time step
clds -D 55 -E -e "err" -G 10000 -g 8000 -H 15000 -h 2000 -k 3. -l 103 -M 1 -m 3.e-6 -N 1.e6 -n 1 -o "cdf" -p 20 -r 5 -S 1.2 -s .4 -x 50
mv err /cgd/data/zender/data/S.T.final.init1.50.err
mv cdf /cgd/data/zender/data/S.T.final.init1.50.cdf
#

