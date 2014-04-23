
python -m cProfile benchmark.py > cprofile.out
cat cprofile.out | sort -k4,4 -nr | less
     #    1    0.502    0.502  382.458  382.458 benchmark.py:2(<module>)
     #  204    0.005    0.000  354.140    1.736 indexing.py:79(__setitem__)
     #  204   43.675    0.214  354.103    1.736 indexing.py:159(_setitem_with_inde
     #  816    0.036    0.000  310.373    0.380 indexing.py:338(setter)
     #  816   42.706    0.052  217.887    0.267 frame.py:1874(__setitem__)
     # 1632  182.750    0.112  182.750    0.112 {method 'copy' of 'numpy.ndarray'