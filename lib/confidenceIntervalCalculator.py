import numpy as np
import scipy as sp
import scipy.stats

def mean_confidence_interval(data, confidence=0.95):
	# 95% confidence interval calculation
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    return h

def mean_confidence_interval2(data,confidence=0.95):
	# my own calculation of 95% confidence interval, using scipy function
	a = 1.0*np.array(data)
	m, se = np.mean(a), scipy.stats.sem(a)
	low,high = scipy.stats.norm.interval(confidence,m,se)
	h1 = high - m
	h2 = m - low

	if h1-h2 != 0.:
		print 'Something went wrong'
		print '{} and {} should be identical distances from mean {}..'.format(h1,h2,m)
	else:
		return h1
