from scipy import stats
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import csv


def normaltest_GLU50unbound(i):

	# write a csv file
	csvout = csv.writer(open("GLU50_dataset"+str(i)+".csv","w"), delimiter=',',quoting=csv.QUOTE_ALL)

	nonRNAbound = ['A','B','C','D','E','F','G','H','I','J','K']

	# for the damage set i GLU50 OEs unbound ring atoms
	GLU50OEatm_densities = []
	for OEnum in [1,2]:
		for chain in nonRNAbound:
			filename = 'TRAPRNAdamage%s_GLU50DENSITYFILES/GLU50OE%s%s.txt' %(str(i),str(OEnum),str(chain)) 
			fileopen = open(filename,'r')
			GLU50OEatm_density = [float(line.split()[0]) for line in fileopen.readlines()]
			fileopen.close()
			GLU50OEatm_densities.append(np.array(GLU50OEatm_density))
			csvout.writerow(GLU50OEatm_density)

	GLU50OEatm_densities.append(np.array(GLU50OEatm_density)+5)
	csvout.writerow(list(np.array(GLU50OEatm_density)+5))

	ANOVA_F,ANOVA_P = stats.f_oneway(GLU50OEatm_densities[0],GLU50OEatm_densities[1],GLU50OEatm_densities[2],
									 GLU50OEatm_densities[3],GLU50OEatm_densities[4],GLU50OEatm_densities[5],
									 GLU50OEatm_densities[6],GLU50OEatm_densities[7],GLU50OEatm_densities[8],
									 GLU50OEatm_densities[9],GLU50OEatm_densities[10],GLU50OEatm_densities[11],
									 GLU50OEatm_densities[12],GLU50OEatm_densities[13],GLU50OEatm_densities[14],
									 GLU50OEatm_densities[15],GLU50OEatm_densities[16],GLU50OEatm_densities[17],
									 GLU50OEatm_densities[18],GLU50OEatm_densities[19],GLU50OEatm_densities[20],
									 GLU50OEatm_densities[21])

	print ANOVA_F
	print ANOVA_P

	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (12, 12)})
	fig = plt.figure()
	for atom in GLU50OEatm_densities:
		sns.kdeplot(atom, shade=True)
	fig.savefig('GLU50kde_test.png',bbox_inches='tight')

