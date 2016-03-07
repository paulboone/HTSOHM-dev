import numpy as np

def FindMissing(file_name, start, stop):

	names = np.genfromtxt(file_name, usecols=0, dtype=str)

#	print names
	missing = []

	for i in np.arange(start, stop + 1):
		mat = 'MAT-' + str(i)

		if mat not in names:
			missing.append(mat)
			print mat

	print 'total :   ', len(missing)
