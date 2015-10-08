import re
sequence = "ATGCTATCGATCGATCGATCATGCGATATGCACATCAGAATGCTCGA"
for startIdx in [match.start() for match in re.finditer('ATG', sequence)]:
	print ('ATG starting at idx:', startIdx)