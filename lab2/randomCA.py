import random
randomStr = ""
for i in range(0, 1000):
	randomStr += 'C' if random.random() <= 0.6 else 'A'
print (randomStr)