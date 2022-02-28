# Fisher-Yates shuffle
from random import randrange
list = ['purple', 'monkey', 'dish', 'washer']
def shuffle(list):  # mutates input list
    i = len(list)
    while i > 1:
        j = randrange(i)  # 0 <= j <= i
        list[j], list[i] = list[i], list[j]
        i = i - 1
    return

import random
items =['a', 'quick', 'brown', 'fox', 'jumped', 'over', 'the', 'fence']
print (random.sample(items, 2))