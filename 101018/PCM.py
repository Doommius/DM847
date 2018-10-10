array = ["cctccacccc", "cctcctcccc", "cctacgcccc", "cctccttgcc", "catcctcccg", "catcctcccg","cctccttgcc", "cctacgcccc",
"cctcctcccc", "cctccacccc","actcatcatc", "cctcctcccc","tatccgcccc" ,"tctcatcctg", "actcatccct" ,"gctcaccctt", "cctcatcctg", "actcctccct"]

dictlist = [dict() for x in range(len(array[1]))]

letters = ["c","t","a","g"]

for index,item in enumerate(dictlist):
    for i in letters:
        item[i] = 0

for index,item in enumerate(array):
    for letternumber, letter in enumerate(item):
            dictlist[letternumber][letter] = dictlist[letternumber][letter] + 1

print()
print()
print("PCM")

matrix = [[0 for x in range(len(array[1])+1)] for x in range(len(letters))]

for row,i in enumerate(letters):
    matrix[row][0] = i
    for index,item in enumerate(dictlist):
        matrix[row][index+1] = item[i]

s = [[str(e) for e in row] for row in matrix]
lens = [max(map(len, col)) for col in zip(*s)]
fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
table = [fmt.format(*row) for row in s]
print ('\n'.join(table))

print()
print()
print("weight matrix")

for index,item in enumerate(dictlist):
    item["c"] = item["c"] *0.2
    item["t"] = item["t"] *0.2
    item["a"] = item["a"] *0.3
    item["g"] = item["g"] *0.3


matrix = [[0 for x in range(len(array[1])+1)] for x in range(len(letters))]

for row,i in enumerate(letters):
    matrix[row][0] = i
    for index,item in enumerate(dictlist):
        matrix[row][index+1] = round(item[i],2)

s = [[str(e) for e in row] for row in matrix]
lens = [max(map(len, col)) for col in zip(*s)]
fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
table = [fmt.format(*row) for row in s]
print ('\n'.join(table))