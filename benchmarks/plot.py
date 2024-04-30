import matplotlib.pyplot as plt
import re

f = open('cilk.txt')


def extract_floats(text):
    # Regular expression to find float numbers
    float_pattern = r'[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?'
    floats = re.findall(float_pattern, text)
    return [float(num) for num in floats]

TreeHeight = dict()

TreeHeight['tri'] = dict()
TreeHeight['banded'] = dict()

currDict = TreeHeight['tri']
currLevel = 0

currExp = 'Tree height experiments:'

Task = {'tri': dict(), 'band': dict()}
TaskUntied = {'tri': dict(), 'band': dict()}
Cilk = {'tri': dict(), 'band': dict()}

currMatType = None
nthrd = 0
for line in f:

    if line.strip() == 'Tree height experiments:':
        currExp = 'Tree height experiments:'

    if line.strip() == 'Comparision experiments:':
        currExp = 'Comparision experiments:'

    if currExp == 'Tree height experiments:':
        if line.strip() == 'Tridiagonal:':
            currDict = TreeHeight['tri']
        if line.strip() == 'Banded:':
            currDict = TreeHeight['banded']

        l = line.strip().split(' ')

        if l[0] == 'Cutoff':
            currLevel = int(l[1])
            continue

        numbers = extract_floats(line)

        if not numbers:
            continue

        currDict[currLevel] = currDict.get(currLevel, 0) + numbers[0]

    # currDict = None
    if currExp == 'Comparision experiments:':
        if line.strip() == 'Parallel Task tied:':
            currDict = Task
        
        if line.strip() == 'Parallel task untied:':
            currDict = TaskUntied
        
        if line.strip() == 'Parallel Cilk:':
            currDict = Cilk

        if line.strip() == 'Tridiagonal:':
            currMatType = 'tri'

        if line.strip() == 'Banded:':
            currMatType = 'band'

        if currDict is not None and currMatType is not None:
            exp = currDict[currMatType]

        s = line.strip().split(' ')
        if s[0] == 'N' and s[1] == 'threads' and s[2] == '=':
            nthrd = int(s[3])
            continue

        numbers = extract_floats(line)

        if not numbers:
            continue

        exp[nthrd] = exp.get(nthrd, 0) + numbers[0]
        
print(TreeHeight)
print(Task)
print('OMP untied Tasks version:', TaskUntied)
print('Cilk version:', Cilk)    

# fig, ax = plt.subplots()
# x = list(TreeHeight['banded'].keys())
# y = list(TreeHeight['banded'].values())
# ax.set_xlabel('Tree Height')
# ax.set_ylabel('half Band 5 matrx time taken')
# ax.plot(x, y)
# plt.savefig('treeheight.png')