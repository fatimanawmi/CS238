# Partial Digest Problem
all_X = []

def checkSymmetry(sol1, sol2, width):
    temp = []
    for i in range(len(sol1)):
        temp.append(abs(sol1[i] - width))

    temp.sort()

    if len(temp) != len(sol2):
        return False

    for i in range(len(sol1)):
        if temp[i] != sol2[i]:
            return False
    return True


def inputFormat(data):
    L = list()
    data = data[1:-1]
    token = data.split(',')
    for i in range(len(token)):
        subToken = token[i].split('(')
        num = int(subToken[0])
        multicity = subToken[1].split(')')[0]
        for j in range(int(multicity)):
            L.append(num)
    return L


def Delete(y, L):
    for i in range(len(y)):
        L.remove(y[i])
    return L

def Add(y, L):
    for i in range(len(y)):
        L.append(y[i])

    return L

def multiset(X):
    deltaX = []
    for i in range(len(X)):
        for j in range(i + 1, len(X)):
            deltaX.append(abs(X[i] - X[j]))
    return deltaX


def multisetDiffInL(y, X, L):
    diff = []
    for i in range(len(X)):
        diff.append(abs(X[i] - y))

    for i in range(len(diff)):
        if diff[i] not in L:
            return False, diff
    return True, diff


def Place(L, X, width):
    if len(L) == 0:
        check = True
        X.sort()
        for sol in all_X:
            if checkSymmetry(X, sol, width):
                check = False
                break
        if check:
            all_X.append(X.copy())
        return
    
    y = max(L)
    flag, diff = multisetDiffInL(y, X, L)

    if flag:
        L = Delete(diff, L)
        X.append(y)
        Place(L, X, width) 
        L = Add(diff, L)
        X.remove(y)

    flag, diff = multisetDiffInL(width-y, X, L)
    if flag:
        L = Delete(diff, L)
        X.append(width-y)
        Place(L, X, width)
        L = Add(diff, L)
        X.remove(width-y)
    return 

def PartialDigest(L):
    width = max(L)
    L.remove(width)
    X = [0, width]
    Place(L, X, width)
  
with open('input.txt', 'r') as file:
    lines = file.readlines()
    for line in lines:
        L = inputFormat(line.strip())
        PartialDigest(L)
        print("Test case:", line.strip())
        if len(all_X) == 0:
            print("No solution")
        else:
            print("The number of distinct solution(s):", len(all_X))
            print("Three distinct solution(s):", all_X[:3])
        all_X.clear()
     





