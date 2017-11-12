#!/home/joey/Apps/SysBio2/bin/python
''' Author: Weizhou Xing
    System biology module 3
    Part A
'''
import numpy as np

MATCH = 2
MISMATCH = -1
GAP = -4


class return_object(object):
    def __init__(self):
        pass


def global_alignment(s="", t=""):
    m = len(s) + 1
    n = len(t) + 1
    A = np.zeros([m, n])

    for i in range(m):
        A[i, 0] = i*GAP
    for j in range(n):
        A[0, j] = j*GAP

    for i in range(m-1):
        for j in range(n-1):
            if s[i] == t[j]:
                A[i+1, j+1] = np.max([A[i,j]+MATCH, A[i,j+1]+GAP, A[i+1,j]+GAP])
            else:
                A[i+1, j+1] = np.max([A[i,j]+MISMATCH, A[i,j+1]+GAP, A[i+1,j]+GAP])

    r = return_object()
    r.A = A
    r.score = A[m-1,n-1]
    r.start = (0, 0)
    r.stop = (m-1, n-1)
    r.s = s
    r.t = t
    return r


def semi_global_alignment(s="", t=""):
    m = len(s) + 1
    n = len(t) + 1
    A = np.zeros([m, n])

    for i in range(m-1):
        for j in range(n-1):
            if s[i] == t[j]:
                A[i+1, j+1] = max(A[i,j]+MATCH, A[i,j+1]+GAP, A[i+1,j]+GAP)
            else:
                A[i+1, j+1] = max(A[i,j]+MISMATCH, A[i,j+1]+GAP, A[i+1,j]+GAP)
    r = return_object()
    r.A = A
    r.score = max(max(A[m-1,1:]), max(A[1:,n-1]))
    r.stop = (1,2)
    return max(max(A[m-1,1:]), max(A[1:,n-1]))


def local_alignment(s="", t=""):
    m = len(s) + 1
    n = len(t) + 1
    A = np.zeros([m, n])

    for i in range(m-1):
        for j in range(n-1):
            if s[i] == t[j]:
                A[i+1, j+1] = max(0, A[i,j]+MATCH, A[i,j+1]+GAP, A[i+1,j]+GAP)
            else:
                A[i+1, j+1] = max(0, A[i,j]+MISMATCH, A[i,j+1]+GAP, A[i+1,j]+GAP)
    r = return_object()
    r.A = A
    r.score = max(max(A[m-1,1:]), max(A[1:,n-1]))
    r.start = (1,2)
    r.stop = (1,2)
    return np.max(A[:,:])


def trace_back(A, i, j, s, t, ei=0, ej=0):
    '''
    A: score matrix
    i: current index of str1 (row)
    j: current index of str2 (column)
    s: sequence 1
    t: sequecne 2
    ei: row of end point
    ej: column of end point
    '''
    if i==ei and j==ej:
        # TODO dimension of insertion_pos
        insertion_pos = [[list(),list()]]
        return insertion_pos

    if s[i-1] == t[j-1]:
        cost = MATCH
    else:
        cost = MISMATCH
    b = list()
    b.append(A[i,j-1]+GAP == A[i,j])     # b[0] from ->
    b.append(A[i-1,j-1]+cost == A[i,j])  # b[1] from  \
    b.append(A[i-1,j]+GAP == A[i,j])     # b[2] from  |

    if b[0] and not b[1] and not b[2]:
        insertion_pos = trace_back(A,i,j-1,s,t)
        for trace in insertion_pos:
            trace[0].append(i)

    elif not b[0] and b[1] and not b[2]:
        insertion_pos = trace_back(A,i-1,j-1,s,t)

    elif not b[0] and not b[1] and b[2]:
        insertion_pos = trace_back(A,i-1,j,s,t)
        for trace in insertion_pos:
            trace[1].append(j)

    elif not b[0] and b[1] and b[2]:
        insertion_pos2 = trace_back(A,i-1,j,s,t)
        for trace in insertion_pos2:
            trace[1].append(j)
        insertion_pos = insertion_pos2 + trace_back(A,i-1,j-1,s,t)

    elif b[0] and b[1] and not b[2]:
        insertion_pos0 = trace_back(A,i,j-1,s,t)
        for trace in insertion_pos0:
            trace[0].append(i)
        insertion_pos = insertion_pos0 + trace_back(A,i-1,j-1,s,t)

    elif b[0] and not b[1] and b[2]:
        insertion_pos0 = trace_back(A,i,j-1,s,t)
        for trace in insertion_pos0:
            trace[1].append(j)
        insertion_pos2 = trace_back(A,i-1,j,s,t)
        for trace in insertion_pos2:
            trace[0].append(i)
        insertion_pos = insertion_pos0 + insertion_pos2

    elif b[0] and b[1] and b[2]:
        insertion_pos0 = trace_back(A,i,j-1,s,t)
        for trace in insertion_pos1:
            trace[1].append(j)
        insertion_pos2 = trace_back(A,i-1,j,s,t)
        for trace in insertion_pos2:
            trace[0].append(i)
        insertion_pos = insertion_pos0 + insertion_pos2 +\
                        trace_back(A,i-1,j-1,s,t)

    return insertion_pos


def special_semi_global_alignment(s="", t=""):
    # maximum of 10 gaps for each sequence
    m = len(s) + 1
    n = len(t) + 1
    # 3rd dim: [cost, num of gaps for t, num of gaps for s]
    A = np.zeros([m, n, 3],dtype=np.int16)

    for i in range(m-1):
        for j in range(n-1):
            if np.abs(i-j) > 10:
                #A[i,j,0] = -
                break
            # b_gap_direction = False
            # MATCH
            if s[i] == t[j]:
                cost = MATCH
            else:
                cost = MISMATCH

            if np.argmax([A[i,j,0]+cost, A[i,j+1,0]+GAP, A[i+1,j,0]+GAP]) == 0:
                A[i+1,j+1,0] = A[i,j,0]+cost
                A[i+1,j+1,1] = A[i,j,1]
                A[i+1,j+1,2] = A[i,j,2]
            #### gap for "t" 
            elif np.argmax([A[i,j,0]+cost, A[i,j+1,0]+GAP, A[i+1,j,0]+GAP]) == 1:
                #  if "gap for t" <= 10, then, increase "gap for t"
                if A[i,j+1,1] + 1 <= 10:
                    A[i+1,j+1,0] = A[i,j+1,0]+GAP
                    A[i+1,j+1,1] = A[i,j+1,1]+1
                    A[i+1,j+1,2] = A[i,j+1,2]

                # if "gap for t" > 10 and "s+gap more than match" and "gap for s" <= 10
                # increase "gap for s"
                elif A[i,j,0]+cost < A[i+1,j,0]+GAP and A[i+1,j,2] + 1 <= 10:
                    A[i+1,j+1,0] = A[i+1,j,0]+GAP
                    A[i+1,j+1,1] = A[i+1,j,1]
                    A[i+1,j+1,2] = A[i+1,j,2]+1

                # only match or mismatch direction is ok
                else:
                    A[i+1,j+1,0] = A[i,j,0]+GAP
                    A[i+1,j+1,1] = A[i,j,1]
                    A[i+1,j+1,2] = A[i,j,2]
            #### gap for "s"
            elif np.argmax([A[i,j,0]+cost, A[i,j+1,0]+GAP, A[i+1,j,0]+GAP]) == 2:
                # if "gap for s" <= 10, then, increase "gap for s"
                if A[i+1,j,2] + 1 <= 10:
                    A[i+1,j+1,0] = A[i+1,j,0]+GAP
                    A[i+1,j+1,1] = A[i+1,j,1]
                    A[i+1,j+1,2] = A[i+1,j,2]+1

                # if "gap for s" > 10 and "t+gap more than (mis)match" and "gap for t" <= 10
                # increase "gap for t"
                elif A[i,j,0]+cost < A[i,j+1,0]+GAP and A[i,j+1,1] + 1 <= 10:
                    A[i+1,j+1,0] = A[i,j+1,0]+GAP
                    A[i+1,j+1,1] = A[i,j+1,1]+1
                    A[i+1,j+1,2] = A[i+1,j,2]

                # only match or mismatch direction is ok
                else:
                    A[i+1,j+1,0] = A[i,j,0]+cost
                    A[i+1,j+1,1] = A[i,j,1]
                    A[i+1,j+1,2] = A[i,j,2]
    print A[:,:,0]
    return max(max(A[m-1,1:,0]), max(A[1:,n-1,0]))


def read_fasta(path="ebolasequences-1.fasta"):
    with open(path) as f:
        fasta = f.read()

    sequence_list = fasta.split("\n>")
    for i in range(len(sequence_list)):
        start_index = sequence_list[i].find("\n")
        sequence_list[i] = sequence_list[i][start_index+1:]
        sequence_list[i] = sequence_list[i].replace("\n", "")

    return sequence_list


if __name__ == '__main__':
    # (a)
    str1 = "ACAAGGA"
    str2 = "ACAGG"
    print global_alignment(s=str1, t=str2)
    # (b)
    str3 = "AGCCATTACCAATTAAGG"
    str4 = "CCAATT"
    print semi_global_alignment(s=str3, t=str4)
    # (c)
    str5 = "AGCCTTCCTAGGG"
    str6 = "GCTTCGTTT"
    print local_alignment(s=str5, t=str6)
    # (d)
    [str7,str8] = read_fasta(path="/home/joey/Work/systembiology/PartA/ebolasequences-1.fasta")
    # str7 = "ACAAGTAGCTA"
    # str8 = "GGTAGCTAG"
    print special_semi_global_alignment(s=str7,t=str8)
