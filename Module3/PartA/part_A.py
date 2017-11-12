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


def myArgmax(l):
    '''
    find all the max index
    i.e. myArgmax([1,2,3,4,5,5,5,3,2,1]) = [4,5,6]
    l: iterable
    '''
    r = [np.argmax(l)]
    for i in range(r,len(l)):
        if l[i] == l[r[0]]:
            r.append(i)


def global_alignment(s="", t=""):
    m = len(s) + 1
    n = len(t) + 1
    A = np.zeros([m, n], dtype=np.int32)

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
    r.alignment = "global"
    r.A = A
    r.s = s
    r.t = t
    r.score = A[-1,-1]
    # r.start = (0, 0)
    r.stop = (m-1, n-1)

    return r


def semi_global_alignment(s="", t=""):
    m = len(s) + 1
    n = len(t) + 1
    A = np.zeros([m, n], dtype=np.int32)

    for i in range(m-1):
        for j in range(n-1):
            if s[i] == t[j]:
                A[i+1, j+1] = max(A[i,j]+MATCH, A[i,j+1]+GAP, A[i+1,j]+GAP)
            else:
                A[i+1, j+1] = max(A[i,j]+MISMATCH, A[i,j+1]+GAP, A[i+1,j]+GAP)
    r = return_object()
    r.alignment = "semi"
    r.A = A
    r.s = s
    r.t = t
    max_r = np.max(A[-1,1:]) # max of last row
    argmax_r = np.argmax(A[-1,1:])+1 # argmax of last row
    max_c = np.max(A[1:,-1]) # max of last column
    argmax_c = np.argmax(A[1:,-1])+1 # argmax of last column
    if max_r > max_c:
        r.score = max_r
        r.stop = (m-1, argmax_r)
    else: 
        r.score = max_c
        r.stop = (argmax_c, n-1)
    # r.start = (0, 0)    
    return r


def local_alignment(s="", t=""):
    m = len(s) + 1
    n = len(t) + 1
    A = np.zeros([m, n], dtype=np.int32)

    for i in range(m-1):
        for j in range(n-1):
            if s[i] == t[j]:
                A[i+1, j+1] = max(0, A[i,j]+MATCH, A[i,j+1]+GAP, A[i+1,j]+GAP)
            else:
                A[i+1, j+1] = max(0, A[i,j]+MISMATCH, A[i,j+1]+GAP, A[i+1,j]+GAP)
    r = return_object()
    r.alignment = "local"
    r.A = A
    r.score = np.max(A[:,:])
    # r.start = (0, 0)
    r.stop = (m-1, n-1)
    r.s = s
    r.t = t
    return r


def trace_back(A, i, j, s, t, ei=0, ej=0, alignment=None):
    '''
    A: score matrix
    i: current index of str1 (row)
    j: current index of str2 (column)
    s: sequence 1
    t: sequecne 2
    ei: row of end point
    ej: column of end point
    return: list of element [insertion pos of s, insertion pos of t, end point]
    '''
    if alignment == "global":
        if i==ei and j==ej: # and!
            insertion_pos = [[list(),list(),(i, j)]]
            return insertion_pos
    elif alignment == "semi":
        if i==ei or j==ej: # or!
            insertion_pos = [[list(),list(), (i, j)]]
            return insertion_pos
    # local alignment
    elif alignment == "local":
        if i==ei and j==ej:
            insertion_pos = [[list(),list(), (i, j)]]
            return insertion_pos
    else:
        print "Function `trace_back` Error: unknow alignment method."
        return

    if s[i-1] == t[j-1]:
        cost = MATCH
    else:
        cost = MISMATCH
    
    insertion_pos = list()
    b = list()
    b.append(A[i,j-1]+GAP == A[i,j])     # b[0] from ->
    b.append(A[i-1,j-1]+cost == A[i,j])  # b[1] from  \
    b.append(A[i-1,j]+GAP == A[i,j])     # b[2] from  |

    if b[0]:
        insertion_pos0 = trace_back(A,i,j-1,s,t,alignment=alignment)
        for trace in insertion_pos0:
            trace[0].append(i)
        insertion_pos = insertion_pos + insertion_pos0

    if b[1]:
        insertion_pos1 = trace_back(A,i-1,j-1,s,t,alignment=alignment)
        insertion_pos = insertion_pos + insertion_pos1

    if b[2]:
        insertion_pos2 = trace_back(A,i-1,j,s,t,alignment=alignment)
        for trace in insertion_pos2:
            trace[1].append(j)
        insertion_pos = insertion_pos + insertion_pos2

    return insertion_pos


def print_sequences(s, t, trace_list, insertSpaces=False):
    '''
    s: sequence 1
    t: sequence 2
    trace_list: insertion_pos list
    insertSpaces: global=false, semi-global=true, local=true
    '''
    for k in range(len(trace_list)):
        trace = trace_list[k]
        
        x = list()
        if trace[0]:
            x.append(s[:trace[0][0]])
            for i in range(len(trace[0])-1):
                x.append(str1[trace[0][i]:trace[0][i+1]])
            x.append(s[trace[0][-1]:])
        else:
            x.append(s)
        sequence1 = "-".join(x)
        
        x = list()
        if trace[1]:
            x.append(t[:trace[1][0]])
            for i in range(len(trace[1])-1):
                x.append(t[trace[1][i]:trace[1][i+1]])
            x.append(t[trace[1][-1]:])
        else:
            x.append(t)
        sequence2 = "-".join(x)

        bar = ""
        bar_length = len(sequence1) if len(sequence1) < len(sequence2) else len(sequence2)
        for i in range(bar_length):
            if sequence1[i] == sequence2[i]:
                bar = bar + "|"
            elif sequence1[i] == "-" or sequence2[i] == "-":
                bar = bar + " "
            else:
                bar = bar + "x"
        print "Trace%d:" % k
        print sequence1+"\n"+ bar + "\n" + sequence2


def special_semi_global_alignment(s="", t=""):
    # maximum of 10 gaps for each sequence
    m = len(s) + 1
    n = len(t) + 1
    # 3rd dim: [cost, num of gaps for t, num of gaps for s]
    # TODO separate A(np.int32) and gap record matrix(np.int8)
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
    # print A[:,:,0]
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


def main():
    print "(a):"
    str1 = "ACAAGGA"
    str2 = "ACAGG"
    # str1 = "ACAAGAGCGTAGA"
    # str2 = "ACAGGFTFCTA"
    ra = global_alignment(s=str1, t=str2)
    print "Score:", ra.score
    insertion_idx_a = trace_back(ra.A, ra.stop[0], ra.stop[1], ra.s, ra.t, alignment=ra.alignment)
    print insertion_idx_a
    print_sequences(ra.s, ra.t, insertion_idx_a, insertSpaces=False)

    print "(b):"
    str3 = "AGCCAATTACCAATTAAGG"
    str4 = "CCAATT"
    rb = semi_global_alignment(s=str3, t=str4)
    print "Score:", rb.score
    # print rb.A
    insertion_idx_b = trace_back(rb.A, rb.stop[0], rb.stop[1], rb.s, rb.t, alignment=rb.alignment)
    print insertion_idx_b
    print_sequences(rb.s, rb.t, insertion_idx_b, insertSpaces=False)

    print "(c):"
    str5 = "AGCCTTCCTAGGG"
    str6 = "GCTTCGTTT"
    rc = local_alignment(s=str5, t=str6)
    print "Score:", rc.score

    # print "(d):"
    # # [str7,str8] = read_fasta(path="/home/joey/Work/systembiology/PartA/ebolasequences-1.fasta")
    # str7 = "ACAAGTAGCTA"
    # str8 = "GGTAGCTAG"
    # print special_semi_global_alignment(s=str7,t=str8)

if __name__ == "__main__":
    main()    
